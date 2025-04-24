function [userParam, Dati, ADati, Flag, Results] = main(configName)
    % Setup paths 
    addpath(genpath('.\src'));

    % run the conficName function if it exists (.\src\config\) contains
    % examples
    if exist([configName '.m'], 'file')
        [userParam,userFlag,userOpt] = feval(configName);  % This passes arguments properly
    else
        error('Script "%s" does not exist.', scriptName);
    end

    % sets all the flags and fsolve options that were not specified by the 
    % user with the default values. The default scripts are contains in
    % .\init\defaults\
    Flag = setSimulationFlags(userFlag);
    Opt = setFsolveOptions(userOpt);
    


    % Initialize full parameters including defaults via initPlasma
    % Note that initPlasma fills in default values for parameters not specified by userParam.
    if Flag.init
        switch Flag.model
            case "diode"
                [Param, Dati, ADati] = initDiode(userParam, Flag);
            case "plasma"
                [Param, Dati, ADati] = initPlasma(userParam, Flag);
        end

        % Loads solution from file or starts the solution vector
        if ~strcmp(Flag.loadSol,"no")
            [Dati, ADati, Results] = loadFile(Dati, ADati, Flag);
        else 
            Results.ASol(:,1) = ADati.x0;       % x0 serves as the prev time step and as inital guess
        end 

        % Print complete configuration: parameters (user-specified vs. defaults), flags, and options
        printConfiguration(userParam, Param, userFlag, Flag, userOpt, Opt);
        input('Press Enter to continue...','s');
    elseif ~Flag.init && ~strcmp(Flag.loadSol,"no") % Only for loading a solution and outputting
        load(fullfile(".\sim\", Flag.loadSol));
        fprintf("\nLoading solution: %s",Flag.loadSol);
        Dati = file.Dati;
        ADati = file.ADati;
        Results = file.Res;
    end
        
       
    if Flag.solve
        % Run the simulation and time it
        tic;
        Results = simulation(Dati, ADati, Flag, Opt, Results);
        Results.elapsedTime = toc;
    end

    if Flag.postProcess
        % Postprocessing of the results
        Results = postProcess(Dati, ADati, Results, Flag);
    end
    
    % Plotting configuration can be adjusted in the config or by running
    % the function by itself an  modifying the flags needed
    plotter(Results, Dati, Flag);

end


function [D, AD, Res] = loadFile(D, AD, Flag)
    load(fullfile(".\sim\", Flag.loadSol));
    loadKf = length(file.Res.ASol(1,:));
    % In case of resuming simulation handles the differences
    if strcmp(Flag.loadSol,"checkpoint")  ... % is named checkpoint
            && loadKf ~= file.Dati.K ... % it hasn't finished 
            && file.Dati.K  == D.K && file.Dati.T == D.T % same simulation  time and K

        D.K = D.K - loadKf;
        AD.K = AD.K - loadKf;

        D.T = D.T - D.tsave(loadKf);
        AD.T = AD.T - AD.tsave(loadKf);

        D.tsave = D.tsave(loadKf:end);
        AD.tsave = AD.tsave(loadKf:end);
        Res.ASol = file.Res.ASol;
        fprintf("\nLoaded from checkpoint at time t = %g\n", D.tsave(1));
    else
        % otherwise just take the last  value of the solution 
        Res = struct();
        AD.x0 = file.Res.ASol(:,end);
        fprintf("\nLoaded initial point from simulation %s\n", Flag.loadSol);
    end
end



function Flag = setSimulationFlags(Flag)
    def = FlagDefaults();

    % Get field names from the default struct
    defaultFields = fieldnames(def);
    % Loop over each field and add it to s_in if it's not already present
    for i = 1:length(defaultFields)
        if ~isfield(Flag, defaultFields{i})
            Flag.(defaultFields{i}) = def.(defaultFields{i});
        end
    end
end


% Running optimoptions('fsolve') already sets all  default values
function defaultOpt = setFsolveOptions(userOpt)
    defaultOpt = optimoptions('fsolve');
    
    defaultFields = fieldnames(defaultOpt);
    for i = 1:length(defaultFields)
        if isfield(userOpt, defaultFields{i})
            defaultOpt.(defaultFields{i}) = userOpt.(defaultFields{i});
        end
    end

end

function printConfiguration(userParam, fullParam, userFlag, fullFlag, userOpt, Opt)
    % PRINTCONFIGURATION displays the current simulation configuration.
    % It differentiates between the user-specified parameters and the defaults.
    %
    % Inputs:
    %   userParam - structure with parameters provided by the user
    %   fullParam - full parameter structure returned by initPlasma (combining user and default values)
    %   Flag      - structure containing simulation flags
    %   Opt       - structure containing fsolve options

    fprintf('========================================\n');
    fprintf('        Simulation Configuration        \n');
    fprintf('========================================\n');
    
    % Display Parameters
    fprintf('\nParameters:\n');
    fullFields = fieldnames(fullParam);
    for i = 1:length(fullFields)
        fieldName = fullFields{i};
        if isfield(userParam, fieldName)
            specType = 'User specified';
        else
            specType = 'Default';
        end
        
        % Format the field value for display
        value = fullParam.(fieldName);
        if isnumeric(value)
            valueStr = num2str(value,'%g');
        elseif isstring(value) || ischar(value) 
            valueStr = char(string(value));
        else
            valueStr = '[non-displayable]';
        end
        
        fprintf('  %25s: %-25s (%s)\n', fieldName, valueStr, specType);
        % The - left alligns the values. the numbers control the width
        % fprintf('  %-25s: %-25s (%s)\n', fieldName, valueStr, specType);

    end

    % Display Flags
    fprintf('\nFlags:\n');
    flagFields = fieldnames(fullFlag);
    for i = 1:length(flagFields)
        fieldName = flagFields{i};
        
        if isfield(userFlag, fieldName)
            specType = 'User specified';
        else
            specType = 'Default';
        end
        
        % Format the field value for display
        value = fullFlag.(fieldName);
        if isnumeric(value)
            valueStr = num2str(value,'%g');
        elseif isstring(value) || ischar(value) || islogical(value)
            valueStr = char(string(value));
        else
            valueStr = '[non-displayable]';
        end
        
        fprintf('  %25s: %-25s (%s)\n', fieldName, valueStr, specType);
    end
    


    % Display Flags
    fprintf('\nFsolve Options:\n');
    optFields = fieldnames(Opt);
    for i = 1:length(optFields)
        fieldName = optFields{i};
        
        if isfield(userOpt, fieldName)
            specType = 'User specified';
        else
            specType = 'Default';
        end
        
        % Format the field value for display
        value = Opt.(fieldName);
        if isnumeric(value)
            valueStr = num2str(value,'%g');
        elseif isstring(value) || ischar(value) || islogical(value)
            valueStr = char(string(value));
        else
            valueStr = '[non-displayable]';
        end
        
        fprintf('  %25s: %-25s (%s)\n', fieldName, valueStr, specType);
    end
    
    fprintf('========================================\n\n');
end
