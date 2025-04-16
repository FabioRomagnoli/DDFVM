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
    switch Flag.model
        case "diode"
            [Param, Dati, ADati] = initDiode(userParam, Flag);
        case "plasma"
            [Param, Dati, ADati] = initPlasma(userParam, Flag);
    end

    % Print complete configuration: parameters (user-specified vs. defaults), flags, and options
    printConfiguration(userParam, Param, userFlag, Flag, userOpt, Opt);
        
    input('Press Enter to continue...','s');

    % Run the simulation and time it
    tic;
    Results.ASol = simulation(Dati, ADati, Flag, Opt);
    Results.elapsedTime = toc;
    
    % Postprocessing of the results
    Results = postProcess(Dati, ADati, Results, Flag);
    
    % Plotting configuration can be adjusted in the config or by running
    % the function by itself an  modifying the flags needed
    plotter(Results, Dati, Flag);

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
