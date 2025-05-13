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
    Results = postProcess(Dati, ADati, Results, Flag, Param);
    
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

