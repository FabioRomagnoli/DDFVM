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