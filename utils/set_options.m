function opt = set_options(F)
    % Start with default options for fsolve
    opt = optimoptions('fsolve');

    % Display iteration output
    if isfield(F, 'Display')
        opt.Display = F.Display;
    end

    % Use user-supplied gradient
    if isfield(F, 'SpecifyObjectiveGradient') && F.SpecifyObjectiveGradient
        opt.SpecifyObjectiveGradient = true;
    end    

    % Optimization algorithm
    if isfield(F, 'Algorithm') && ischar(F.Algorithm)
        opt.Algorithm = F.Algorithm;
    end

    % Optimization algorithm
    if isfield(F, 'FiniteDifferenceType') && ischar(F.FiniteDifferenceType)
        opt.FiniteDifferenceType = F.FiniteDifferenceType;
    end

    % Max iterations
    if isfield(F, 'MaxIterations') && isnumeric(F.MaxIterations)
        opt.MaxIterations = F.MaxIterations;
    end

    % Max function evaulation
    if isfield(F, 'MaxFunctionEvaluations') && isnumeric(F.MaxFunctionEvaluations)
        opt.MaxFunctionEvaluations = F.MaxFunctionEvaluations;
    end


    % Optimality tolerance
    if isfield(F, 'OptimalityTolerance') && isnumeric(F.OptimalityTolerance)
        opt.OptimalityTolerance = F.OptimalityTolerance;
    end

    % Step tolerance (newer versions)
    if isfield(F, 'StepTolerance') && isnumeric(F.StepTolerance)
        opt.StepTolerance = F.StepTolerance;
    end

    % Tolerance on function value
    if isfield(F, 'FunctionTolerance') && isnumeric(F.FunctionTolerance)
        opt.FunctionTolerance = F.FunctionTolerance;  % FunctionTolerance replaces TolFun
    end

end
