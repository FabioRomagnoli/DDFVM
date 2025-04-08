
function [solution] = solve(AD, Flag, Opt)
    options = set_options(Opt);
    t = 0;
    dt = AD.dt;
    dtBest = 0;

    if ~strcmp(Flag.loadSol,"no")
        load(fullfile(".\sim\", Flag.loadSol));
        fprintf("Loaded Solution %s\n",Flag.loadSol);
        xOld = file.Res.ASol(:,end);
    else
        xOld = AD.x0;       % xOld serves as the prev time step and as inital guess
    end

    % indexes for interior and boundary elements of full vector x 
    intIdxs = [2:AD.lr-1, AD.lr+2:2*AD.lr-1, 2*AD.lr+2:3*AD.lr-1];
    bcsIdxs = [1, AD.lr,  AD.lr+1,  2*AD.lr, 2*AD.lr+1,  3*AD.lr];

    solution = xOld;    % saves first solution ie x0 
    
    % Make sure that it executes exact times steps of tsave
    for it = 2 : numel(AD.tsave)  
        while t < AD.tsave(it) - 1e-12
            if exist('STOP_NOW.txt', 'file')
                fprintf('Interrupt requested. Saving and exiting...\n');
                delete('STOP_NOW.txt');  % optional: clean up after yourself
                return
            end

           %  resume from best dt 
           if dtBest ~= 0 
                dt = dtBest;
                dtBest = 0;
           end 

            % clean leftovers
            if (t + dt - 1e-12> AD.tsave(it))
                dtBest = dt;
	            dt = AD.tsave(it) - t;
            end

            xOldRed = xOld(intIdxs);    % reduced vector, solve  w/out BCs
            BCs =  xOld(bcsIdxs);       % BCs go in known RHS

            % + dt
            fun = @(x) assembler(x,xOldRed, BCs,  AD, Flag, t, dt);

             if Flag.CheckGradients 
                [isValid, diffMatrix] = checkJacobian(fun, xOldRed, 1e-3);
                % isValid = checkGradients(fun,xOldRed,optimoptions("fmincon",FiniteDifferenceType="central"),"Display","on");
                if ~isValid
                    % disp(diffMatrix)
                    error("Failed check gradient\n"); 
                end
            end

            [xNew, info] = solver(fun, xOldRed, options, Flag, AD);
            
            
            if isfield(info, 'exitflag') && info.exitflag <= 0
                fprintf("fsolve didn't solve\n"); 
                return;
            end 

            if Flag.adaptive
                % + dt/2
                funhalf = @(x) assembler(x,xOldRed, BCs,  AD, Flag, t, dt/2);
                [xTemp, info] = solver(funhalf, xOldRed, options, Flag, AD);

                if isfield(info, 'exitflag') && info.exitflag <= 0
                    fprintf("fsolve didn't solve\n"); 
                    return;
                end 

                 % + dt/2
                funhalf = @(x) assembler(x,xTemp, BCs,  AD, Flag, t+dt/2, dt/2);
                [xTemp, ~] = solver(funhalf, xTemp, options, Flag, AD);

                relativeError = computeRelErrors(xNew, xTemp, AD);
                timeScaling = AD.stepBuffer * (AD.tolError / relativeError)^AD.scalingPower;
                fprintf('t = %12.4e \t dt = %12.4e, \t Vt = %12.6e \t scale = %g\n', t*AD.tbar, dt*AD.tbar, AD.Vt(t+dt)*AD.Vbar, timeScaling);
                dtNew = dt * timeScaling;
                if relativeError < AD.tolError
                        t = t + dt;
                        fprintf("(+) Step accepted\n\n");

                        % Update boundary and interior solution
                        xOld(1) = AD.Vt(t);           % Update boundary condition
                        xOld(intIdxs) = xNew;         % Update the interior solution
                    
                        % Increase time step for next step
                        dt = min(dtNew, AD.dtMax);
                else
                        if dt == AD.dtMin, return; end      % avoid repetitive loop  
                        fprintf("(-) Step rejected\n\n");	
                        dt = max(dtNew, AD.dtMin);
                        dtBest = 0;
                end
            else % of adaptive if
                t = t + dt;
                fprintf('t = %g \t Vt = %g\n', t*AD.tbar, AD.Vt(t)*AD.Vbar);
                xOld(1) =  AD.Vt(t);
                xOld(intIdxs) = xNew;
            end
        end

        % save solution
        fprintf("Saved solution n = %g, at time = %g\n", it, AD.tsave(it)*AD.tbar)
        disp("----------------------------------------------")
        solution(:,it) = xOld;
    end
end



function [xNew, info] = solver(fun, x0, options,Flag,AD)
    info = struct();
    switch lower(Flag.method)
        case 'fsolve'
            [xNew, fval, exitflag, output] = fsolve(fun, x0, options);
            info.fval = fval;
            info.exitflag = exitflag;
            info.output = output;
        case 'newton'
            [xNew, it] = newtonsys(fun, x0, AD.Nmaxit, AD.Ntoll, AD.Nverbose);
            info.iterations = it;

        otherwise
            error('Unknown method: %s', method);
    end
end


function relativeError = computeRelErrors(xNew, xTemp, AD)
    relErrV = norm(xNew(AD.vIR) - xTemp(AD.vIR)) / max(norm(xTemp(AD.vIR)), eps);
    relErrN = norm(xNew(AD.nIR) - xTemp(AD.nIR)) / max(norm(xTemp(AD.nIR)), eps);
    relErrP = norm(xNew(AD.pIR) - xTemp(AD.pIR)) / max(norm(xTemp(AD.pIR)), eps);
    relativeError = relErrV + relErrN + relErrP;
    fprintf("||X||=%12.4e \t ||V||=%12.4e \t ||N||=%12.4e \t ||P||=%12.4e\n",  relativeError, relErrV,relErrN,relErrP);
end 


function [isValid, diffMatrix, relDiffMatrix]  = checkJacobian(fun, x, tol)
    [Jfd,err] = jacobianest(fun,x);
    [F0, Janalytic] = fun(x);

    % If the analytic Jacobian is sparse, convert to full
    if issparse(Janalytic) , Janalytic = full(Janalytic);  end
    
    % Print the number of nonzero entries in the analytic Jacobian
    numNonZeroAnalytic = nnz(Janalytic);
    fprintf('[checkJacobianTermByTerm] Nonzero entries in analytic Jacobian: %d\n', numNonZeroAnalytic);
    % Print the number of nonzero entries in the FD Jacobian
    fprintf('[checkJacobianTermByTerm] Nonzero entries in FD Jacobian: %d\n', nnz(Jfd));

    % Define an absolute threshold below which differences are considered machine noise.
    absTol = 1e-10;

    % Compute the absolute difference between analytic and FD Jacobians
    diffMatrix = Janalytic - Jfd;
    
    % Zero out differences that are smaller than absTol
    diffMatrix(abs(diffMatrix) < absTol) = 0;

    % Compute relative differences elementwise.
    % Adding eps in the denominator avoids division by zero.
    relDiffMatrix = abs(diffMatrix) ./ (max(abs(Janalytic), abs(Jfd)) + eps);

    % Find entries where the relative difference exceeds the tolerance
    [rowIdx, colIdx] = find(relDiffMatrix > tol);

    if isempty(rowIdx)
        isValid = true;
        fprintf('[checkJacobianTermByTerm] All Jacobian entries match within relative tol = %.2e.\n', tol);
    else
        isValid = false;
        fprintf('[checkJacobianTermByTerm] Mismatch in %d entries (relative tol = %.2e):\n', numel(rowIdx), tol);
        % Use fixed width formatting for neat columns.
        header = sprintf('%-8s %-8s %-22s %-22s %-22s %-22s\n', 'Row', 'Col', 'Analytic', 'FD', 'Abs Diff', 'Rel Diff');
        fprintf('%s', header);
        for k = 1:numel(rowIdx)
            r = rowIdx(k);
            c = colIdx(k);
            anaVal = Janalytic(r, c);
            fdVal  = Jfd(r, c);
            absDiff = anaVal - fdVal;
            relDiff = relDiffMatrix(r, c);
            fprintf('%-8d %-8d %-22.6e %-22.6e %-22.6e %-22.6e\n', r, c, anaVal, fdVal, absDiff, relDiff);
        end
        fprintf('Consider lowering tol or re-checking your derivative implementation.\n');
    end
end


function [isValid, diffMatrix, relDiffMatrix] = checkJacobianTermByTerm(fun, x, tol)
% CHECKJACOBIANTERMBYTERM
%   Compares each entry of the analytic Jacobian to a central finite-difference
%   approximation and checks the relative differences.
%
%   Syntax:
%       [isValid, diffMatrix, relDiffMatrix] = checkJacobianTermByTerm(fun, x, tol)
%
%   Inputs:
%       fun : function handle that returns [Fval, Jval] for a given x.
%             The function should have the signature:
%                [Fval, Jval] = fun(x)
%       x   : current solution vector.
%       tol : relative tolerance for mismatch between analytic and FD Jacobian.
%
%   Outputs:
%       isValid       : logical flag, true if all entries match within tol.
%       diffMatrix    : matrix of absolute differences (analytic minus FD).
%       relDiffMatrix : matrix of relative differences.
%
%   Note:
%       If the analytic Jacobian is sparse, it is converted to full.
    
    % Evaluate function at x to get baseline residual and analytic Jacobian
    [F0, Janalytic] = fun(x);
    
    % If the analytic Jacobian is sparse, convert to full
    if issparse(Janalytic)
        Janalytic = full(Janalytic);
    end
    
    % Print the number of nonzero entries in the analytic Jacobian
    numNonZeroAnalytic = nnz(Janalytic);
    fprintf('[checkJacobianTermByTerm] Nonzero entries in analytic Jacobian: %d\n', numNonZeroAnalytic);
    
    % Dimensions
    m = numel(F0);   % number of residuals
    n = numel(x);    % dimension of x

    % Preallocate finite-difference Jacobian
    Jfd = zeros(m, n);

    % Loop over each variable to compute central differences
    for j = 1:n
        % Choose step size relative to the magnitude of x(j)
        h = 1e-8 * (1 + abs(x(j)));

        % Evaluate at x + h in the jth direction
        xFwd = x;
        xFwd(j) = x(j) + h;
        [Ffwd, ~] = fun(xFwd);

        % Evaluate at x - h in the jth direction
        xBwd = x;
        xBwd(j) = x(j) - h;
        [Fbwd, ~] = fun(xBwd);

        % Central difference approximation for the jth column
        Jfd(:, j) = (Ffwd - Fbwd) / (2 * h);
    end

    % Print the number of nonzero entries in the FD Jacobian
    fprintf('[checkJacobianTermByTerm] Nonzero entries in FD Jacobian: %d\n', nnz(Jfd));
    
    % Define an absolute threshold below which differences are considered machine noise.
    absTol = 1e-10;

    % Compute the absolute difference between analytic and FD Jacobians
    diffMatrix = Janalytic - Jfd;
    
    % Zero out differences that are smaller than absTol
    diffMatrix(abs(diffMatrix) < absTol) = 0;
    
    % Compute relative differences elementwise.
    % Adding eps in the denominator avoids division by zero.
    relDiffMatrix = abs(diffMatrix) ./ (max(abs(Janalytic), abs(Jfd)) + eps);
    
    % For entries where both analytic and FD values are very small, set relative diff to 0.
    smallEntries = (abs(Janalytic) < absTol) & (abs(Jfd) < absTol);
    relDiffMatrix(smallEntries) = 0;

    % Find entries where the relative difference exceeds the tolerance
    [rowIdx, colIdx] = find(relDiffMatrix > tol);

    if isempty(rowIdx)
        isValid = true;
        fprintf('[checkJacobianTermByTerm] All Jacobian entries match within relative tol = %.2e.\n', tol);
    else
        isValid = false;
        fprintf('[checkJacobianTermByTerm] Mismatch in %d entries (relative tol = %.2e):\n', numel(rowIdx), tol);
        % Use fixed width formatting for neat columns.
        header = sprintf('%-8s %-8s %-22s %-22s %-22s %-22s\n', 'Row', 'Col', 'Analytic', 'FD', 'Abs Diff', 'Rel Diff');
        fprintf('%s', header);
        for k = 1:numel(rowIdx)
            r = rowIdx(k);
            c = colIdx(k);
            anaVal = Janalytic(r, c);
            fdVal  = Jfd(r, c);
            absDiff = anaVal - fdVal;
            relDiff = relDiffMatrix(r, c);
            % Adjust the field widths (8 for row/col, 22 for numbers) as needed.
            fprintf('%-8d %-8d %-22.6e %-22.6e %-22.6e %-22.6e\n', r, c, anaVal, fdVal, absDiff, relDiff);
        end
        fprintf('Consider lowering tol or re-checking your derivative implementation.\n');
    end
end
 