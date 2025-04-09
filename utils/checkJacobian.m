function [isValid, diffMatrix, relDiffMatrix] = checkJacobian(fun, x, tol)
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
 