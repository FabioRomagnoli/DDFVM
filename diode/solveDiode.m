    function [solution] = solveDiode(AD, F, O)
    options = set_options(O);
    t = 0;
    dt = AD.dt;
    dtBest = 0;

    if ~strcmp(F.loadSol,"no")
        load(fullfile(".\sim\", F.loadSol));
        xold = saveFile.ASol(:,end);
    else
        xold = AD.x0;       % xold serves as the prev time step and as inital guess
    end

    % indexes for interior and boundary elements of full vector x 
    intIdxs = [2:AD.lr-1, AD.lr+2:2*AD.lr-1, 2*AD.lr+2:3*AD.lr-1];
    bcsIdxs = [1, AD.lr,  AD.lr+1,  2*AD.lr, 2*AD.lr+1,  3*AD.lr];

    solution = xold;    % saves first solution ie x0 
    
    % Make sure that it executes exact times steps of tsave
    for it = 2 : numel(AD.tsave)  
        while t < AD.tsave(it) - 1e-10
            %  resume from  best dt
            if dtBest ~= 0 
                dt = dtBest;
                dtBest = 0;
            end 

            % clean leftovers
            if (t + dt - 1e-10 > AD.tsave(it)) 
                dtBest = dt;
	            dt = AD.tsave(it) - t;
            end
            
            xOldRed = xold(intIdxs);    % reduced vector, solve  w/out BCs
            BCs =  xold(bcsIdxs);       % BCs go in known RHS

            % + dt
            fun = @(x) assembleDiode(x,xOldRed, BCs,  AD, F, t, dt);


            if F.CheckGradients
                [isValid, diffMatrix] = checkJacobianTermByTerm(fun, xOldRed, 1e-3);
                % valid = checkGradients(fun,xOldRed,optimoptions("fmincon",FiniteDifferenceType="central"),"Display","on");
                if ~isValid
                    % disp(diffMatrix)
                    error("Failed check gradient\n"); 
                end
            end
            [xnew, ~] = solver(fun, xOldRed, options, F, AD);

            if F.adaptive
                 % + dt/2
                funhalf = @(x) assembleDiode(x, xOldRed, BCs,  AD, F, t, dt/2);
                [xtemp, info] = solver(funhalf, xOldRed, options, F, AD);

                if isfield(info, 'exitflag') && info.exitflag <= 0
                    error("fsolve didn't solve"); 
                end 

                 % + dt/2
                funhalf = @(x) assembleDiode(x, xtemp, BCs,  AD, F, t+dt/2, dt/2);
                [xtemp, ~] = solver(funhalf, xtemp, options, F, AD);

                relErr = computeRelErrors(xnew, xtemp, AD);
                fprintf('t = %12.4e \t dt = %12.4e, \t Vt = %12.6e \n', t*AD.tbar, dt*AD.tbar, AD.Vt(t+dt)*AD.Vbar);
                if relErr < AD.tol
                        t = t + dt;             % update t and dt
                        dt = dt*1.2;
                        fprintf("(+) Increasing dt \n\n");
                        xold(1) = AD.Vt(t);     % update boundary conditon 
                        xold(intIdxs) = xnew;   % update vector at prev dt 
                else
                        fprintf("(-) Halving dt \n\n")
	                    dt = dt/2;
                        dtBest = 0;
                end
            else % of adaptive if
                t = t + dt;
                fprintf('t = %g \t Vt = %g\n', t*AD.tbar, AD.Vt(t)*AD.Vbar);
                xold(1) =  AD.Vt(t);
                xold(intIdxs) = xnew;
            end
        end

        % save solution
        fprintf("Saved solution n = %g, at time = %g\n", it, AD.tsave(it)*AD.tbar)
        disp("----------------------------------------------")
        solution(:,it) = xold;
    end
end



function [F,jac]=assembleDiode(x, x0, BCs,  AD, Flags, t, dt)
    % splitting x into its components
    lrr = AD.lr - 2;

    % get boundary
    v_bc = BCs(1:2);
    v_bc(1) = AD.Vt(t+dt);
    n_bc = BCs(3:4);
    p_bc = BCs(5:6);

    % Extract vectors
    v = [v_bc(1); x(1:lrr); v_bc(2)];
    n = [n_bc(1); x(lrr+1:2*lrr); n_bc(2)];
    p = [p_bc(1); x(2*lrr+1:end); p_bc(2)];
    
    % Matrix definitions
    A_full =  ax_laplacian(AD.r,AD.eps);
    M_full =  ax_mass(AD.r, 1);
    An_full = ax_dd(AD.r, v, AD.mu, AD.Vth, -1);
    Ap_full = ax_dd(AD.r, v, AD.mu, AD.Vth, 1); % Così il numero di valenza è corretto
   
    % Get reduced matrix
    A = A_full(2:end-1,2:end-1);
    M = M_full(2:end-1,2:end-1);
    An = An_full(2:end-1,2:end-1);
    Ap = Ap_full(2:end-1,2:end-1);
    
    % getting the first and last elements of the matrices 
    A_bc = A_full(2:end-1,[1 end]);
    An_bc = An_full(2:end-1,[1 end]);
    Ap_bc = Ap_full(2:end-1,[1 end]);
    
    % Generation term
    R = @(n,p) (AD.ni^2 - n.*p)./(AD.tau*(n+p));

    % Full matrix construction
    zeri = zeros(lrr);
    NL = [A M -M; 
        zeri (M+dt*An) zeri; 
        zeri zeri (M+dt*Ap)];
    
    % previous time step (no bounds)
    n0 = x0(lrr+1:2*lrr);
    p0 = x0(2*lrr+1:end);

    % system elements
    bounds = [A_bc*v_bc; dt*An_bc*n_bc; dt*Ap_bc*p_bc];
    rhs = [M*AD.N(2:end-1); M*n0; M*p0];
    gen = [zeros(lrr,1); (-dt*M*R(n(2:end-1),p(2:end-1))); (-dt*M*R(n(2:end-1),p(2:end-1)))];

    % Build system
    F = NL*x + gen + bounds - rhs;

    
    % JACOBIAN 
    if nargout>1
        dR_dn=@(n,p) -(p.^2 + AD.ni^2)./(AD.tau*((n+p).^2));
        dR_dp=@(n,p) -(n.^2 + AD.ni^2)./(AD.tau*((n+p).^2));

        dRn = diag(dt*M*dR_dn(n(2:end-1),p(2:end-1)));  
        dRp = diag(dt*M*dR_dp(n(2:end-1),p(2:end-1))); 

        dAn = zeros(AD.lr, AD.lr);    
        dAp = zeros(AD.lr, AD.lr);    
    
        for i=2:AD.lr-1
            log_pos = 1/log(AD.r(i)/AD.r(i+1));
            log_neg = 1/log(AD.r(i-1)/AD.r(i));
            
            up_neg = (v(i+1) - v(i))/AD.Vth;
            up_pos = (v(i) - v(i+1))/AD.Vth;
            dw_neg = (v(i) - v(i-1))/AD.Vth;
            dw_pos = (v(i-1) - v(i))/AD.Vth;
    
    
            dAn(i,i-1) = log_neg * (n(i)*DB(dw_neg) + n(i-1)*DB(dw_pos));
            dAn(i,i) = log_pos * (-n(i+1)*DB(up_neg) - n(i)*DB(up_pos)) +...
                            log_neg * (-n(i)*DB(dw_neg) - n(i-1)*DB(dw_pos));
            dAn(i,i+1) = log_pos * (n(i+1)*DB(up_neg) + n(i)*DB(up_pos));
    
    
            dAp(i,i-1) = log_neg * (-p(i)*DB(-dw_neg) - p(i-1)*DB(-dw_pos));
            dAp(i,i) = log_pos * (p(i+1)*DB(-up_neg) + p(i)*DB(-up_pos)) +...
                            log_neg * (p(i)*DB(-dw_neg) + p(i-1)*DB(-dw_pos));
            dAp(i,i+1) = log_pos * (-p(i+1) * DB(-up_neg) - p(i)*DB(-up_pos));
        end
        
        J11 = A;
        J12 = M;
        J13 = -M;
        J21 = dt*AD.mu*dAn(2:end-1, 2:end-1);
        J22 = M + dt*An - dRn;
        J23 = -dRp;
        J31 = dt*AD.mu*dAp(2:end-1, 2:end-1); 
        J32 = -dRn;
        J33 = M + dt*Ap - dRp;

        jac = [J11, J12, J13;
               J21, J22, J23;
               J31, J32, J33];

    end
    
end

function [xnew, info] = solver(fun, x0, options,F,AD)
    info = struct();
    switch lower(F.method)
        case 'fsolve'
            [xnew, fval, exitflag, output] = fsolve(fun, x0, options);
            info.fval = fval;
            info.exitflag = exitflag;
            info.output = output;

        case 'newton'
            [xnew, it] = newtonsys(fun, x0, AD.Nmaxit, AD.Ntoll, AD.Nverbose);
            info.iterations = it;

        otherwise
            error('Unknown method: %s', method);
    end
end

function relErr = computeRelErrors(xnew, xtemp, AD)
    % helper local function to make code more readable
    % compute the realative errors separately

    relErrV = norm(xnew(AD.vIR) - xtemp(AD.vIR)) / max(norm(xtemp(AD.vIR)), eps);
    relErrN = norm(xnew(AD.nIR) - xtemp(AD.nIR)) / max(norm(xtemp(AD.nIR)), eps);
    relErrP = norm(xnew(AD.pIR) - xtemp(AD.pIR)) / max(norm(xtemp(AD.pIR)), eps);
    relErr = relErrV + relErrN + relErrP;
    fprintf("||X||=%12.4e \t ||V||=%12.4e \t ||N||=%12.4e \t ||P||=%12.4e\n",  relErr, relErrV,relErrN,relErrP);
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
    fprintf('[checkJacobianTermByTerm] Size of analytic Jacobian: %d x %d\n', length(Janalytic(:,1)),length(Janalytic(1,:)));

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

    % Compute the absolute difference between analytic and FD Jacobians
    diffMatrix = Janalytic - Jfd;
    
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
            % Adjust the field widths (8 for row/col, 22 for numbers) as needed.
            fprintf('%-8d %-8d %-22.6e %-22.6e %-22.6e %-22.6e\n', r, c, anaVal, fdVal, absDiff, relDiff);
        end
        fprintf('Consider lowering tol or re-checking your derivative implementation.\n');
    end
end
