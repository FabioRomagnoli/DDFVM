function [solution] = solve(AD, Flag, Opt)
    dt = AD.dt;

    % Load solution as starting point 
    if ~strcmp(Flag.loadSol,"no")
        load(fullfile(".\sim\", Flag.loadSol));
        fprintf("Loaded Solution %s\n",Flag.loadSol);
        x0 = file.Res.ASol(:,end);
    else
        x0 = AD.x0;       % x0 serves as the prev time step and as inital guess
    end
 

    solution(:,1) = x0;    % saves first solution ie x0 
    % actual computation of solution
    % Make sure that it executes exact times steps of tsave
    for it = 2 : numel(AD.tsave)  
        x0 = solution(:,it-1);
        t = AD.tsave(it-1);
        t1 = AD.tsave(it);
        [x1, dt] = timeSteppingLoop(AD, Flag, Opt, x0, t, t1, dt);
        if ~Flag.adaptive, dt = AD.dt; end

        % save solution
        fprintf("Saved solution n = %g, at time = %g\n", it, AD.tsave(it)*AD.tbar)
        disp("----------------------------------------------")
        solution(:,it) = x1;
    end
end


function [x0, dtBest] = timeSteppingLoop(AD, Flag, Opt, x0, t, t1, dt)

    while t < t1 - 1e-12
        % to stop execution create a file "STOP_NOW.txt" in DDFVM folder
        if exist('STOP_NOW.txt', 'file')
            fprintf('Interrupt requested. Saving and exiting...\n');
            delete('STOP_NOW.txt');  % optional: clean up after yourself
            return
        end

        % clean leftovers before saving step
        if (t + dt - 1e-12> t1)
            dtBest = dt;
            dt = t1 - t;
        end

        x0r = x0(AD.intIdxs);        % reduced vector, solve  w/out BCs
        BCs =  x0(AD.bcsIdxs);       % BCs go in known RHS

        % + dt
        [xNew, info] = scheme(x0r, BCs, AD, Flag, Opt, t, dt);
        if Flag.verbose, solverOutput("1) dt", info); end
        if info.exitflag == 0, dt= dt*0.75; fprintf("\t Reducing dt\n\n"); continue; end

        if Flag.adaptive
            % + dt/2
            [xTemp, info] = scheme(x0r, BCs, AD, Flag, Opt, t, dt/2);
            if Flag.verbose, solverOutput("\t 2) dt/2", info); end
            if info.exitflag <= 0, return; end

             % + dt/2
            [xTemp, info] = scheme(xTemp, BCs, AD, Flag, Opt, t+dt/2, dt/2);
            if Flag.verbose, solverOutput("\t 3) dt/2", info); end
            if info.exitflag <= 0, return; end
    
            relativeError = computeRelErrors(xNew, xTemp, AD, Flag);
            timeScaling = max(0.5, min(1.3, AD.stepBuffer * (AD.tolError / relativeError)^AD.scalingPower));
            if Flag.verbose, fprintf('  t =%12.4e \t  dt =%12.4e, \t Vt = %12.6e \t scale = %g\n', t*AD.tbar, dt*AD.tbar, AD.Vt(t+dt)*AD.Vbar, timeScaling); end
            dtNew = dt * timeScaling;
            if relativeError < AD.tolError
                    t = t + dt;
                    fprintf("(+) Step accepted\n\n");

                    % Update boundary and interior solution
                    x0(1) = AD.Vt(t);           % Update boundary condition
                    x0(AD.intIdxs) = xNew;         % Update the interior solution
                
                    % Increase time step for next step
                    dt = min(dtNew, AD.dtMax);
            else
                    if dt == AD.dtMin, return; end      % avoid repetitive loop  
                    fprintf("(-) Step rejected\n\n");	
                    dt = max(dtNew, AD.dtMin);
            end
        else % of adaptive if
            t = t + dt;
            if Flag.verbose, fprintf('\n t = %g \t Vt = %g\n', t*AD.tbar, AD.Vt(t)*AD.Vbar); end
            x0(1) =  AD.Vt(t);
            x0(AD.intIdxs) = xNew;
        end
    end
end



function relativeError = computeRelErrors(xNew, xTemp, AD, Flag)
    relErrV = norm(xNew(AD.vIR) - xTemp(AD.vIR)) / max(norm(xTemp(AD.vIR)), eps);
    relErrN = norm(xNew(AD.nIR) - xTemp(AD.nIR)) / max(norm(xTemp(AD.nIR)), eps);
    relErrP = norm(xNew(AD.pIR) - xTemp(AD.pIR)) / max(norm(xTemp(AD.pIR)), eps);
    relativeError = relErrV + relErrN + relErrP;
    if Flag.verbose, fprintf("\n||X||=%12.4e \t ||V||=%12.4e \t ||N||=%12.4e \t ||P||=%12.4e\n",  relativeError, relErrV,relErrN,relErrP); end
end 


function  solverOutput(step,info)
    if isfield(info, 'exitflag')
        fprintf("%s: ", sprintf(step));
        switch info.exitflag
            case 1
                % fprintf("Equation solved. First-order optimality is small.\n");
                % fprintf("Equation solved. Exit flag = %d.\n", info.exitflag);
                fprintf("Solved, EF=%d", info.exitflag)
            case 2
                % fprintf("Equation solved. Change in x smaller than the specified tolerance, or Jacobian at x is undefined.\n");
                % fprintf("Equation solved. Exit flag = %d.\n", info.exitflag);
                fprintf("Solved, EF=%d", info.exitflag)
            case 3
                % fprintf("Equation solved. Change in residual smaller than the specified tolerance.\n");
                % fprintf("Equation solved. Exit flag = %d.\n", info.exitflag);
                fprintf("Solved, EF=%d", info.exitflag)
            case 4
                % fprintf("Equation solved. Magnitude of search direction smaller than specified tolerance.\n");
                % fprintf("Equation solved. Exit flag = %d.\n", info.exitflag);
                fprintf("Solved, EF=%d", info.exitflag)
            case 0
                fprintf("Failed, EF=%d", info.exitflag)
                % fprintf("Number of iterations exceeded or number of function evaluations exceeded.\n");
                % disp(info.output.message);
            case -1
                fprintf("Output function or plot function stopped the algorithm.\n");
                disp(info.output.message);
            case -2
                fprintf("Equation not solved. The exit message can have more information.\n");
                disp(info.output.message);
            case -3
                fprintf("Equation not solved. Trust region radius became too small.\n");
                disp(info.output.message);
            otherwise
                % Unrecognized exitflag
                fprintf("Unknown exitflag: %d\n", info.exitflag);
                disp(info.output.message);
        end
    end
end
