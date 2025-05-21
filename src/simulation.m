function [solution] = simulation(D, AD, Flag, Opt)
    dt = AD.dt;

    file.Dati = D;
    file.ADati = AD;
    file.Flag = Flag;
    
    % Load solution as starting point 
    if ~strcmp(Flag.loadSol,"no")
        loadedFile = load(fullfile(".\sim\", Flag.loadSol)).file;
        loadKf = length(loadedFile.Res.ASol(1,:));
        % In case of resuming simulation handles the differences
        if strcmp(Flag.loadSol,"checkpoint")  ... % is named checkpoint
                && loadKf ~= loadedFile.Dati.K ... % it hasn't finished 
                && loadedFile.Dati.K  == D.K && loadedFile.Dati.T == D.T % same simulation  time and K

            D.K = D.K - loadKf;
            AD.K = AD.K - loadKf;

            D.T = D.T - D.tsave(loadKf);
            AD.T = AD.T - AD.tsave(loadKf);

            D.tsave = D.tsave(loadKf:end);
            AD.tsave = AD.tsave(loadKf:end);
            solution = loadedFile.Res.ASol;
            fprintf("\nLoaded from checkpoint at time t = %g\n", D.tsave(1));
        % elseif loadKf  == file.Dati.K + 1   % Already completed just loadin
        %     solution = file.Res.ASol;
        %     fprintf("\nLoaded full Solution: %s\n",Flag.loadSol);
        %     return
        else
            % otherwise just take the last  value of the solution properly
            % adimensionalized
            solution(:,1) = loadedFile.Res.Sol(:,end) ./ AD.xBarVec;
            fprintf("\nLoaded initial point from simulation %s\n", Flag.loadSol);
        end
    else
        solution(:,1) = AD.x0;       % x0 serves as the prev time step and as inital guess
    end
 

    % actual computation of solution
    for it = 2 : numel(AD.tsave)      % Make sure that it executes exact times steps of tsave
        x0 = solution(:,end);
        t = AD.tsave(it-1);
        t1 = AD.tsave(it);
        [x1, dt] = timeSteppingLoop(AD, Flag, Opt, x0, t, t1, dt);
        if dt <= AD.dtMin, fprintf("\n dt lower bound reached. Exiting..."); return; end      % avoids repetitive loop  
        if ~Flag.adaptive, dt = AD.dt; end              % edge case for t - t1 not divisible by dt

        % save solution
        fprintf("\nSaved solution n = %g, at time = %g\n", length(solution(1,:))+1, AD.tsave(it)*AD.tbar)
        disp("-----------------------------------------------------")
        solution(:,end+1) = x1;

        % Saves checkpoint
        file.Res.ASol = solution;
        save(fullfile(".\sim\", "checkpoint"), 'file');

        % to stop execution create a file "STOP_NOW.txt" in DDFVM folder
        % cleanly exits loop
        if exist('STOP_NOW.txt', 'file')
            fprintf('\nInterrupt requested. Saving and exiting...\n');
            delete('STOP_NOW.txt');  % optional: clean up after yourself
            return
        end

    end
end


function [x0, dtBest] = timeSteppingLoop(AD, Flag, Opt, x0, t, t1, dt)
    while t + 1e-8 < t1
        dtBest = dt;
        % clean leftovers before saving step
        if (t + dt > t1)
            dt = t1 - t;
        end

        x0r = x0(AD.intIdxs);        % reduced vector, solve  w/out BCs
        BCs =  x0(AD.bcsIdxs);       % BCs go in known RHS

        if Flag.verbose >  1, fprintf("\n1) dt: "); end         % + dt
        [xNew, info] = scheme(x0r, BCs, AD, Flag, Opt, t, dt);
        if isfield(info, 'exitflag') && info.exitflag == 0 && Flag.verbose, dt= dt*0.75; fprintf("\t Reducing dt to %g \n\n", dt); continue; end

        if Flag.adaptive
            if Flag.verbose > 1, fprintf("2) dt/2: "); end   % + dt/2
            [xTemp, info] = scheme(x0r, BCs, AD, Flag, Opt, t, dt/2);
            if isfield(info, 'exitflag') && info.exitflag < 0, error(info.output.message); end

            if Flag.verbose > 1, fprintf("3) dt/2: "); end   % + dt/2
            [xTemp, info] = scheme(xTemp, BCs, AD, Flag, Opt, t+dt/2, dt/2);
            if isfield(info, 'exitflag') && info.exitflag < 0, error(info.output.message); end
    
            relativeError = computeRelErrors(xNew, xTemp, AD, Flag);
            timeScaling = max(0.5, min(1.3, AD.stepBuffer * (AD.tolError / relativeError)^AD.scalingPower));
            if Flag.verbose > 1, fprintf('\n  t =%12.4e \t  dt =%12.4e, \t Vt = %12.6e \t EndVt = %12.6e \t scale = %g', t*AD.tbar, dt*AD.tbar, AD.Vt(t+dt)*AD.Vbar, AD.EndVt(t+dt)*AD.Vbar,  timeScaling); end
            dtNew = dt * timeScaling;
            if relativeError < AD.tolError
                    t = t + dt;
                    if Flag.verbose > 1, fprintf("\n(+) Step accepted\n\n"); end

                    % Update boundary and interior solution
                    x0(1) = AD.Vt(t);           % Update boundary condition
                    x0(AD.intIdxs) = xNew;         % Update the interior solution
                
                    % Increase time step for next step
                    dt = min(dtNew, AD.dtMax);
            else
                    if dt <= AD.dtMin, return; end      % avoids repetitive loop  
                    if Flag.verbose > 1, fprintf("\n(-) Step rejected\n\n"); end
                    dt = max(dtNew, AD.dtMin);
            end
        else % of adaptive if
            t = t + dt;
            if Flag.verbose > 1, fprintf('\t\t t = %-10g \t Vt = %g \t EndVt = %g \n', t*AD.tbar, AD.Vt(t)*AD.Vbar, AD.EndVt(t)*AD.Vbar); end
            x0(1) =  AD.Vt(t);
            x0(AD.lr) =  AD.EndVt(t);
            x0(AD.intIdxs) = xNew;
        end
    end
end



function relativeError = computeRelErrors(xNew, xTemp, AD, Flag)
    relErrV = norm(xNew(AD.vIR) - xTemp(AD.vIR)) / max(norm(xTemp(AD.vIR)), eps);
    relErrN = norm(xNew(AD.nIR) - xTemp(AD.nIR)) / max(norm(xTemp(AD.nIR)), eps);
    relErrP = norm(xNew(AD.pIR) - xTemp(AD.pIR)) / max(norm(xTemp(AD.pIR)), eps);
    relativeError = relErrV + relErrN + relErrP;
    if Flag.verbose > 1, fprintf("\n||X||=%12.4e \t ||V||=%12.4e \t ||N||=%12.4e \t ||P||=%12.4e",  relativeError, relErrV,relErrN,relErrP); end
end 
