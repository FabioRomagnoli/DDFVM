function [xSol, info] = scheme(xPrev, BCs, AD, Flag, Opt, t, dt)
% handles the coupled or splitting scheme
    BCs(1) = AD.Vt(t+dt);
    BCs(2) = AD.EndVt(t+dt);
    
    switch lower(Flag.scheme)
        case 'coupled'
            % returns the correct function depending of the chosen problem
            fun = pickAssembler(xPrev, BCs,  AD, t, dt, Flag);
    
            % Check gradient 
            if Flag.CheckGradients 
                checkJacobian(fun, xPrev, 1e-3);
                % checkGradients(fun,xPrev,optimoptions("fmincon",FiniteDifferenceType="central"),"Display","on");
            end
    
            % Compute solution 
            [xSol, info] = solver(fun, xPrev, Opt, Flag, AD);

        case 'split'
            if Flag.verbose > 2; fprintf("\n\t Reaction:\t"); end
            switch lower(Flag.model)
                case 'diode'
                    funReact = @(x) assemblerDiodeReaction(x, xPrev, BCs, AD, dt);
                case 'plasma'
                    Opt.SpecifyObjectiveGradient = false;
                    funReact = @(x) assemblerPlasmaReaction(x, xPrev, BCs, AD, dt);
            end
            [xSol, info] = solver(funReact, xPrev, Opt, Flag, AD);

            if Flag.verbose > 2; fprintf("\n\t Transport:\t "); end

            % % ATTEMPT at Internal substepping (e.g. to achieve smaller dt)
            % t_sub = t;
            % dt_sub = 1e-11/AD.tbar;  % desired smaller timestep
            % n_substeps = round(dt / dt_sub);  % number of substeps
            % for i_sub = 1:n_substeps
            %     BCs(1) = AD.Vt(t_sub+dt_sub);
            % 	BCs(2) = AD.EndVt(t_sub+dt_sub);
            %     % Create function for the smaller timestep
            %     funTrans = @(x) assemblerTransport(x, xSol, BCs, AD, dt_sub, Flag);
            % 
            %     % Call the solver for each substep
            %     Opt.SpecifyObjectiveGradient = true;
            %     [xSol, info] = solver(funTrans, xSol, Opt, Flag, AD);
            %     t_sub = t_sub + dt_sub;
            %     if Flag.verbose > 2
            %         fprintf("Substep %d/%d completed\n", i_sub, n_substeps);
            %     end
            % end

            % Works properly
            funTrans = @(x) assemblerTransport(x, xSol, BCs, AD, dt,Flag);
            Opt.SpecifyObjectiveGradient = true;
            [xSol, info] = solver(funTrans, xSol, Opt, Flag, AD);
            if Flag.verbose > 2; fprintf("\n\n"); end
        otherwise
            error('Unknown scheme: %s', Flag.scheme);
    end
end


function fun = pickAssembler(xPrev, BCs,  AD, t, dt, Flag)
    switch lower(Flag.model)
        case 'diode'
            % Simple diode system (jacobian exists)
            fun = @(x) assemblerDiode(x, xPrev, BCs,  AD, dt);
        case 'plasma'
            switch lower(Flag.genterm)
                case 'const'
                    % constant generation in ionization length, not
                    % physical, mainly for testin
                    fun = @(x) assemblerPlasmaConstGen(x,xPrev, BCs,  AD, dt);
                case 'non-const'
                    switch lower(Flag.alpha)
                        case 'const'
                            % generation with alpha*Jn, alpha is a constant
                            fun = @(x) assemblerPlasmaJGenAlphaConst(x,xPrev, BCs, AD, dt);
                        case 'exp'
                            % generation with beta*exp(-Ei/E)*Jn 
                            fun = @(x) assemblerPlasmaJGenAlphaExp(x,xPrev, BCs,  AD, dt);
                        otherwise
                            error('Unknown alpha: %s', Flag.alpha);
                    end
                otherwise
                    error('Unknown generation term: %s', Flag.genterm);
            end
        otherwise
            error('Unknown model: %s', Flag.model);
    end 
end


function [xNew, info] = solver(fun, x0, options,Flag,AD)
% decides between fsolve and newton
    info = struct();
    switch lower(Flag.method)
        case 'fsolve'
            % fsolve method
            [xNew, fval, exitflag, output] = fsolve(fun, x0, options);
            info.fval = fval;
            info.exitflag = exitflag;
            info.output = output;
        case 'newton'
            [xNew, it] = newtonsys(fun, x0, AD.Nmaxit, AD.Ntoll, Flag.Nverbose);
            info.iterations = it;
        otherwise
            error('Unknown method: %s', Flag.method);
    end

   if Flag.verbose > 2, solverOutput(info);   fprintf("\t"); end

end


function  solverOutput(info)
% can be modified to return whatever output
    if isfield(info, 'exitflag')
        switch info.exitflag
            case 1
                % fprintf("Equation solved. First-order optimality is small.\n");
                fprintf("Solved, EF=%d", info.exitflag)
            case 2
                % fprintf("Equation solved. Change in x smaller than the specified tolerance, or Jacobian at x is undefined.\n");
                fprintf("Solved, EF=%d", info.exitflag)
            case 3
                % fprintf("Equation solved. Change in residual smaller than the specified tolerance.\n");
                fprintf("Solved, EF=%d", info.exitflag)
            case 4
                % fprintf("Equation solved. Magnitude of search direction smaller than specified tolerance.\n");
                fprintf("Solved, EF=%d", info.exitflag)
            case 0
                % fprintf("Number of iterations exceeded or number of function evaluations exceeded.\n");
                fprintf("Failed, EF=%d", info.exitflag)
            case -1
                fprintf("Output function or plot function stopped the algorithm.\n");
                % disp(info.output.message);
            case -2
                fprintf("Equation not solved. The exit message can have more information.\n");
                % disp(info.output.message);
            case -3
                fprintf("Equation not solved. Trust region radius became too small.\n");
                % disp(info.output.message);
            otherwise
                fprintf("Unknown exitflag: %d\n", info.exitflag);
                % disp(info.output.message);
        end
    end
end
