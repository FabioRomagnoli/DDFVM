function [xSol, info] = scheme(xPrev, BCs, AD, Flag, Opt, t, dt)
% handles the coupled or splitting scheme
    if strcmp(Flag.scheme, "coupled")
        BCs(1) = AD.Vt(t+dt);
        if strcmp(Flag.model,"diode")
            fun = @(x) assemblerDiode(x, xPrev, BCs,  AD, dt);
        elseif strcmp(Flag.model,"plasma")
            if strcmp(Flag.genterm, 'non-const')
                if strcmp(Flag.alpha,"const")
                    fun = @(x) assemblerPlasmaJGenAlphaConst(x,xPrev, BCs, AD, dt);
                elseif strcmp(Flag.alpha,"exp")
                    fun = @(x) assemblerPlasmaJGenAlphaExp(x,xPrev, BCs,  AD, dt);
                end
            elseif strcmp(Flag.genterm, 'const')
                fun = @(x) assemblerPlasmaConstGen(x,xPrev, BCs,  AD, dt);
            end
        end
        
        % Check gradient 
        if Flag.CheckGradients 
            % [isValid, diffMatrix] = checkJacobian(fun, xPrev, 1e-3);
            isValid = checkGradients(fun,xPrev,optimoptions("fmincon",FiniteDifferenceType="central"),"Display","on");
            if ~isValid
                % disp(diffMatrix)
                error("Failed check gradient\n"); 
            end
        end
        
        [xSol, info] = solver(fun, xPrev, Opt, Flag, AD);
        solverOutput(info);

    elseif strcmp(Flag.scheme,"split")
        Flag.operator = "transport";
        fprintf("\n\t Transport: ")
        funTrans = @(x) assemblerTransport(x, xPrev, BCs, AD, dt);
        [xSol, info] = solver(funTrans, xPrev, Opt, Flag, AD);
        solverOutput(info);


        Flag.operator = "reaction";
        fprintf("\n\t Reaction: ")
        funReact = @(x) assemblerDiodeReaction(x, xSol, BCs, AD, dt);
        [xSol, info] = solver(funReact, xSol, Opt, Flag, AD);
        solverOutput(info);
    end
end



function [xNew, info] = solver(fun, x0, options,Flag,AD)
% decides between fsolve and newton
    info = struct();
    switch lower(Flag.method)
        case 'fsolve'
            [xNew, fval, exitflag, output] = fsolve(fun, x0, options);
            info.fval = fval;
            info.exitflag = exitflag;
            info.output = output;
        case 'newton'
            [xNew, it] = newtonsys(fun, x0, AD.Nmaxit, AD.Ntoll, Flag.Nverbose);
            info.iterations = it;
        otherwise
            error('Unknown method: %s', method);
    end
end



function  solverOutput(info)
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
