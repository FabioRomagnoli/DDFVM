function [xSol, info] = scheme(xPrev, BCs, AD, Flag, Opt, t, dt)
% handles the coupled or splitting scheme
    if strcmp(Flag.scheme, "coupled")
        fun = @(x) assembler(x,xPrev, BCs,  AD, Flag, t, dt);
        [xSol, info] = solver(fun, xPrev, Opt, Flag, AD);
       
    elseif strcmp(Flag.scheme,"split")
        % check options to see if i can change jacobian on the fly
        Flag.operator = "transport";
        funTrans = @(x) assembler(x, xPrev, BCs, AD, Flag, t, dt);
        [xSol, info] = solver(funTrans, xPrev, Opt, Flag, AD);

        Flag.operator = "reaction";
        funReact = @(x) assembler(x, xSol, BCs, AD, Flag, t, dt);
        [xSol, info] = solver(funReact, xSol, Opt, Flag, AD);
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
