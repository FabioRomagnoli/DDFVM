implemented 
- newton
- piecewise voltage increase
- formatted output
- adaptive time step 
- unified compcurrent (RHS of plasma)
- saving and loading files
- dynamic plotting

missing
- unification of solvePlasma and solveDiode? worth it?
- alpha exponential (possibly with linear combination of grad(phi)
- operator splitting for diodes  (for plasma too)?
- proper check gradient



The major drawback is that it might converge very
slowly or even fail to converge if the starting guess is too far from the actual solution.
Damping – multiplying the update with a factor less than 1 – is known to increase the
convergence region. Another remedy is parameter embedding where one slowly changes
a parameter, always using the old solution as a new starting guess. We have already
employed this embedding technique in Section 2.6 when gradually increasing the bias
voltage in the solution procedure.