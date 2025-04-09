implemented 
- newton
- piecewise voltage increase
- formatted output
- adaptive time step 
- unified compcurrent (RHS of plasma)
- saving and loading files
- dynamic plotting
- unification of solvePlasma and solveDiode? worth it?
- alpha exponential (possibly with linear combination of grad(phi)
- proper check gradient
- improved adaptive time stepping
- Exit flag output for better control of solve
- logger
- dynamic stopping

missing
- operator splitting for diodes  (for plasma too)?
    - assembler needs to  build the different operators
    - can be done with mix and match and a method for coupling them
- the timestepping scheme can be abstracted to be general enough to be used
  in the splitting too
- set_optiosn can be removed?
- 
- fix jacobian for plasmas and check the diode one
 - 



## Features
This code implements a solver for a Drift diffusion system coupled with the poisson equation. Known as Van Roesenboek system
It works for simulating the charges and holes in a diode, or for ions and electrons in plasmas. For differentiating between the two systems there are two mains.
To run the code input the parameters that you want to modify in the respective main, the ones not specified will be taken from default values specified in the corresponding init file. The solver can run with fsolve or with newton (if the jacobian is implemented). The code if specified will use an adaptive time stepping scheme, the parameters of which can be tuned. The output will be in the form:

1) dt: Solved, EF=2	 2) dt/2: Solved, EF=2	 3) dt/2: Solved, EF=2
||X||=  7.6053e-05 	 ||V||=  5.3665e-11 	 ||N||=  4.5074e-05 	 ||P||=  3.0979e-05
  t =  5.0139e-03 	  dt =  1.4790e-04, 	 Vt = 7.866478e+01 	 scale = 1.56351
(+) Step accepted

The first line rapresents the three step in the time stepping scheme, the full step dt and the two subsequent dt/2, and their fsolve exit flag 
The second line offers the norms for relative error of the full vector (v,n,p) and also the sigular elemnets
The third line gives the current time, the attempted dt, the potential that its solving for (which is a function of time) and the time stepping scaling
The last line says whether the full time step was good enough to be accepted or not.


% Random comments  
The major drawback is that it might converge very
slowly or even fail to converge if the starting guess is too far from the actual solution.
Damping – multiplying the update with a factor less than 1 – is known to increase the
convergence region. Another remedy is parameter embedding where one slowly changes
a parameter, always using the old solution as a new starting guess. We have already
employed this embedding technique in Section 2.6 when gradually increasing the bias
voltage in the solution procedure.