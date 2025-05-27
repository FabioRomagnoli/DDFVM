## Features
This code implements a solver for a Drift diffusion system coupled with the poisson equation.
It works for simulating the charges and holes in a diode, or for ions and electrons in plasmas.
For differentiating between the two systems there is one main which is run with different configs
depending on the case that one wants solved. 
To run the code input the parameters that you want to input a config file (a function that retuns userParams, userFlags,  and fsolve  options),
the ones not specified will be taken from default values specified in the corresponding init file in config defaults folder and are also
outputted at every run. A few example configs are found in the folder.
can be run in command line with:
      [userParam, Dati, ADati, Flag, Results] = main(%NAMEOFCONFIGFILE%);
i.e.  [userParam, Dati, ADati, Flag, Results] = main('diodeCase1');

The solver can run with fsolve or with newton (if the jacobian is implemented). 
The code if specified will use an adaptive time stepping scheme, the parameters of which can be tuned. The output (if verbose)
will be in the form:

1) dt: Solved, EF=2	 2) dt/2: Solved, EF=2	 3) dt/2: Solved, EF=2
||X||=  7.6053e-05 	 ||V||=  5.3665e-11 	 ||N||=  4.5074e-05 	 ||P||=  3.0979e-05
  t =  5.0139e-03 	  dt =  1.4790e-04, 	 Vt = 7.866478e+01 	 scale = 1.56351
(+) Step accepted

The first line rapresents the three step in the time stepping scheme, the full step dt and the two subsequent dt/2, and their fsolve exit flag 
The second line offers the norms for relative error of the full vector (v,n,p) and also the sigular elemnets
The third line gives the current time, the attempted dt, the potential that its solving for (which is a function of time) and the time stepping scaling
The last line says whether the full time step was good enough to be accepted or not.

