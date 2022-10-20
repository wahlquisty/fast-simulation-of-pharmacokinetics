# Fast simulation of pharmacokinetics

This repo generates the results of the paper "Fast simulation of pharmacokinetics", submitted to IFAC World Congress 2023 by authors Ylva Wahlquist, Fredrik Bagge Carlson and Kristian Soltesz.

We demonstrate fast simulation of the three-compartment PK model on a large propofol data set.
The data set contains 1033 patients published in [1]. Simulations in Julia are compared to MATLAB simulations. The fast simulator can be found in the Julia package ``FastPKSim.jl'' (https://github.com/wahlquisty/FastPKSim.jl).

In Julia:
``eleveld_pksim.jl'' implements the fast simulator using ``FastPKSim.jl''.
``eleveld_diffeq.jl'' simulates the data set using DifferentialEquations.jl.
``eleveld_pksimAdBd.jl'' simulates the data set by computing the result directly from $x(k+1) = \Phi x + \Gamma u$.

In Matlab:
``pksim_eleveld.m'' simulates the dataset using ``lsim''.


[1] Eleveld, D.J., Colin, P., Absalom, A.R., Struys, M.M.R.F., 2018. Pharmacokinetic–pharmacodynamic model for propofol for broad application in anaesthesia and sedation. British Journal of Anaesthesia 120, 942–959. https://doi.org/10.1016/j.bja.2018.01.018



