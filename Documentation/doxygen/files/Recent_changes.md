# Recent changes in the case folder {#CHANGES}

* On September 18th 2022, the `reactorState` dictionary has been moved from constant to dictionary. We are aware that changes in the input deck should be minimized. However, we find that having reactorState under the uniform folder is much more logical and that it greatly simlifies the use of GeN-Foam. In particular, this allows to have a separate reactorState dictionary for every time step, which in turns allows for a natural restart of time-dependent simulations, notably for the point-kinetics solver.

* In 2024 a very large restructring of GeN-Foam has taken place. This has led to many modifications to different dictionaries. While we are aware that these changes can make the input deck more complicated, the new features offer very large advantages in terms of code fexibility. Examples of these will be presented in future publications.

* Addition of XS variable selection with the `xsVariables` dictionary in `constant/neutroRegion/nuclearData`. This subdict is mandatory even if the XS are not perturbed, it can be assign empty (`xsVariables {};`).