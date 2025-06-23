# 2D ESFR UTOP: UTOP Transient

In this task, a converged steady-state TH solution of the pool of the ESFR-SMART is coupled to a neutronics PK solver to simulate an Unprotected Transient OverPower, by increasing the reactivity linearly from 0 to 400 pcm in 10 seconds. 

The TH solution contains a simplified sub-scale power model which can be used to retrieve the fuel and cladding temperatures. However, this model gives no information about more complex phenomena, such as pores formation, fission gas release etc. On the other hand, the fuel behaviour model does contain this information, but due to the different scales of the two physics, coupling with the other physics is not straight-forward. 

To overcome this limitation, antoher region, channelRegion, is introduced. This consists of a buffer region where no equations are solved, but where fluid properties are mapped and heat transfer coefficient is computed to set up a BC for the fuel solution.

The workflow of the simulation is the following: the hydraulics and PK solution are coupled together conventionally (using the sub-scale fuel model, which, given the speed of the transient, is accurate enough to reproduce correctly the temperatures to calculate the feedback effects). Parallely, a one way coupling is performed between the fluid and the channel, by using a CHTLoop (Conjugate heat transfer), which allows to set up a BC for the fuel which is computed from:
$$
k\cdot \nabla T = h (T_s - T_f)
$$
Moreover, the PK solution is used to update the LHR in the fuel pin, so that both the power increase and the coolant temperature changes effects are accounted for.

---

# Tasks

1. **Complete regionsDict.** In constant/regionsDict:
    - Complete "regionSolvers/Level_0" to correctly set up the physics that need to be solved
    - Complete "regionSolvers/Level_1" to create a CHTLoop
    - Complete "mappings" to couple the physics together
2. **Complete nuclearData.** In constant/neutroRegion/nuclearData:
    - Define the different reactivity coefficients for fuel/cladding temperature and coolant density
    - Create the reactity map to simulate a 400 pcm insertion in 10 seconds, starting from 5 s.
    - (Optional) Add the CR reactivity map to simulate CR driveline expansion
3. **Run case.** Run the simulation with `GeN-Foam | tee log.GeN-Foam` (20 minutes)
