# Testing for physics coupling

## Description

This a very simple 2-D case that one can use to test and play around with physics coupling. It consists of a steady state based on diffusion and a transient based on point-kinetics

Each physics is simulated on a 2-D rectangular domain in the x-z plain
- the N (neutronics) domain has one cellZone and is between x=0 and x=0.75
- the TH (thermal-hydraulics) domain has one cellZone and is between x=0.25 and x=0.5
- the TM (thermal-mechanics) domain has 3 cellZones: one from x=0.25 and x=0.5; one from x=0.5 and x=0.75; one from x=0.75 and x=1

|---|---|---|---|
| N   N   N |
    | TH|
    | TM  TM  TM|
  
This combination provides the TM domain with: 
- a zone overlapped to both N and TH, where it can get the temperature from the TH
- one region overlapped only with N where it can get the a powerDenisity from 
  N and calculated its own temperature
- one region not overlapped with anything where it will calculate T without sources

Similarly, the N domain will have: 
- a region not overlapped with anything, where it will use default values of XS
- a region overlapped with both TM and TH, where it will get field from both
- a region overlapped only with TM

In TH, coolant (properties ~ sodium) flows from the bottom of the domain to
the top (along z). 
In N, cross-sections have been arbitrarily chosen so to obtain a critical
system in 2 groups. They are used in an initial steady-state diffusion
calculation to obtain the power shape to be used in the point kinetics
transient calculations. 
The TM feature a fixed boundary on the right side.
In case an initial diffusion calculation is not performed, the point kinetic solver 
will assume a spatially uniform power.
