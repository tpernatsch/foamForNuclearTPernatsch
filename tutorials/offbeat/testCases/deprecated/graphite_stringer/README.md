## OFFBEAT Graphite Stringer Tutorial

This folder contains a tutorial for simulating the graphite stringer of a MSR graphite core block using the OFFBEAT code. 
The tutorial is designed to showcase how OFFBEAT can be utilized for simple but relevant thermal analysis cases.

### Graphite Stringer Overview

The graphite stringer represents a portion of the graphite block found in the MSRE (Molten Salt Reactor Experiment) core. 
It has a fairly complex cross-sectional shape which is then extruded along the z-direction. 
While the cross-section corresponds closely to that of the MSRE stringer, the height (1m in this model) is only an approximation.

### OFFBEAT Case Description

The OFFBEAT case consists of a model inspired by the MSRE's graphite geometry, though not an exact replica. As we are interested in the graphite temperature and as we assume that we can neglect at a first approximation the thermal expansion effect on the fluid temperature, only the thermal solver is activated, while the mechanics solver is deactivated.

The thermal behavior is modeled with the `solidConduction` thermal solver, while the volumetric heat source is directly set in the `Q` file using the `fromLatestTime` heatSource class. The volumetric heat source follows half of a cosine shape (just an assumption), making the case half-symmetrical from a thermal standpoint. To ensure significant temperature gradients within the graphite block, the volumetric heat source has been increased by a factor of 2 compared to the low-power MSRE configuration.

### Boundary Conditions

The boundary condition for the temperature field is a simple fixed-value, representing the temperature of the fluid in the MSRE. However, the tutorial encourages experimentation with other complex boundary conditions available in OFFBEAT. For example, exploring time-dependent temperature or time-dependent heat exchange coefficient, or coupling the code with an OpenFOAM fluid dynamic simulation would be interesting options.