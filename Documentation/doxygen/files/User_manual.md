# User manual {#USERMAN}

**Work in progress!!**

Most of the capabilities of GeN-Foam are encapsulated into  different C++ classes for *neutronics*, *thermalHydraulics* and *thermoMechanics*.

In addition, a 4th class called *multiPhysicsControl* is employed to streamline the handling of various multi-physics options.

GeN-Foam itself is nothing but a fairly coupling loop that calls various functionalities of the *neutronics*, *thermalHydraulics*, and *thermoMechanics* classes and takes care of transfering coupling fields among them. 

The classes *neutronics*, *thermalHydraulics*, and *thermoMechanics*  translate for the user into the 3 different regions (and meshes). There is no requirement for the three meshes to occupy the same region of space. Consistent mapping of fields is performed and a reference value is given to a field if no correspondence is found in the mesh where its value is being projected from.

The following sections describe the use of the 3 main classes of GeN-Foam, incuding some theory and useful references:
* The [*neutronics* class](@ref NEUTRONICS)
* The [*thermalHydraulics* class](@ref TH)
* The [*thermoMechanics* class](@ref TM)

The following section describes instead the coupling strategy and the general GeN-Foam options.

* [*GeN-Foam*](@ref GF)

N.B.: Sub-solvers and models are normally coded as classes, with usage information included in the header (.H) file. These files can be searched using the seearch function at the top right of the page. 


© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021