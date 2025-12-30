===============
Tips and tricks
===============

- When unsure about the functionalities of a class, take a look at the corresponding C++ documentation.
- When unsure about the meaning of an input parameters, take a look at the tutorials. Many of them are commented. Another option is the famous "banana method". Type in anything (for instance "banana"). In many cases, the solver will tell you that "banana" is not a valid option, and it will suggest valid options.


.. warning::
    ``decomposePar`` and ``reconstructPar`` are not capable to work on FieldFields (i.e., fields of fields, such as fields of multi-node temperatures emplyed in sub-scale power models for porous-medium simulations). This affect the parallel decomposition and reconstruction of objects such as ``T.lumpedNuclearStructure`` and ``T.nuclearFuelPin``. The main consequence is that a user cannot automatically restart a case that has been reconstructed and re-decomposed. The best is simply not to reconstruct/decompose before a restart. If, for whatever reason, that is needed, a user may need to do that manually or through a dedicated python script.
