# Source code {#SOURCE}

The source code is subdivided into 3 main folders:
* main: containing the main  GeN-Foam.C source file (and other files directly employed by it);
* classes: containing the 3 main classed employed to solve for neutronics, thermal-hydraulics and thermal-mechanics, as well as a class for multi-physics controls;
* include: containing specialized versions of some OpenFOAM base files.

As a general rule, most classes have an *include* folder that is used to store all the .H files that are included via "#include" in the definition of the respective class (normally in the .C file).
