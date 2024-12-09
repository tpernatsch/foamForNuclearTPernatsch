# Tutorials

This folder contains a set of tutorials that have been selected to provide a relatively comprehensive outlook on the functionalities of GeN-Foam. In `3D_SmallESFR`, essentially all non-obvious or non-OpenFOAM standard input parameters are commented, or they include an indication on the tutorial to look at for additional comments.

In each tutorial:

- a `README.md` file provides a general description of the tutorials

- the input files (dictionaries) are normally commented to allow
identifying the role of non-obvious parameters

- an `Allrun` bash script is provided that can be used to run the tutorial.
The steps of the `Allrun` script can also give an understanding of the steps
to take to run other simulations.

- in computational-intensive cases, an `Allrun_parallel` bash script is
provided to run the tutorial using multiple cores.

- in complex cases, an `Allclean` script is provided to clean up a case after
running it and before another simulation


## Regression test suite

It is highly recommended to run the `Alltest` or `regressionTest` scripts when developing to make sure that GeN-Foam behaves as expected. This script runs all the tutorials contained in this folder, except `testing` and `toBeUpdated`. It is an extensive process that might require a day of calculation depending on the hardware. 

The results of the test are summarized in the [regressionResults.md](regressionResults.md) file.
