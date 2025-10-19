# Tutorials

A set of tutorials is distributed to provide a relatively comprehensive outlook on the functionalities of FFN. In each tutorial:

- a `README.md` file provides a general description of the tutorials

- the input files (dictionaries) are normally commented to allow identifying the role of non-obvious parameters

- an `Allrun` bash script is provided that can be used to run the tutorial. The steps of the `Allrun` script can also give an understanding of the steps to take to run other simulations

- in computational-intensive cases, an `Allrun_parallel` bash script is provided to run the tutorial using multiple cores

- in complex cases, an `Allclean` script is provided to clean up a case after running it and before another simulation

N.B.: In [`3D_SmallESFR`](./3D_SmallESFR_NewSolverVerification/), many non-obvious or non-OpenFOAM standard input parameters are commented, or they include an indication on the tutorial to look at for additional comments.

## Cases

The Tutorials folder is divided in four types of tutorials:
- [guidedCases](./guidedCases/): simple cases to teach the user how to use FFN
- [reactorCases](./reactorCases/): complex cases representing partial or complete nuclear reactors
- [featureCases](./featureCases/): tests and validation cases for specific features
- [fmuCases](./fmuCases/): examples of use of the FMI standard


## Regression test suite

The [`Alltest`](./Alltest) script can be used when developing to make sure that FFN behaves as expected. This script runs all the tutorials contained in this folder. It is an extensive process that might require a day of calculation depending on the hardware.

The results of the test are summarized in the [regressionResults.md](regressionResults.md) file. In `--fast` mode, the results are summarized in [regressionResults-fast.md](./regressionResults-fast.md).

```bash
./Alltest

# Parallelize the cases
./Alltest -j<num_jobs>

# To add fmi cases to be run, works in parallel
./Alltest --fmi

# To run only fast cases, works in parallel (< 1 min per case)
./Alltest --fast

# To run with the minimum information printed on the terminal
./Alltest --quiet
```
