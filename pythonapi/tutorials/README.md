# Tutorials

A set of tutorials is distributed to provide a relatively comprehensive outlook
on the functionalities of FFN and its Python API. In each tutorial:

- a `README.md` file provides a general description of the tutorials
- an `Allrun.py` Python script is provided that contains the tutorial and that can be used to run the simulation. The steps of the `Allrun` script can also give an understanding of the steps to take to run other simulations
- in computational-intensive cases, an `Allrun_parallel.py` bash script is provided to run the tutorial using multiple cores
- in complex cases, an `Allclean.py` script is provided to clean up a case after running it and before another simulation


## Run integration tests

Some tutorials can be tested using the `pytest` Python package when developing to make sure that
FFN and its Python API behave as expected. This script runs all the tutorials contained in this
folder. It is an extensive process that might require a day of calculation
depending on the hardware.

```bash
# Run regression test in a working folder called "test"
pytest --basetemp=test

# Run multiple tests in parallel
pytest --basetemp=test -n 4
```
