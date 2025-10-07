# Developer's Guide

**Work in progress**

All merge request must be precede by an `Alltest` execution:

```bash
./Alltest
# or
./Alltest -j<N>
```

This `Alltest` script runs all the tutorials and compare the results of the simulation to reference values. All the tests should be passed. Depending on the computer, this Alltest can take up to several hours.

If any deviation or error occurs, please make sure to find the source of the bugs in the development recently made on your branch. Any merge request with an incomplete Alltest will be rejected. If this error persist, please contact the GeN-Foam develop team by opening a [GitLab Issue](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/issues) explaining the bug and how to reproduce it.

The `--fast` flag can be added to speed up tests while developing:

```bash
./Alltest --fast
```

With this flag, only the tutorials if an execution of a minute or less are run.
