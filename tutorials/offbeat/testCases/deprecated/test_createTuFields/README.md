# OFFBEAT test Case: test_createTuField

## Overview

This case is designed to test the capabilities of the "createTuFields" application. The goal is to generate a 3D mesh representing a section of a 1.5-D rod run with TRANSURANUS, and subsequently map the transuranus (TU) fields onto the OFFBEAT 3D mesh.

- **system/TuOFFBEATdic**
    is the dictionary, created with a fortran routine, from where the TU fields and geometrical information are read by the `createTuFields` utility.
