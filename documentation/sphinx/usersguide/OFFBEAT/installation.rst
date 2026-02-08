Code Installation
=================

Requirements
------------

OFFBEAT has been developed and tested using Ubuntu and Red Hat Enterprise Linux
distributions and should compile and run without problems on any Linux system
supported by OpenFOAM.

The code supports compilation on both the OpenFOAM Foundation version
(``.org``) and the ESI version (``.com``) of OpenFOAM. Specifically, it is
compatible with:

- `OpenFOAM-9 <https://openfoam.org/release/9>`_ from the OpenFOAM Foundation
- `OpenFOAM-v2512 <https://www.openfoam.com/news/main-news/openfoam-v2512>`_
  from the ESI Group

The OFFBEAT compilation script automatically identifies which version of the
OpenFOAM environment is currently sourced.


Installation
------------

Ensure that OpenFOAM is installed and that the environment is properly
configured for compiling OpenFOAM applications. See the
`OpenFOAM User Guide <https://cfd.direct/openfoam/user-guide/v9-compiling-applications>`_
for more information.

To clone the OFFBEAT repository and install the *master* branch version of the
code, create a folder to contain the source files. Typically, this is inside::

  OpenFOAM/userName-v2512/applications/solvers

(where ``userName`` must be adapted accordingly).

Then execute the following commands:

.. code-block:: bash

   git clone https://gitlab.com/foam-for-nuclear/offbeat.git
   cd offbeat
   git checkout master

Use the ``git checkout`` command to change branch. For instance, to switch to
the ``develop`` branch:

.. code-block:: bash

   git checkout develop


Compilation
-----------

Compile OFFBEAT with either:

.. code-block:: bash

   make

for serial compilation, or:

.. code-block:: bash

   make -j <n>

for parallel compilation, where ``n`` is the number of CPU cores to use.

By default, OFFBEAT is compiled and installed in the OpenFOAM user applications
directory ``$FOAM_USER_APPBIN``. Once compilation is complete, OFFBEAT should be
available from the command line:

.. code-block:: bash

   offbeat


Verification
------------

If OFFBEAT is installed correctly, you should see the typical OpenFOAM splash
screen when launching the solver, for example:

.. code-block:: none

   /*---------------------------------------------------------------------------*\
     =========                 |
     \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
      \\    /   O peration     | Website:  https://openfoam.org
       \\  /    A nd           | Version:  9
        \\/     M anipulation  |
   \*---------------------------------------------------------------------------*/
   Build  : <build_number>
   Exec   : offbeat
   Date   : <date>
   Time   : 11:54:30
   Host   : "pc10115"
   PID    : 319486
   I/O    : uncollated
   Case   : <case_directory>
   nProcs : 1
   sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
   fileModificationChecking : Monitoring run-time modified files using timeStampMaster
                              (fileModificationSkew 10)
   allowSystemOperations : Allowing user-supplied system call operations

   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
   Create time
   ...
