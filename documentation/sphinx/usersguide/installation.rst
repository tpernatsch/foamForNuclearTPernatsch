
.. _installation:

------------
Installation
------------

We recommend building OpenFOAM and foamForNuclear on a Linux system or using Windows Subsystem for Linux (WSL). Among Linux distributions, Ubuntu is generally preferred as it tends to minimize installation issues. The instuctions below refer to the use of Ubuntu. For other operating systems, please refer to the `official OpenFOAM documentation <https://www.openfoam.com/>`_. Once OpenFOAM is installed, foamForNuclear should be compiled as a standard OpenFOAM application. 

Building OpenFOAM on Ubuntu
===========================

The foamForNuclear project is built on the `OpenFOAM <https://www.openfoam.com/>`_
open-source CFD software.

To use the OpenFOAM early-release debian/ubuntu repository, you will need to add
the signing key and the repository location(s) to your system. This typically
only needs to be done once per distribution. The most convenient way is to use
the installer script with one of these commands:

.. code-block:: bash

   # Add the repository
   curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash

   # Update the repository information
   sudo apt-get update

   # Install preferred package. Eg,
   sudo apt-get install openfoam2606-default

Add in your ``~/.bashrc`` the following command to access all the functionalities of OpenFOAM.

.. code-block:: bash

   source /usr/lib/openfoam/openfoam2606/etc/bashrc

If installing from source, follow the `instructions here <https://develop.openfoam.com/Development/openfoam/-/blob/master/doc/Build.md>`_.
Check for the appropriate version of OpenFOAM specified in the readme file of the `GitLab repository
<https://gitlab.com/foamForNuclear/foamForNuclear>`_!

Building foamForNuclear on Ubuntu
=================================

All foamForNuclear source code is hosted on `GitLab <https://gitlab.com/foamForNuclear/foamForNuclear>`_.
If you have `git <https://git-scm.com>`_, a modern C++ compiler, and the right version of OpenFOAM installed, you can
download and install foamForNuclear by entering the following commands in a terminal:

.. code-block:: bash

   # Clone the repo
   git clone --recursive https://gitlab.com/foamForNuclear/foamForNuclear.git

   # Compile the project and build the foamForNuclear Python API
   ./Allwmake -j<N> --api

Building the Python API
=======================

The :mod:`foamForNuclear <pythonapi>` Python API must be installed separately. The easiest way-+
to install it is using `pip <https://pip.pypa.io/en/stable/>`_.
From the root directory of the foamForNuclear repository, run:

.. code-block:: bash

   python3 -m pip install .
