.. _installation:

============
Installation
============


-----------------
Building OpenFOAM
-----------------

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
   sudo apt-get install openfoam2506-default


Add in your ``~/.bashrc`` the following command to access all the functionalities
of OpenFOAM.

.. code-block:: bash

   source /usr/lib/openfoam/openfoam2506/etc/bashrc


------------------------
Building Source on Linux
------------------------

All foamForNuclear source code is hosted on `GitLab
<https://gitlab.com/foamForNuclear/foamForNuclear>`_. If you have `git
<https://git-scm.com>`_, and a modern C++ compiler installed, you can
download and install foamForNuclear by entering the following commands in a terminal:

.. code-block:: bash

   # Clone the repo
   git clone --recursive https://gitlab.com/foamForNuclear/foamForNuclear.git

   # Compile the project and build the foamForNuclear Python API
   ./Allwmake -j<N> --api


The :mod:`foamForNuclear <pythonapi>` Python package must be installed separately. The easiest way
to install it is using `pip <https://pip.pypa.io/en/stable/>`_.
From the root directory of the foamForNuclear repository, run:

.. code-block:: bash

   python3 -m pip install .
