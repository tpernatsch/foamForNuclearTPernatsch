
FMU coupling
============

foamForNuclear provides several interface points to communicate with
`Functional Mock-up Units <https://fmi-standard.org/>`_ (FMUs). FMUs are
containers of software and data that are based on a widely employed
communication standard called Functional Mockup Interface (FMI). The FMI is
developed by an industrial consortium led by the Modelica Association.


Compiling
---------

Follow the instructions for

- ECI4FOAM: https://gitlab.com/foam-for-nuclear/ECI4FOAM
- FMU4FOAM: https://gitlab.com/foam-for-nuclear/FMU4FOAM

To include into foamForNuclear you have to export the `LIB_ECI4FOAM` environment variable such as:

.. code :: bash

    # In your .bashrc, must end with "ECI4FOAM"
    export LIB_ECI4FOAM="/home/.../path/to/ECI4FOAM"


Then execute:

.. code :: bash

    ./Allwmake --fmi
    # or
    ./Allwmake --fmi -j<N>


Alternative installation from the original FMU4FOAM repository
--------------------------------------------------------------

To use the FMI coupling interface in foamForNuclear, the user has to install the
`FMU4FOAM <https://github.com/DLR-RY/FMU4FOAM>`_ project developed by the DLR
using the following commands:

First of all, you have to install an old version of **conan**. Nothing works with 2.x

.. code :: bash

    pip install oftest conan==1.58.0


Now clone the repository into the src folder:

.. code :: bash

    cd src
    git clone https://github.com/DLR-RY/FMU4FOAM.git


Now run:

.. code :: bash

    cd FMU4FOAM
    ./build-ECI4FOAM.sh


Then run the ``Allwmake`` in ``ECI4FOAM`` and in the ``FMU4FOAM`` root folder:

.. code :: bash

    cd ECI4FOAM
    ./Allwmake
    # Get to the root folder of FMU4FOAM
    cd ..
    ./Allwmake


Finally, install the Python packages for FMUs:

.. code :: bash

    pip install fmu4foam OMSimulator pythonfmu


To test:

.. code :: bash

    cd examples/heatedRoom
    ./Allrun


To include into foamForNuclear you have to export the ``LIB_ECI4FOAM`` environment
variable such as:

.. code :: bash

    # In your .bashrc, must end with "ECI4FOAM"
    export LIB_ECI4FOAM="/home/.../path/to/ECI4FOAM"

Then the foamForNuclear project can be built as usual.


Features
--------

In this section, a list of FMI inputs/outputs in foamForNuclear is provided with a
link to each class.

Inputs from FMUs
~~~~~~~~~~~~~~~~

.. list-table:: Input from FMUs
    :widths: 50 50
    :header-rows: 1

    * - Feature
      - Description
    * - :ref:`pump <pump>`
      - Momentum source
    * - :ref:`fixedTemperature <fixedTemperature>`
      - Fix-temperature structure
    * - :ref:`fixedPower <fixedPower>`
      - Fix-power structure
    * - :ref:`NusseltWallAndHfromFMU <NusseltWallAndHfromFMUFSHeatTransferCoefficient>`
      - Additional heat transfer coefficient in series
    * - :ref:`pointKinetics <pointKineticNeutronics>`
      - External reactivity in the point-kinetics solver
    * - :ref:`pointKinetics <pointKineticNeutronics>`
      - External neutron source modulation in the point-kinetics solver
    * - :ref:`pointKinetics <pointKineticNeutronics>`
      - Boron reactivity in the point-kinetics solver
    * - :ref:`pointKinetics <pointKineticNeutronics>`
      - Decay power in the point-kinetics solver
    * - :ref:`timeProfile <timeProfile>`
      - Time profile object


Outputs to FMUs
~~~~~~~~~~~~~~~

.. list-table:: Output from FMUs
    :widths: 50 50
    :header-rows: 1

    * - Feature
      - Description
    * - FMU4FOAM
      - External sensor
    * - :ref:`fieldIntegralToFMU <fieldIntegralToFMU>`
      - Field integral over a cellZone


Inputs and outputs from/to an FMU
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Inputs and outputs from/to an FMU
    :widths: 50 50
    :header-rows: 1

    * - Feature
      - Description
    * - :ref:`nuclearFuelFMU <nuclearFuelFMU>`
      - Nuclear fuel structure based on an FMU


Tips
----

If your simulation crashes because of an unexpected error from externalComm
json. Make sure that you have provided all the FMI ports in both codes with the
correct spelling.



.. warning::

   The FMI capabilities in **foamForNuclear** are advanced features that require the coordinated use of multiple software packages, each with critical dependencies. Installation and configuration issues are relatively common. Therefore, this functionality is recommended only for experienced users who are highly familiar with:

   - **Linux** system administration
   - **OpenFOAM**
   - **Python**
   - The **FMI interface**
   - Any FMI-compliant software intended for coupling with **foamForNuclear**
