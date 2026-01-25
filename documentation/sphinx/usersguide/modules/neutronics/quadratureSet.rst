.. _userguide_neutronics_quadratureSet:

The *quadratureSet* dictionary
------------------------------

The *quadratureSet* dictionary is used to provide the quadrature set when performing
discrete ordinate calculations. The *quadratureSet* dictionary is found under
*constant/(neutronicsRegionName)/*. It contains the quadrature set for discrete ordinate
calculations. One can find examples of three different quadrature sets in the
tutorial `Godiva_SN
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/constant/neutroRegion/>`_.
S4 and S8 Chebyshev-Legendre quadrature sets can be found in `Godiva_SN
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tools/chebichevLegendreQuadratureSets/>`_.

Here follows the list of entries to set the quadratureSet:

:discreteDirectionsOct: Integer representing the number of direction per octant.
:directionUxOct: A list of scalar representing the direction cosines along the 
                 X-axis. The length of the list should be equal to the value 
                 provided in ``discreteDirectionsOct``.
:directionUyOct: A list of scalar representing the direction cosines along the 
                 Y-axis. The length of the list should be equal to the value 
                 provided in ``discreteDirectionsOct``.
:directionUzOct: A list of scalar representing the direction cosines along the 
                 Z-axis. The length of the list should be equal to the value 
                 provided in ``discreteDirectionsOct``.
directionWeightsOct: A list of scalar representing the weight of the direction. 
                     The length of the list should be equal to the value 
                     provided in ``discreteDirectionsOct``.

