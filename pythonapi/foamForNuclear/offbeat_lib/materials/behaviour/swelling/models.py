# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class SwellingModel(OffbeatDict):
    """
    Parent class for swelling model.

    The parent class typeName is "none". If selected in the solverDict the
    swelling is not considered, but the swelling field 'epsilonSwelling' is
    created and added to registry.

    If the 'epsilonSwelling' file is present in the starting time folder,
    the internal field and boundary conditions are set according to the file.

    Otherwise, the field is set to 0 in each cell, with zeroGradient BCs in all
    non-empty/-wedge patches.
    """
    TYPE: ClassVar[str] = 'none'

Swelling = SwellingModel  # alias

@offbeat_define
class ConstantRate(Swelling):
    """
    Class caculating the swelling with a provided constant swelling rate.


    Options
    -------
    swellingRate : scalar
        Constant swelling rate coefficient.
        (required: True)

    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)
    """
    TYPE: ClassVar[str] = 'constantRate'
    swellingRate: float | int
    fastFluenceName: str = 'fastFluence'

@offbeat_define
class PyCCorrelation(Swelling):
    """
    Class handling irradiation-induced dimensional change eigenstrain phenomenon
    for PyC with the correlation below for both radial and tangential directions.


    Options
    -------
    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    sphereCoordinate : bool
        If true, interpret coordinates as spherical.
        (default: True; required: False)

    fluxConversionFactor : scalar
        Conversion factor applied to the flux/fluence input.
        (default: 1.0; required: False)

    radialCoefficients : dict
        Model data table.
        (required: True)

    tangentialCoefficients : dict
        Model data table.
        (required: True)
    """
    TYPE: ClassVar[str] = 'PyCCorrelation'
    fastFluenceName: str = 'fastFluence'
    sphereCoordinate: bool = True
    fluxConversionFactor: float | int = 1.0
    radialCoefficients: dict[str, Any]
    tangentialCoefficients: dict[str, Any]

@offbeat_define
class UPuO2Fbr(Swelling):
    """
    Class modelling swelling strain for FBR MOX fuel from Dienst at al.
    Sources :
    - Dienst, W., Muelle-Lyda, I., Zimmerman, H. Swelling, densification and
    creep of oxide and carbide fuels under irradiation, in: International
    conference on fast breeder reactor performance, 5-8 March 1979, Monterey,
    California, USA, 1979, pp. 166–175.
    - "Validazione e analisi dei gas di fissione in combustibili MOX ad
    elevato burnup", D.Rozzia, N.Forgione, A. Ardizzone, A. Del Nevo.


    Options
    -------
    burnupName : word
        Name of the burnup field.
        (default: 'Bu'; required: False)

    gapWidthName : word
        Name of the gap-width field.
        (default: 'gapWidth'; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 0.02; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 0.012; required: False)

    par3 : scalar
        Correlation parameter.
        (default: 0.0065; required: False)

    outerPatch : word
        Name of the outer boundary patch.
        (required: True)
    """
    TYPE: ClassVar[str] = 'UPuO2Fbr'
    burnupName: str = 'Bu'
    gapWidthName: str = 'gapWidth'
    par1: float | int = 0.02
    par2: float | int = 0.012
    par3: float | int = 0.0065
    outerPatch: str

@offbeat_define
class FeCrAl(Swelling):
    """
    Class handling swelling phenomenon for FeCrAl derived from Bison manual.


    Options
    -------
    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 4.5e-29; required: False)
    """
    TYPE: ClassVar[str] = 'FeCrAl'
    fastFluenceName: str = 'fastFluence'
    par1: float | int = 4.5e-29

@offbeat_define
class UO2Frapcon(Swelling):
    """
    Class modelling swelling phenomenon derived from Frapcon.


    Options
    -------
    burnupName : word
        Name of the burnup field.
        (default: 'Bu'; required: False)

    densityName : word
        Name of the density field.
        (default: 'rho'; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 6000; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 80000; required: False)

    par3 : scalar
        Correlation parameter.
        (default: 29740000000.0; required: False)

    par4 : scalar
        Correlation parameter.
        (default: 2.315e-23; required: False)

    par5 : scalar
        Correlation parameter.
        (default: 86.4; required: False)

    par6 : scalar
        Correlation parameter.
        (default: 3.211e-23; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Frapcon'
    burnupName: str = 'Bu'
    densityName: str = 'rho'
    par1: float | int = 6000.0
    par2: float | int = 80000.0
    par3: float | int = 29740000000.0
    par4: float | int = 2.315e-23
    par5: float | int = 86.4
    par6: float | int = 3.211e-23

@offbeat_define
class Steel1515TiAim1(Swelling):
    """
    Class modelling the void swelling growth phenomenon for 15-15 Ti cladding
    material according to the 'AIM1 15-15 Ti' correlation.
    Source : "Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative
    Lead Cooled FBR - L.Luzzi et al. - 2014"


    Options
    -------
    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 1.3e-05; required: False)

    par2 : scalar
        Correlation parameter.
        (default: -2.5; required: False)

    par3 : scalar
        Correlation parameter.
        (default: 490; required: False)

    par4 : scalar
        Correlation parameter.
        (default: 100; required: False)

    par5 : scalar
        Correlation parameter.
        (default: 3.9; required: False)
    """
    TYPE: ClassVar[str] = 'Steel1515TiAim1'
    fastFluenceName: str = 'fastFluence'
    par1: float | int = 1.3e-05
    par2: float | int = -2.5
    par3: float | int = 490.0
    par4: float | int = 100.0
    par5: float | int = 3.9

@offbeat_define
class ZircaloyBison(Swelling):
    """
    Class modelling the irradiation growth phenomenon derived from Moose documentation.


    Options
    -------
    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    A : scalar
        Correlation parameter.
        (default: '3e-20'; required: False)

    n : scalar
        Correlation parameter.
        (default: 0.794; required: False)

    cladType : word
        Cladding type selector used by the model.
        (default: 'ESCORE'; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyBison'
    fastFluenceName: str = 'fastFluence'
    A: float | int = 3e-20
    n: float | int = 0.794
    cladType: str = 'ESCORE'

@offbeat_define
class Steel1515TiGeneralized(Swelling):
    """
    Class modelling the void swelling growth phenomenon for 15-15 Ti cladding
    material according to the 'Generalized 15-15 Ti' correlation.
    Source : "Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative
    Lead Cooled FBR - L.Luzzi et al. - 2014"


    Options
    -------
    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 0.0015; required: False)

    par2 : scalar
        Correlation parameter.
        (default: -2.5; required: False)

    par3 : scalar
        Correlation parameter.
        (default: 450; required: False)

    par4 : scalar
        Correlation parameter.
        (default: 100; required: False)

    par5 : scalar
        Correlation parameter.
        (default: 2.75; required: False)
    """
    TYPE: ClassVar[str] = 'Steel1515TiGeneralized'
    fastFluenceName: str = 'fastFluence'
    par1: float | int = 0.0015
    par2: float | int = -2.5
    par3: float | int = 450.0
    par4: float | int = 100.0
    par5: float | int = 2.75

@offbeat_define
class ZircaloyMatpro(Swelling):
    """
    Class modelling the irradiation growth phenomenon according to Matpro
    correlations.


    Options
    -------
    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    coldWork : scalar
        Cold work fraction/indicator used by the model.
        (default: 0.0; required: False)

    A : scalar
        Correlation parameter.
        (default: 1.407e-16; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 240.8; required: False)

    fz : scalar
        Correlation parameter.
        (default: 0.05; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyMatpro'
    fastFluenceName: str = 'fastFluence'
    coldWork: float | int = 0.0
    A: float | int = 1.407e-16
    par1: float | int = 240.8
    fz: float | int = 0.05

@offbeat_define
class UO2Matpro(Swelling):
    """
    Class modelling swelling phenomenon derived from Frapcon.


    Options
    -------
    burnupName : word
        Name of the burnup field.
        (default: 'Bu'; required: False)

    densityName : word
        Name of the density field.
        (default: 'rho'; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 5.577e-05; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 1.96e-31; required: False)

    par3 : scalar
        Correlation parameter.
        (default: 2800; required: False)

    par4 : scalar
        Correlation parameter.
        (default: 11.73; required: False)

    par5 : scalar
        Correlation parameter.
        (default: -0.0162; required: False)

    par6 : scalar
        Correlation parameter.
        (default: -0.0178; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Matpro'
    burnupName: str = 'Bu'
    densityName: str = 'rho'
    par1: float | int = 5.577e-05
    par2: float | int = 1.96e-31
    par3: float | int = 2800.0
    par4: float | int = 11.73
    par5: float | int = -0.0162
    par6: float | int = -0.0178

@offbeat_define
class BufferParfume(Swelling):
    """
    Class handling irradiation-induced dimensional change phenomenon for Buffer
    derived from Parfume manual.


    Options
    -------
    densityName : word
        Name of the density field.
        (default: 'rho'; required: False)

    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    fastFluxName : word
        Name of the fast flux field.
        (default: 'fastFlux'; required: False)

    fluxConversionFactor : scalar
        Conversion factor applied to the flux/fluence input.
        (default: 1.0; required: False)

    tempValues : scalar
        Model data table.
        (required: True)

    BAFValues_r : scalar
        Model data table.
        (required: True)

    densityValues : scalar
        Model data table.
        (required: True)

    strain_iso : scalar
        Model option.
        (required: True)

    a1_r : scalar
        Model data table.
        (required: True)

    a2_r : scalar
        Model data table.
        (required: True)

    a3_r : scalar
        Model data table.
        (required: True)

    a4_r : scalar
        Model data table.
        (required: True)
    """
    TYPE: ClassVar[str] = 'BufferParfume'
    densityName: str = 'rho'
    fastFluenceName: str = 'fastFluence'
    fastFluxName: str = 'fastFlux'
    fluxConversionFactor: float | int = 1.0
    tempValues: float | int
    BAFValues_r: float | int
    densityValues: float | int
    strain_iso: float | int
    a1_r: float | int
    a2_r: float | int
    a3_r: float | int
    a4_r: float | int

@offbeat_define
class PyCParfume(Swelling):
    """
    Class handling irradiation-induced dimensional change eigenstrain phenomenon
    for PyC derived from Parfume manual. This class is the general situation for
    the class: swellingParfumeBuffer (BAF=1).


    Options
    -------
    densityName : word
        Name of the density field.
        (default: 'rho'; required: False)

    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    fastFluxName : word
        Name of the fast flux field.
        (default: 'fastFlux'; required: False)

    asFabricatedAnisotropy : scalar
        Model option.
        (default: 1.0; required: False)

    sphereCoordinate : bool
        If true, interpret coordinates as spherical.
        (default: True; required: False)

    fluxConversionFactor : scalar
        Conversion factor applied to the flux/fluence input.
        (default: 1.0; required: False)

    tempValues : scalar
        Model data table.
        (required: True)

    BAFValues_r : scalar
        Model data table.
        (required: True)

    BAFValues_t : scalar
        Model data table.
        (required: True)

    densityValues : scalar
        Model data table.
        (required: True)

    strain_iso : scalar
        Model option.
        (required: True)

    strain_r_t : scalar
        Model data table.
        (required: True)

    a1_r : scalar
        Model data table.
        (required: True)

    a2_r : scalar
        Model data table.
        (required: True)

    a3_r : scalar
        Model data table.
        (required: True)

    a4_r : scalar
        Model data table.
        (required: True)

    a1_t : scalar
        Model data table.
        (required: True)

    a2_t : scalar
        Model data table.
        (required: True)

    a3_t : scalar
        Model data table.
        (required: True)

    a4_t : scalar
        Model data table.
        (required: True)
    """
    TYPE: ClassVar[str] = 'PyCParfume'
    densityName: str = 'rho'
    fastFluenceName: str = 'fastFluence'
    fastFluxName: str = 'fastFlux'
    asFabricatedAnisotropy: float | int = 1.0
    sphereCoordinate: bool = True
    fluxConversionFactor: float | int = 1.0
    tempValues: float | int
    BAFValues_r: float | int
    BAFValues_t: float | int
    densityValues: float | int
    strain_iso: float | int
    strain_r_t: float | int
    a1_r: float | int
    a2_r: float | int
    a3_r: float | int
    a4_r: float | int
    a1_t: float | int
    a2_t: float | int
    a3_t: float | int
    a4_t: float | int

@offbeat_define
class HastelloyNWrightSham(Swelling):
    """
    Class modelling the void swelling phenomenon for HastelloyN cladding
    material according to Wright and Sham correlations.


    Options
    -------
    fastFluenceName : word
        Name of the fast fluence field.
        (default: 'fastFluence'; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 0.9845; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 0.4385; required: False)

    par3 : scalar
        Correlation parameter.
        (default: -0.981; required: False)

    par4 : scalar
        Correlation parameter.
        (default: 490.0; required: False)
    """
    TYPE: ClassVar[str] = 'HastelloyNWrightSham'
    fastFluenceName: str = 'fastFluence'
    par1: float | int = 0.9845
    par2: float | int = 0.4385
    par3: float | int = -0.981
    par4: float | int = 490.0
