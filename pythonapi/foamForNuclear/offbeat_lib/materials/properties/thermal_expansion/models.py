# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class ThermalExpansionModel(OffbeatDict):
    """
    Mother class for thermalExpansion models.


    Options
    -------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)
    """
    TYPE: ClassVar[str] = 'thermalExpansion'
    Tref: float | int = 293.0

ThermalExpansion = ThermalExpansionModel  # alias

@offbeat_define
class Constant(ThermalExpansion):
    """
    Model for constant thermal expansion coefficient. The value is read from
    dictionary.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    value : scalar
        Constant thermal expansion value [-/K].
        (required: True)
    """
    TYPE: ClassVar[str] = 'constant'
    value: float | int

@offbeat_define
class Steel1515TiGehr(ThermalExpansion):
    """
    Class modelling the thermal expansion of 15-15 Ti cladding through the
    correlation proposed by Gehr (1973).

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: -0.0003101; required: False)

    par2 : scalar
        Option read by the model.
        (default: 1.545e-05; required: False)

    par3 : scalar
        Option read by the model.
        (default: 2.75e-09; required: False)
    """
    TYPE: ClassVar[str] = 'Steel1515TiGehr'
    par1: float | int = -0.0003101
    par2: float | int = 1.545e-05
    par3: float | int = 2.75e-09

@offbeat_define
class UPuO2Lemehov(ThermalExpansion):
    """
    Class modelling thermal expansion of (U,Pu)O2 MOX fuel from Lemehov(2020).
    This correlation is valid for:
    - O/M ratios between 1.94 and 2;
    - Pu/HM ratios lower than 60 %.

    Source :
    https://re.public.polimi.it/retrieve/handle/11311/1172415/619939/Technical%20Report_INSPYRE_WP7-D7.2_%282020%29.pdf

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    b0 : scalar
        Option read by the model.
        (default: -0.308; required: False)

    b1 : scalar
        Option read by the model.
        (default: 3.4303; required: False)

    b2 : scalar
        Option read by the model.
        (default: -1.9157; required: False)

    b3 : scalar
        Option read by the model.
        (default: 3.4636; required: False)

    by : scalar
        Option read by the model.
        (default: 3.98; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Lemehov'
    b0: float | int = -0.308
    b1: float | int = 3.4303
    b2: float | int = -1.9157
    b3: float | int = 3.4636
    by: float | int = 3.98

@offbeat_define
class MaUPuO2Kato(ThermalExpansion):
    """
    Class modelling  isotropic thermal expansion in MA-MOX fuel (for Pu = 0.3)
    Source : JNM 469 (2016) 223-227.
    url : https://doi.org/10.1016/j.jnucmat.2015.11.048

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    oxygenMetalRatio : scalar
        Oxygen to metal ratio.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)
    """
    TYPE: ClassVar[str] = 'MaUPuO2Kato'
    oxygenMetalRatio: float | int | None = None

@offbeat_define
class UPuO2Martin(ThermalExpansion):
    """
    Class modelling thermal expansion of MOX fuel according to Martin
    correlation.
    Source : "The thermal expansion of solid UO2 and (U, Pu) mixed oxides — a
    review and recommendations", D.G.Martin

    url : https://doi.org/10.1016/0022-3115(88)90315-7

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    oxygenMetalRatio : scalar
        Oxygen to metal ratio.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)

    par1 : scalar
        Option read by the model.
        (default: 9.828e-06; required: False)

    par2 : scalar
        Option read by the model.
        (default: -6.39e-10; required: False)

    par3 : scalar
        Option read by the model.
        (default: 1.33e-12; required: False)

    par4 : scalar
        Option read by the model.
        (default: -1.757e-17; required: False)

    par5 : scalar
        Option read by the model.
        (default: 1.1833e-05; required: False)

    par6 : scalar
        Option read by the model.
        (default: -5.013e-09; required: False)

    par7 : scalar
        Option read by the model.
        (default: 3.756e-12; required: False)

    par8 : scalar
        Option read by the model.
        (default: -6.125e-17; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Martin'
    oxygenMetalRatio: float | int | None = None
    par1: float | int = 9.828e-06
    par2: float | int = -6.39e-10
    par3: float | int = 1.33e-12
    par4: float | int = -1.757e-17
    par5: float | int = 1.1833e-05
    par6: float | int = -5.013e-09
    par7: float | int = 3.756e-12
    par8: float | int = -6.125e-17

@offbeat_define
class UPuO2Matpro(ThermalExpansion):
    """
    Class modelling thermal expansion of (U,Pu)O2 MOX fuel from Matprov11.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: -0.00039735; required: False)

    par2 : scalar
        Option read by the model.
        (default: 8.4955e-06; required: False)

    par3 : scalar
        Option read by the model.
        (default: 2.1513e-09; required: False)

    par4 : scalar
        Option read by the model.
        (default: 3.7143e-16; required: False)

    par5 : scalar
        Option read by the model.
        (default: -0.0004972; required: False)

    par6 : scalar
        Option read by the model.
        (default: 7.107e-06; required: False)

    par7 : scalar
        Option read by the model.
        (default: 2.581e-09; required: False)

    par8 : scalar
        Option read by the model.
        (default: 1.14e-13; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Matpro'
    par1: float | int = -0.00039735
    par2: float | int = 8.4955e-06
    par3: float | int = 2.1513e-09
    par4: float | int = 3.7143e-16
    par5: float | int = -0.0004972
    par6: float | int = 7.107e-06
    par7: float | int = 2.581e-09
    par8: float | int = 1.14e-13

@offbeat_define
class ZircaloyMatpro(ThermalExpansion):
    """
    Class modelling thermal expansion of UO2 fuel from Matprov11.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 4.441e-06; required: False)

    par2 : scalar
        Option read by the model.
        (default: 0.001238; required: False)

    par3 : scalar
        Option read by the model.
        (default: 6.721e-06; required: False)

    par4 : scalar
        Option read by the model.
        (default: 0.002073; required: False)

    par5 : scalar
        Option read by the model.
        (default: 9.7e-06; required: False)

    par6 : scalar
        Option read by the model.
        (default: 0.011; required: False)

    par7 : scalar
        Option read by the model.
        (default: 0.00945; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyMatpro'
    par1: float | int = 4.441e-06
    par2: float | int = 0.001238
    par3: float | int = 6.721e-06
    par4: float | int = 0.002073
    par5: float | int = 9.7e-06
    par6: float | int = 0.011
    par7: float | int = 0.00945

@offbeat_define
class Molybdenum(ThermalExpansion):
    """
    Class modelling thermal expansion of Molybdenum.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 4.985e-06; required: False)

    par2 : scalar
        Option read by the model.
        (default: 6.667e-10; required: False)
    """
    TYPE: ClassVar[str] = 'Molybdenum'
    par1: float | int = 4.985e-06
    par2: float | int = 6.667e-10

@offbeat_define
class BufferParfume(ThermalExpansion):
    """
    Class modelling the thermal expansion of Buffer through the code Parfume.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 5.0; required: False)

    par2 : scalar
        Option read by the model.
        (default: 0.11; required: False)

    par3 : scalar
        Option read by the model.
        (default: 400.0; required: False)

    par4 : scalar
        Option read by the model.
        (default: 700.0; required: False)
    """
    TYPE: ClassVar[str] = 'BufferParfume'
    par1: float | int = 5.0
    par2: float | int = 0.11
    par3: float | int = 400.0
    par4: float | int = 700.0

@offbeat_define
class PyCParfume(ThermalExpansion):
    """
    Class modelling the transversly isotropic thermal expansion of PyC through
    the code Parfume.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 30.0; required: False)

    par2 : scalar
        Option read by the model.
        (default: 37.5; required: False)

    par3 : scalar
        Option read by the model.
        (default: 0.11; required: False)

    par4 : scalar
        Option read by the model.
        (default: 673.0; required: False)

    par5 : scalar
        Option read by the model.
        (default: 700.0; required: False)

    par6 : scalar
        Option read by the model.
        (default: 36.0; required: False)
    """
    TYPE: ClassVar[str] = 'PyCParfume'
    par1: float | int = 30.0
    par2: float | int = 37.5
    par3: float | int = 0.11
    par4: float | int = 673.0
    par5: float | int = 700.0
    par6: float | int = 36.0

@offbeat_define
class SiCParfume(ThermalExpansion):
    """
    Class modelling the thermal expansion of SiC with a constant value 4.9×10^−6 /K from
    Parfume.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)
    """
    TYPE: ClassVar[str] = 'SiCParfume'

@offbeat_define
class UO2Relap(ThermalExpansion):
    """
    Class modelling thermal expansion of UO2 fuel derived from Relap.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    ED : scalar
        Option read by the model.
        (default: 1.32e-19; required: False)

    K1 : scalar
        Option read by the model.
        (default: 9.8e-06; required: False)

    K2 : scalar
        Option read by the model.
        (default: 0.00261; required: False)

    K3 : scalar
        Option read by the model.
        (default: 0.316; required: False)

    k : scalar
        Option read by the model.
        (default: 1.38e-23; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Relap'
    ED: float | int = 1.32e-19
    K1: float | int = 9.8e-06
    K2: float | int = 0.00261
    K3: float | int = 0.316
    k: float | int = 1.38e-23

@offbeat_define
class SiCSnead(ThermalExpansion):
    """
    Class modelling the thermal expansion of SiC with a value from Snead. the mean
    thermal expansion coefficient is calculated and the The method of Niffenegger and
    Reichlin (2012)
    is employed to convert the mean thermal expansion values into instantaneous values.

    Snead model was used for UN TRIOS by Collin(2014).

    -M. Niffenegger and K. Reichlin. The proper use of thermal expansion coefficients
    in finite element calculations. Nuclear Engineering and Design, 243:356–359, 2012.
    -L. L. Snead, T. Nozawa, Y. Katoh, T.-S. Byun, S. Kondo, and D. A. Petti.
    Handbook of sic properties for fuel performance modeling. Journal of Nuclear Materials,
    371:329–377, 2007.
    -Blaise P. Collin. Modeling and analysis of UN TRISO fuel for LWR application
    using the Parfume code. Journal of Nuclear Materials, 451(1):65–77, 2014.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    Treference : scalar
        Option read by the model.
        (default: 298.15; required: False)

    par1 : scalar
        Option read by the model.
        (default: -1.8276; required: False)

    par2 : scalar
        Option read by the model.
        (default: 0.0178; required: False)

    par3 : scalar
        Option read by the model.
        (default: -1.5544e-05; required: False)

    par4 : scalar
        Option read by the model.
        (default: 4.5246e-09; required: False)

    par5 : scalar
        Option read by the model.
        (default: '5e-6'; required: False)
    """
    TYPE: ClassVar[str] = 'SiCSnead'
    Treference: float | int = 298.15
    par1: float | int = -1.8276
    par2: float | int = 0.0178
    par3: float | int = -1.5544e-05
    par4: float | int = 4.5246e-09
    par5: float | int = 5e-6

@offbeat_define
class HastelloyNSwindeman(ThermalExpansion):
    """
    Class modelling thermal expansion of H-N from Swinderman.

    Options (inherited)
    -------------------
    Tref : scalar
        Reference temperature for thermal strain.
        (default: 293.0; required: False)

    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 0.005291; required: False)

    par2 : scalar
        Option read by the model.
        (default: 9.682; required: False)

    par3 : scalar
        Option read by the model.
        (default: 107.8; required: False)
    """
    TYPE: ClassVar[str] = 'HastelloyNSwindeman'
    par1: float | int = 0.005291
    par2: float | int = 9.682
    par3: float | int = 107.8
