# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class ConductivityModel(OffbeatDict):
    """
    Mother class for conductivity models.
    """
    TYPE: ClassVar[str] = 'conductivity'

Conductivity = ConductivityModel  # alias

@offbeat_define
class Constant(Conductivity):
    """
    The `constant` conductivity model assigns a constant thermal conductivity
    value to all cells of the addressed material region.

    The preferred input style is a `conductivity` sub-dictionary with a `type`
    entry and a `value` keyword. The legacy keyword `k` in the material dictionary
    is still accepted for retro-compatibility but is deprecated.


    Options
    -------
    value : scalar
        Constant thermal conductivity value [W/m/K].
        (required: True)
    """
    TYPE: ClassVar[str] = 'constant'
    value: float | int

@offbeat_define
class UO2Matpro(Conductivity):
    """
    The `UO2Matpro` conductivity model implements a UO$_2$ thermal conductivity
    correlation inspired by MATPRO v11.

    The model depends on the local temperature and burnup field. It also accounts
    for porosity effects using the `densityFraction` read from the base material
    dictionary (and may additionally correct for porosity evolution if a `porosity`
    field is present), as well as for the presence of Gadolinium.


    Options
    -------
    burnupName : word
        Name of the burnup field looked up in the mesh.
        (default: 'Bu'; required: False)

    GdContent : scalar
        Initial gadolinium content used by the correlation.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)

    densityFraction : scalar
        Fraction of theoretical density.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)

    par1 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 0.0452; required: False)

    par2 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 0.000246; required: False)

    par3 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 0.00187; required: False)

    par4 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 1.1599; required: False)

    par5 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 0.9; required: False)

    par6 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 0.04; required: False)

    par7 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 0.038; required: False)

    par8 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 0.28; required: False)

    par9 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 396.0; required: False)

    par10 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 6380.0; required: False)

    par11 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: '3.5e9'; required: False)

    par12 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 2.0; required: False)

    par13 : scalar
        Optional correlation parameter override (retro-compatibility).
        (default: 16360.0; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Matpro'
    burnupName: str = 'Bu'
    GdContent: float | int | None = None
    densityFraction: float | int | None = None
    par1: float | int = 0.0452
    par2: float | int = 0.000246
    par3: float | int = 0.00187
    par4: float | int = 1.1599
    par5: float | int = 0.9
    par6: float | int = 0.04
    par7: float | int = 0.038
    par8: float | int = 0.28
    par9: float | int = 396.0
    par10: float | int = 6380.0
    par11: float | int = 3.5e9
    par12: float | int = 2.0
    par13: float | int = 16360.0

@offbeat_define
class Molybdenum(Conductivity):
    """
    The `Molybdenum` conductivity model provides a temperature-dependent thermal
    conductivity correlation for molybdenum.

    Parameters can optionally be provided in the `conductivity` sub-dictionary.
    If not provided, the default values embedded in the model are used.


    Options
    -------
    par1 : scalar
        Quadratic coefficient in the correlation.
        (default: 9.128e-06; required: False)

    par2 : scalar
        Linear coefficient in the correlation.
        (default: 0.04945; required: False)

    par3 : scalar
        Constant term in the correlation.
        (default: 152; required: False)
    """
    TYPE: ClassVar[str] = 'Molybdenum'
    par1: float | int = 9.128e-06
    par2: float | int = 0.04945
    par3: float | int = 152.0

@offbeat_define
class UO2Nfir(Conductivity):
    """
    The `UO2Nfir` conductivity model implements a UO$_2$ thermal conductivity
    correlation from NFIR.

    The model depends on the local temperature and on a burnup field. It also
    applies a porosity-related correction based on `densityFraction` and includes
    a Gd-content effect through `GdContent`.

    `densityFraction` and `GdContent` can be given as optional overrides in the
    `conductivity` sub-dictionary. If not provided there, they are expected in the
    material base/composition dictionary.


    Options
    -------
    burnupName : word
        Name of the burnup field looked up in the mesh.
        (default: 'Bu'; required: False)

    densityFraction : scalar
        Optional override of the fuel fraction of theoretical density used for the
        porosity correction.

        If not provided in the `conductivity` sub-dictionary, it is read from the
        material base/composition dictionary (recommended place).
        (required: False)

    GdContent : scalar
        Optional override of the initial gadolinium content used by the model.

        If not provided in the `conductivity` sub-dictionary, it is read from the
        material base/composition dictionary (recommended place).
        (required: False)
    """
    TYPE: ClassVar[str] = 'UO2Nfir'
    burnupName: str = 'Bu'
    densityFraction: float | int | None = None
    GdContent: float | int | None = None

@offbeat_define
class BufferParfume(Conductivity):
    """
    The `BufferParfume` conductivity model provides the thermal conductivity of
    the TRISO buffer layer as a function of its density.

    The model reads the density field (by default `rho`) and uses an initial and a
    theoretical density, together with an initial conductivity and a conductivity
    at theoretical density.

    Theoretical density can be provided as an optional override in the
    `conductivity` sub-dictionary. If not provided there, it is expected in the
    material base/composition dictionary.


    Options
    -------
    densityName : word
        Name of the density field looked up in the mesh.
        (default: 'rho'; required: False)

    initialDensity : scalar
        Initial density of the buffer material used by the model.
        (default: 1000.0; required: False)

    theoreticalDensity : scalar
        Theoretical density of the buffer material.

        Optional override: if not provided in the `conductivity` sub-dictionary,
        it is read from the material base/composition dictionary (recommended).
        (required: False)

    initialConductivity : scalar
        Initial thermal conductivity of the buffer material.
        (default: 0.5; required: False)

    theoreticalConductivity : scalar
        Thermal conductivity at theoretical density.
        (default: 4.0; required: False)
    """
    TYPE: ClassVar[str] = 'BufferParfume'
    densityName: str = 'rho'
    initialDensity: float | int = 1000.0
    theoreticalDensity: float | int | None = None
    initialConductivity: float | int = 0.5
    theoreticalConductivity: float | int = 4.0

@offbeat_define
class SiCParfume(Conductivity):
    """
    The `SiCParfume` conductivity model provides a temperature-dependent thermal
    conductivity correlation for SiC based on PARFUME.

    Parameters can optionally be provided in the `conductivity` sub-dictionary.
    If not provided, the default values embedded in the model are used.


    Options
    -------
    par1 : scalar
        Correlation parameter.
        (default: 17885.0; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 2.0; required: False)
    """
    TYPE: ClassVar[str] = 'SiCParfume'
    par1: float | int = 17885.0
    par2: float | int = 2.0

@offbeat_define
class ZircaloyRelap(Conductivity):
    """
    The `ZircaloyRelap` conductivity model provides a temperature-dependent
    thermal conductivity correlation for Zircaloy derived from RELAP.

    Parameters can optionally be provided in the `conductivity` sub-dictionary.
    If not provided, the default values embedded in the model are used.


    Options
    -------
    par1 : scalar
        Correlation parameter.
        (default: 7.51; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 0.0209; required: False)

    par3 : scalar
        Correlation parameter.
        (default: 1.45e-05; required: False)

    par4 : scalar
        Correlation parameter.
        (default: 7.67e-09; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyRelap'
    par1: float | int = 7.51
    par2: float | int = 0.0209
    par3: float | int = 1.45e-05
    par4: float | int = 7.67e-09

@offbeat_define
class HastelloyNSwindeman(Conductivity):
    """
    Thermal conductivity correlation for Hastelloy-N derived from the work of
    Swindeman.

    The correlation is linear in temperature and is expressed as:

    $$
    k = \mathrm{par1} + \mathrm{par2} \cdot (T - 273.15)
    $$


    Options
    -------
    par1 : scalar
        Constant term of the conductivity correlation.
        (default: 8.431; required: False)

    par2 : scalar
        Linear temperature coefficient of the conductivity correlation.
        (default: 0.0205; required: False)
    """
    TYPE: ClassVar[str] = 'HastelloyNSwindeman'
    par1: float | int = 8.431
    par2: float | int = 0.0205

@offbeat_define
class Steel1515TiTobbe(Conductivity):
    """
    Thermal conductivity correlation for 15-15Ti steel based on the formulation
    by Tobbe.

    The model provides a temperature-dependent thermal conductivity for
    austenitic 15-15Ti steel.


    Options
    -------
    par1 : scalar
        Correlation constant term.
        (default: 13.95; required: False)

    par2 : scalar
        Linear temperature coefficient of the correlation.
        (default: 0.01163; required: False)
    """
    TYPE: ClassVar[str] = 'Steel1515TiTobbe'
    par1: float | int = 13.95
    par2: float | int = 0.01163

@offbeat_define
class UO2Ifa601(Conductivity):
    """
    Thermal conductivity model for UO₂ fuel based on the IFA-601 formulation.

    The correlation depends on temperature and burnup and is commonly used for
    LWR fuel applications.


    Options
    -------
    burnupName : word
        Name of the burnup field used by the correlation.
        (default: 'Bu'; required: False)

    par1 : scalar
        Correlation parameter.
        (default: 40; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 5; required: False)

    par3 : scalar
        Correlation parameter.
        (default: 0.24; required: False)

    par4 : scalar
        Correlation parameter.
        (default: 0.457; required: False)

    par5 : scalar
        Correlation parameter.
        (default: 0.00433; required: False)

    par6 : scalar
        Correlation parameter.
        (default: 0.11; required: False)

    par7 : scalar
        Correlation parameter.
        (default: 1.32e-05; required: False)

    par8 : scalar
        Correlation parameter.
        (default: 0.00188; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Ifa601'
    burnupName: str = 'Bu'
    par1: float | int = 40.0
    par2: float | int = 5.0
    par3: float | int = 0.24
    par4: float | int = 0.457
    par5: float | int = 0.00433
    par6: float | int = 0.11
    par7: float | int = 1.32e-05
    par8: float | int = 0.00188

@offbeat_define
class UPuO2Brancheria(Conductivity):
    """
    Thermal conductivity correlation for MOX fuel (U-Pu-O₂) based on the
    formulation by Brancheria et al.

    The model accounts for porosity effects. If a porosity field is available,
    it is used directly; otherwise an initial porosity is derived from the
    density fraction.


    Options
    -------
    densityFraction : scalar
        Initial fraction of theoretical density used to derive porosity when no
        porosity field is available.

        Optional override: if not provided here, it is expected in the material
        base/composition dictionary.
        (default: 0.945; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Brancheria'
    densityFraction: float | int = 0.945

@offbeat_define
class UPuO2Kato(Conductivity):
    """
    Thermal conductivity correlation for MOX fuel (U-Pu-O₂) based on the
    formulation by Kato.


    Options
    -------
    densityFraction : scalar
        Fraction of theoretical density used by the correlation.

        Optional override: if not provided here, it is expected in the material
        base/composition dictionary.
        (default: 0.95; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Kato'
    densityFraction: float | int = 0.95

@offbeat_define
class UPuO2LanningBeyer(Conductivity):
    """
    Thermal conductivity correlation for MOX fuel (U-Pu-O₂) based on the
    Lanning–Beyer formulation.


    Options
    -------
    densityFraction : scalar
        Fraction of theoretical density used by the correlation.

        Optional override: if not provided here, it is expected in the material
        base/composition dictionary.
        (default: 0.95; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2LanningBeyer'
    densityFraction: float | int = 0.95

@offbeat_define
class MaUPuO2Magni(Conductivity):
    """
    Thermal conductivity correlation for MA-bearing MOX fuel based on the
    formulation by Magni.

    The model is intended for fuels containing minor actinides.


    Options
    -------
    densityFraction : scalar
        Fraction of theoretical density used by the correlation.

        Optional override: if not provided here, it is expected in the material
        base/composition dictionary.
        (default: 0.95; required: False)
    """
    TYPE: ClassVar[str] = 'MaUPuO2Magni'
    densityFraction: float | int = 0.95

@offbeat_define
class UPuO2Philipponeau(Conductivity):
    """
    Thermal conductivity correlation for MOX fuel (U-Pu-O₂) based on the
    formulation by Philipponneau.


    Options
    -------
    densityFraction : scalar
        Fraction of theoretical density used by the correlation.

        Optional override: if not provided here, it is expected in the material
        base/composition dictionary.
        (default: 0.95; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Philipponeau'
    densityFraction: float | int = 0.95
