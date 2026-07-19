from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import OpenFOAMDict


class CoolantChannel(Patch):
    """
    Full multi-phase coolant channel boundary condition
    (``coolantChannel`` fvPatchScalarField).

    Two operating modes, selected by *enthalpy_model*:

    - ``enthalpy_model=True`` (default): bulk coolant temperature is computed
      by integrating the homogeneous equilibrium model (HEM) along the pin.
      Requires *inlet_temperature*.

    - ``enthalpy_model=False``: coolant bulk temperature is prescribed from a
      time × axial-location table (*time_points*, *axial_locations*, *T0_data*).
      Does not need *inlet_temperature*.

    Correlations default to the C++ defaults (see below); only the ones passed
    explicitly are written to the dict — unset keys fall back to C++ defaults.

    Options
    -------
    value : float
        Initial temperature on the patch [K].

    coolant_pressure : float
        Coolant pressure [Pa] (uniform, constant in time).

    mass_flow_rate : float
        Mass flow rate [kg/s] (uniform, constant in time).

    flow_area : float
        Cross-sectional flow area of the coolant channel [m²].

    hydraulic_diameter : float
        Hydraulic diameter of the coolant channel [m].

    inlet_temperature : float
        Coolant inlet temperature [K]. Required when *enthalpy_model=True* or
        when *chf=``"EPRI"``* (default CHF model).

    enthalpy_model : bool
        Use HEM enthalpy model for bulk T evolution (default True).

    time_points : list[float]
        Time points for the prescribed T0 table [s user-time].
        Required when *enthalpy_model=False*.

    axial_locations : list[float]
        Axial positions for the prescribed T0 table.
        Required when *enthalpy_model=False*.

    T0_data : list[list[float]]
        Coolant temperature table [K]. Shape (len(time_points), len(axial_locations)).
        Required when *enthalpy_model=False*.

    flow_direction : float
        Sign convention for axial flow direction (default 1.0).

    max_picard_iter : int
        Maximum Picard iterations for HEM convergence (default 10).

    picard_tol : float
        Convergence tolerance for Picard iterations (default 1.0).

    kappa : str
        Name of the solid thermal-conductivity field (default ``"k"``).

    relax : float
        Under-relaxation factor for HTC update (default 1.0).

    Correlation keys (all optional — if None, C++ default is used):
        single_phase        : default ``"Gnielinski"``
        onb                 : default ``"Basu"``
        subcooled_boiling   : default ``"Rohsenow"``
        saturated_boiling   : default ``"Rohsenow"``
        film_boiling        : default ``"Frederking"``
        chf                 : default ``"EPRI"``
        leidenfrost         : default ``"GroeneveldStewartCorrected"``
        transition_boiling  : default ``"McDonoughMilichKing"``
    """

    def __init__(
            self,
            value: float,
            coolant_pressure: float,
            mass_flow_rate: float,
            flow_area: float,
            hydraulic_diameter: float,
            inlet_temperature: float = None,
            enthalpy_model: bool = True,
            time_points: list = None,
            axial_locations: list = None,
            T0_data: list = None,
            flow_direction: float = 1.0,
            max_picard_iter: int = 10,
            picard_tol: float = 1.0,
            kappa: str = "k",
            relax: float = 1.0,
            single_phase: str = None,
            onb: str = None,
            subcooled_boiling: str = None,
            saturated_boiling: str = None,
            film_boiling: str = None,
            chf: str = None,
            leidenfrost: str = None,
            transition_boiling: str = None,
        ):
        super().__init__(type="coolantChannel", value=value)

        check_type("coolant_pressure", coolant_pressure, (int, float))
        check_type("mass_flow_rate", mass_flow_rate, (int, float))
        check_type("flow_area", flow_area, (int, float))
        check_type("hydraulic_diameter", hydraulic_diameter, (int, float))
        check_type("enthalpy_model", enthalpy_model, bool)
        check_type("kappa", kappa, str)
        check_type("relax", relax, (int, float))

        if enthalpy_model and inlet_temperature is None:
            raise ValueError("inlet_temperature is required when enthalpy_model=True")
        if not enthalpy_model and (time_points is None or axial_locations is None or T0_data is None):
            raise ValueError(
                "time_points, axial_locations and T0_data are required when enthalpy_model=False"
            )

        self.__setitem__("kappa", kappa)
        self.__setitem__("flowArea", flow_area)
        self.__setitem__("hydraulicDiameter", hydraulic_diameter)
        self.__setitem__("flowDirection", flow_direction)
        self.__setitem__("enthalpyModel", enthalpy_model)
        self.__setitem__("coolantPressure", f"uniform {coolant_pressure}")
        self.__setitem__("massFlowRate", f"uniform {mass_flow_rate}")
        self.__setitem__("relax", relax)

        if inlet_temperature is not None:
            self.__setitem__("inletTemperature", inlet_temperature)

        if max_picard_iter != 10:
            self.__setitem__("maxPicardIter", max_picard_iter)
        if picard_tol != 1.0:
            self.__setitem__("picardTol", picard_tol)

        # Correlations subdictionary — only write explicitly set values
        corr_items = {}
        if single_phase is not None:
            corr_items["singlePhase"] = single_phase
        if onb is not None:
            corr_items["onb"] = onb
        if subcooled_boiling is not None:
            corr_items["subcooledBoiling"] = subcooled_boiling
        if saturated_boiling is not None:
            corr_items["saturatedBoiling"] = saturated_boiling
        if film_boiling is not None:
            corr_items["filmBoiling"] = film_boiling
        if chf is not None:
            corr_items["chf"] = chf
        if leidenfrost is not None:
            corr_items["leidenfrost"] = leidenfrost
        if transition_boiling is not None:
            corr_items["transitionBoiling"] = transition_boiling

        # correlations subDict is always required by C++
        self.__setitem__("correlations", OpenFOAMDict(items=corr_items, name="correlations"))

        if not enthalpy_model:
            t0_rows = " ".join(
                f"( {' '.join(str(v) for v in row)} )" for row in T0_data
            )
            axial = OpenFOAMDict(name="axialProfileDict")
            axial["timePoints"] = time_points
            axial["axialLocations"] = axial_locations
            axial["T0Data"] = f"( {t0_rows} )"
            self.__setitem__("axialProfileDict", axial)
