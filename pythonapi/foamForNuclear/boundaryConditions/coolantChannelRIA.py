from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import OpenFOAMDict


class CoolantChannelRIA(Patch):
    """
    Coolant channel boundary condition for Reactivity Initiated Accident (RIA)
    transients (``coolantChannelRIA`` fvPatchScalarField).

    Applies a Bessiron (2007) boiling-crisis criterion: crisis is triggered when
    the vaporized water layer thickness exceeds *thicknessBessironRIA*, not when
    wall heat flux exceeds CHF. Coolant bulk temperature is prescribed via a
    time × axial-location table (``axialProfileDict``).

    Options
    -------
    value : float
        Initial temperature on the patch [K].

    coolant_pressure : float
        Coolant pressure [Pa] (uniform, constant in time).

    mass_flow_rate : float
        Mass flow rate [kg/s] (uniform, constant in time).

    time_points : list[float]
        Time points for the axial coolant-temperature table [s user-time].

    axial_locations : list[float]
        Axial positions of the table rows, normalised 0–1 or in metres
        depending on how the profile is specified.

    T0_data : list[list[float]]
        Coolant bulk temperature table [K]. Outer index = time, inner index =
        axial location. Shape must be (len(time_points), len(axial_locations)).

    hydraulic_diameter : float
        Hydraulic diameter of the coolant channel [m].

    flow_area : float
        Cross-sectional flow area of the coolant channel [m²].

    kappa : str
        Name of the solid thermal-conductivity field (default ``"k"``).

    relax : float
        Under-relaxation factor for the HTC update (default 0.15).

    pulse_start_time : float
        User time at which the RIA pulse begins [s]. The transient conduction
        clock only starts from this moment; single-phase natural convection is
        used before it (default 0.0).

    single_phase_model : str
        Single-phase HTC correlation (default ``"ChurchillChu"``).

    film_boiling_model : str
        Film-boiling HTC correlation (default ``"Sakurai"``).

    film_boiling_ria_correction : float
        Film-boiling HTC multiplier at the moment of boiling crisis (default 5.0).
        Linearly decays to 1 over *film_boiling_ria_correction_duration* seconds.

    film_boiling_ria_correction_duration : float
        Duration over which *film_boiling_ria_correction* decays to 1 [s]
        (default 15.0).

    leidenfrost_dt_bessiron : float
        Temperature increment above T_sat defining the Leidenfrost point [K]
        (default 450.0).

    onb_dt_bessiron : float
        Temperature increment above T_sat defining the ONB / CHF point [K]
        (default 20.0).

    thickness_bessiron_ria : float
        Vaporized layer thickness that triggers boiling crisis [m]
        (default 3e-5).

    reset_delta_t_bessiron : float
        Wall–coolant temperature margin for post-crisis reset [K] (default 1.0).
    """

    def __init__(
            self,
            value: float,
            coolant_pressure: float,
            mass_flow_rate: float,
            time_points: list,
            axial_locations: list,
            T0_data: list,
            hydraulic_diameter: float,
            flow_area: float,
            kappa: str = "k",
            relax: float = 0.15,
            pulse_start_time: float = 0.0,
            single_phase_model: str = "ChurchillChu",
            film_boiling_model: str = "Sakurai",
            film_boiling_ria_correction: float = 5.0,
            film_boiling_ria_correction_duration: float = 15.0,
            leidenfrost_dt_bessiron: float = 450.0,
            onb_dt_bessiron: float = 20.0,
            thickness_bessiron_ria: float = 3e-5,
            reset_delta_t_bessiron: float = 1.0,
        ):
        super().__init__(type="coolantChannelRIA", value=value)

        check_type("coolant_pressure", coolant_pressure, (int, float))
        check_type("mass_flow_rate", mass_flow_rate, (int, float))
        check_type("time_points", time_points, list)
        check_type("axial_locations", axial_locations, list)
        check_type("T0_data", T0_data, list)
        check_type("hydraulic_diameter", hydraulic_diameter, (int, float))
        check_type("flow_area", flow_area, (int, float))
        check_type("kappa", kappa, str)
        check_type("relax", relax, (int, float))
        check_type("pulse_start_time", pulse_start_time, (int, float))
        check_type("single_phase_model", single_phase_model, str)
        check_type("film_boiling_model", film_boiling_model, str)

        self.__setitem__("kappa", kappa)
        self.__setitem__("hydraulicDiameter", hydraulic_diameter)
        self.__setitem__("flowArea", flow_area)
        self.__setitem__("relax", relax)
        self.__setitem__("pulseStartTime", pulse_start_time)
        self.__setitem__("coolantPressure", f"uniform {coolant_pressure}")
        self.__setitem__("massFlowRate", f"uniform {mass_flow_rate}")
        self.__setitem__("filmBoilingRIACorrection", film_boiling_ria_correction)
        self.__setitem__("filmBoilingRIACorrectionDuration", film_boiling_ria_correction_duration)
        self.__setitem__("leidenfrostdTBessiron", leidenfrost_dt_bessiron)
        self.__setitem__("onbdTBessiron", onb_dt_bessiron)
        self.__setitem__("thicknessBessironRIA", thickness_bessiron_ria)
        self.__setitem__("resetDeltaTBessiron", reset_delta_t_bessiron)

        corr = OpenFOAMDict(
            items={"singlePhase": single_phase_model, "filmBoiling": film_boiling_model},
            name="correlations",
        )
        self.__setitem__("correlations", corr)

        t0_rows = " ".join(
            f"( {' '.join(str(v) for v in row)} )" for row in T0_data
        )
        axial = OpenFOAMDict(name="axialProfileDict")
        axial["timePoints"] = time_points
        axial["axialLocations"] = axial_locations
        axial["T0Data"] = f"( {t0_rows} )"
        self.__setitem__("axialProfileDict", axial)
