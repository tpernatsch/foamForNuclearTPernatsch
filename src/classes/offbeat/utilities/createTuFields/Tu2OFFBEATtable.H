// Define useful dimensionSet
dimensionSet dimKappa(dimPower/dimLength/dimTemperature);
dimensionSet dimCp(dimEnergy/dimMass/dimTemperature);
dimensionSet dimPower(dimEnergy/dimTime);
dimensionSet dimBurnup(dimEnergy/dimMass);
dimensionSet dimDensity(dimMass/dimVolume);

// Define Hash table for conversion Tu->OFFBEAT
HashTable<TuField, word> Tu2OFFBEATtable
{
    // Stress quantities
    {"Sigma_Radial_MPa",        TuField("sigma_radial", dimPressure, 1e6)},
    {"Sigma_Tangential_MPa",    TuField("sigma_theta", dimPressure, 1e6)},
    {"Sigma_Axial_MPa",         TuField("sigma_axial", dimPressure, 1e6)},
    {"Sigma_Eq_MPa",            TuField("sigmaEq", dimPressure, 1e6)},
    {"Strain_radial",           TuField("epsilon_radial", dimPressure, 1e6)},
    {"Strain_tangential",       TuField("epsilon_theta", dimPressure, 1e6)},
    {"Strain_axial",            TuField("epsilon_axial", dimPressure, 1e6)},
    {"HydrostaticStress_MPa",   TuField("sigmaHyd", dimPressure, 1e6)},

    // Strain quantities
    {"SwellingStrain_radial",       TuField("epsilonSwelling_radial", dimless, 1.0)},
    {"SwellingStrain_tangential",   TuField("epsilonSwelling_theta", dimless, 1.0)},
    {"SwellingStrain_axial",        TuField("epsilonSwelling_axial", dimless, 1.0)},
    {"CreepStrain_radial",          TuField("epsilonCreep_radial", dimless, 1.0)},
    {"CreepStrain_tangential",      TuField("epsilonCreep_theta", dimless, 1.0)},
    {"CreepStrain_axial",           TuField("epsilonCreep_axial", dimless, 1.0)},
    {"RelocationStrain_radial",     TuField("epsilonRelocation_radial", dimless, 1.0)},
    {"RelocationStrain_tangential", TuField("epsilonRelocation_theta", dimless, 1.0)},
    {"RelocationStrain_axial",      TuField("epsilonRelocation_axial", dimless, 1.0)},
    {"SinteringStrain_radial",      TuField("epsilonDensification_radial", dimless, 1.0)},
    {"SinteringStrain_tangential",  TuField("epsilonDensification_theta", dimless, 1.0)},
    {"SinteringStrain_axial",       TuField("epsilonDensification_axial", dimless, 1.0)},
    {"PlasticStrain_radial",        TuField("epsilonP_radial", dimless, 1.0)},
    {"PlasticStrain_tangential",    TuField("epsilonP_theta", dimless, 1.0)},
    {"PlasticStrain_axial",         TuField("epsilonP_axial", dimless, 1.0)},

    // Thermal quantities
    {"Temperature_C",               TuField("T", dimTemperature, 1.0)},
    {"PowerDensity_W_mm3",          TuField("Q", dimPower, 1e9)},
    {"ThermalConductivity_W_mmK",   TuField("k", dimKappa, 1e3)},
    {"SpecificHeat_J_gK",           TuField("Cp", dimCp, 1e3)},
    {"Density_g_mm3",               TuField("rho", dimDensity, 1e6)},
   
    // Burnup
    {"Burnup_MWd_tHM",              TuField("Bu", dimBurnup, 1e-6)}, 
    {"BurnupEffective_MWd_tHM", TuField("SCIANTIX::Effective_burn_up", dimBurnup, 1e-6)},

    // Material
    {"OM_ratio",                    TuField("oxygenMetalRatio", dimless, 1.0)},
    {"grainSize_mm",                TuField("grainRadius", dimLength, 1e-3)},
    {"totalPorosity", TuField("porosity", dimless, 1.0)},
    {"gaseousPorosity", TuField("", dimless, 0.0)},     // ?
    {"sinteringPorosity", TuField("", dimless, 0.0)},   // ?
    {"PorosityHBS", TuField("", dimless, 0.0)},         // ?

    // FGR quantities
    {"FractionalFGR",                       TuField("", dimless, 0.0)}, // ?
    {"FissionGasesCreated_umol_mm3",        TuField("SCIANTIX::Gas_produced", dimless, 1.0)},
    {"FissionGasesGrain_umol_mm3",          TuField("SCIANTIX::Gas_grain", dimless, 1.0)},
    {"FissionGasesGrainSolution_umol_mm3",  TuField("SCIANTIX::Gas_grain_solution", dimless, 1.0)},
    {"FissionGasesGrainBubbles_umol_mm3",   TuField("SCIANTIX::Gas_grain_bubbles", dimless, 1.0)},
    {"IntragranularBubbleConcentration_bbl_mm3", TuField("SCIANTIX::Intragranular_bubble_concentration", dimless, 1.0)},
    {"IntragranularBubbleRadius_nm",    TuField("SCIANTIX::Intragranular_bubble_radius", dimless, 1.0)},
    {"IntergranularBubbleRadius_nm", TuField("SCIANTIX::Intergranular_bubble_radius", dimless, 1.0)},
    {"IntergranularBubbleArea_mm2", TuField("SCIANTIX::Intergranular_bubble_area", dimless, 1.0)},
    {"IntergranularBubbleVolume_mm3", TuField("SCIANTIX::Intergranular_bubble_volume", dimless, 1.0)},
    {"IntergranularBubbleConcentration_bubbl_mm2", TuField("SCIANTIX::Intergranular_bubble_concentration", dimless, 1.0)},
    {"IntergranularAtomsPerBubble_at_bbl", TuField("SCIANTIX::Intergranular_atoms_per_bubble", dimless, 1.0)},
    {"IntergranularVacanciesPerBubble_vac_mm2", TuField("SCIANTIX::Intergranular_vacancies_per_bubble", dimless, 1.0)},
    {"IntergranularFractionalCoverage", TuField("SCIANTIX::Intergranular_fractional_coverage", dimless, 1.0)},
    {"IntergranularSaturationFractionalCoverage", TuField("SCIANTIX::Intergranular_saturation_fractional_coverage", dimless, 1.0)},
    {"IntergranularFractionalIntactness", TuField("SCIANTIX::Intergranular_fractional_intactness", dimless, 1.0)},
    {"ModesInitialConditions_1",  TuField("SCIANTIX::grainModes_1", dimless, 1)},
    {"ModesInitialConditions_2",  TuField("SCIANTIX::grainModes_2", dimless, 1)},
    {"ModesInitialConditions_3",  TuField("SCIANTIX::grainModes_3", dimless, 1)},
    {"ModesInitialConditions_4",  TuField("SCIANTIX::grainModes_4", dimless, 1)},
    {"ModesInitialConditions_5",  TuField("SCIANTIX::grainModes_5", dimless, 1)},
    {"ModesInitialConditions_6",  TuField("SCIANTIX::grainModes_6", dimless, 1)},
    {"ModesInitialConditions_7",  TuField("SCIANTIX::grainModes_7", dimless, 1)},
    {"ModesInitialConditions_8",  TuField("SCIANTIX::grainModes_8", dimless, 1)},
    {"ModesInitialConditions_9",  TuField("SCIANTIX::grainModes_9", dimless, 1)},
    {"ModesInitialConditions_10", TuField("SCIANTIX::grainModes_10", dimless, 1)},
    {"ModesInitialConditions_11", TuField("SCIANTIX::grainModes_11", dimless, 1)},
    {"ModesInitialConditions_12", TuField("SCIANTIX::grainModes_12", dimless, 1)},
    {"ModesInitialConditions_13", TuField("SCIANTIX::grainModes_13", dimless, 1)},
    {"ModesInitialConditions_14", TuField("SCIANTIX::grainModes_14", dimless, 1)},
    {"ModesInitialConditions_15", TuField("SCIANTIX::grainModes_15", dimless, 1)},
    {"ModesInitialConditions_16", TuField("SCIANTIX::grainModes_16", dimless, 1)},
    {"ModesInitialConditions_17", TuField("SCIANTIX::grainModes_17", dimless, 1)},
    {"ModesInitialConditions_18", TuField("SCIANTIX::grainModes_18", dimless, 1)},
    {"ModesInitialConditions_19", TuField("SCIANTIX::grainModes_19", dimless, 1)},
    {"ModesInitialConditions_20", TuField("SCIANTIX::grainModes_20", dimless, 1)},
    {"ModesInitialConditions_21", TuField("SCIANTIX::grainModes_21", dimless, 1)},
    {"ModesInitialConditions_22", TuField("SCIANTIX::grainModes_22", dimless, 1)},
    {"ModesInitialConditions_23", TuField("SCIANTIX::grainModes_23", dimless, 1)},
    {"ModesInitialConditions_24", TuField("SCIANTIX::grainModes_24", dimless, 1)},
    {"ModesInitialConditions_25", TuField("SCIANTIX::grainModes_25", dimless, 1)},
    {"ModesInitialConditions_26", TuField("SCIANTIX::grainModes_26", dimless, 1)},
    {"ModesInitialConditions_27", TuField("SCIANTIX::grainModes_27", dimless, 1)},
    {"ModesInitialConditions_28", TuField("SCIANTIX::grainModes_28", dimless, 1)},
    {"ModesInitialConditions_29", TuField("SCIANTIX::grainModes_29", dimless, 1)},
    {"ModesInitialConditions_30", TuField("SCIANTIX::grainModes_30", dimless, 1)},
    {"FissionGasesGrainBoundary_umol_mm3", TuField("SCIANTIX::Gas_boundary", dimless, 1e3)},
    {"FissionGasesReleased_umol_mm3", TuField("SCIANTIX::Gas_released", dimless, 1e3)},
    {"FissionGasesHBS_umol_mm3", TuField("", dimless, 0.0)},// ?
    {"HeCreated_umol_mm3", TuField("SCIANTIX::Helium_produced", dimless, 1e3)},
    {"HeGrain_umol_mm3", TuField("SCIANTIX::Helium_grain", dimless, 1e3)},
    {"HeGrainBoundary_umol_mm3", TuField("SCIANTIX::Helium_boundary", dimless, 1e3)},
    {"GaseousSwellingIntra", TuField("intragranularGasSwelling", dimless, 1.0)},
    {"GaseousSwellingInter", TuField("intergranularGasSwelling", dimless, 1.0)},
    {"XeCreated_umol_mm3",                  TuField("", dimless, 0.0)},// ?
    {"XeGrain_umol_mm3",                    TuField("", dimless, 0.0)},// ?
    {"XeGrainBoundary_umol_mm3",            TuField("", dimless, 0.0)},// ?
    {"XeHBS_umol_mm3",                      TuField("", dimless, 0.0)},// ?
    {"KrCreated_umol_mm3",                  TuField("", dimless, 0.0)},// ?
    {"KrGrain_umol_mm3",                    TuField("", dimless, 0.0)},// ?
    {"KrGrainBoundary_umol_mm3",            TuField("", dimless, 0.0)},// ?
    {"KrHBS_umol_mm3",                      TuField("", dimless, 0.0)},// ?
    {"CsCreated_umol_mm3",                  TuField("", dimless, 0.0)},// ?
    {"CsGrain_umol_mm3",                    TuField("", dimless, 0.0)},// ?
    {"CsGrainBoundary_umol_mm3",            TuField("", dimless, 0.0)},// ?
    {"NdCreated_umol_mm3",                  TuField("", dimless, 0.0)},// ?
    {"NdGrain_umol_mm3",                    TuField("", dimless, 0.0)},// ?
    {"NdGrainBoundary_umol_mm3",            TuField("", dimless, 0.0)},// ?
    {"XeKrRatioGrain",                      TuField("", dimless, 0.0)},// ?
    {"XeKrRatioGrainBoundary",              TuField("", dimless, 0.0)},// ?
    {"XeKrHBS",                             TuField("", dimless, 0.0)},// ?
    {"XeKrHBXeKrFuelS",                     TuField("", dimless, 0.0)},// ?
    {"GaseousSwelling",                     TuField("", dimless, 0.0)}, // ?
    {"gFormasXe_1", TuField("", dimless, 0.0)}, // ?
    {"gFormasXe_2", TuField("", dimless, 0.0)}, // ?
    {"gFormasXe_3", TuField("", dimless, 0.0)}, // ?
    {"gFormasXe_4", TuField("", dimless, 0.0)}, // ?
    {"gFormasXe_5", TuField("", dimless, 0.0)}, // ?
    {"gFormasXe_6", TuField("", dimless, 0.0)}, // ?
    {"gFormasKr_1", TuField("", dimless, 0.0)}, // ?
    {"gFormasKr_2", TuField("", dimless, 0.0)}, // ?
    {"gFormasKr_3", TuField("", dimless, 0.0)}, // ?
    {"gFormasKr_4", TuField("", dimless, 0.0)}, // ?
    {"gFormasKr_5", TuField("", dimless, 0.0)}, // ?
    {"gFormasKr_6", TuField("", dimless, 0.0)}, // ?
    {"gFormasCs_1", TuField("", dimless, 0.0)}, // ?
    {"gFormasCs_2", TuField("", dimless, 0.0)}, // ?
    {"gFormasCs_3", TuField("", dimless, 0.0)}, // ?
    {"gFormasCs_4", TuField("", dimless, 0.0)}, // ?
    {"gFormasCs_5", TuField("", dimless, 0.0)}, // ?
    {"gFormasCs_6", TuField("", dimless, 0.0)}, // ?
    {"gFormasNd_1", TuField("", dimless, 0.0)}, // ?
    {"gFormasNd_2", TuField("", dimless, 0.0)}, // ?
    {"gFormasNd_3", TuField("", dimless, 0.0)}, // ?
    {"gFormasNd_4", TuField("", dimless, 0.0)}, // ?
    {"gFormasNd_5", TuField("", dimless, 0.0)}, // ?
    {"gFormasNd_6", TuField("", dimless, 0.0)}, // ?
    {"gFormasXeKr_1", TuField("", dimless, 0.0)}, // ?
    {"gFormasXeKr_2", TuField("", dimless, 0.0)}, // ?
    {"gFormasXeKr_3", TuField("", dimless, 0.0)}, // ?
    {"gFormasXeKr_4", TuField("", dimless, 0.0)}, // ?
    {"gFormasXeKr_5", TuField("", dimless, 0.0)}, // ?
    {"gFormasXeKr_6", TuField("", dimless, 0.0)}, // ?

    // Nuclides fields
    // --- Gd
    {"Utot_at_cm3",  TuField("N_U", dimless/dimVolume, 1.0)},
    {"Gd155_at_cm3", TuField("N_Gd155", dimless/dimVolume, 1.0)},
    {"Gd157_at_cm3", TuField("N_Gd157", dimless/dimVolume, 1.0)},
    // --- Th
    {"Th232_at_cm3", TuField("N_Th232", dimless/dimVolume, 1.0)},
    // --- U
    {"U233_at_cm3", TuField("N_U233", dimless/dimVolume, 1.0)},
    {"U234_at_cm3", TuField("N_U234", dimless/dimVolume, 1.0)},
    {"U235_at_cm3", TuField("N_U235", dimless/dimVolume, 1.0)},
    {"U236_at_cm3", TuField("N_U236", dimless/dimVolume, 1.0)},
    {"U238_at_cm3", TuField("N_U238", dimless/dimVolume, 1.0)},
    // --- Np
    {"Np237_at_cm3", TuField("N_Np237", dimless/dimVolume, 1.0)},
    // --- Pu
    {"PuTot_at_cm3", TuField("N_Pu", dimless, 1.0)},
    {"Pu238_at_cm3", TuField("N_Pu238", dimless/dimVolume, 1.0)},
    {"Pu239_at_cm3", TuField("N_Pu239", dimless/dimVolume, 1.0)},
    {"Pu240_at_cm3", TuField("N_Pu240", dimless/dimVolume, 1.0)},
    {"Pu241_at_cm3", TuField("N_Pu241", dimless/dimVolume, 1.0)},
    {"Pu242_at_cm3", TuField("N_Pu242", dimless/dimVolume, 1.0)},
    // --- Am
    {"Am241_at_cm3", TuField("N_Am241", dimless/dimVolume, 1.0)},
    {"Am243_at_cm3", TuField("N_Am243", dimless/dimVolume, 1.0)},
    // --- Cm
    {"Cm242_at_cm3", TuField("N_Cm242", dimless/dimVolume, 1.0)},
    {"Cm243_at_cm3", TuField("N_Cm243", dimless/dimVolume, 1.0)},
    {"Cm244_at_cm3", TuField("N_Cm244", dimless/dimVolume, 1.0)},
    {"Cm245_at_cm3", TuField("N_Cm245", dimless, 1.0)},
    // --- He
    {"He4Tot_at_cm3", TuField("N_He4", dimless, 1.0)},
    {"He4AlphaDecay_at_cm3", TuField("N_He4AlphaDecay", dimless, 1.0)},
    {"He4NAlpha_at_cm3", TuField("N_He4Alpha", dimless, 1.0)},
    {"He4TernFission_at_cm3", TuField("N_He4TernFission", dimless, 1.0)}
};