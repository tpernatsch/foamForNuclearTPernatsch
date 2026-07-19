"""
Fuchs Experiment - Reactivity Insertion Transient
"""
#==============================================================================*
# Imports

import numpy as np
import matplotlib
import tabulate
matplotlib.use('Agg')  # Non-interactive backend to avoid threading issues
import matplotlib.pyplot as plt
import scipy.integrate

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# Constants

innerRadius = 0.00635 / 2  # Inner radius [m]
outerRadius = 0.03645 / 2  # Outer radius [m]
fuelHeight = 0.38099  # Height [m]
nr, nz = 1, 1  # Mesh cells
wedgeAngle = 2  # degrees
numberOfElements = 59  # Number of fuel elements

initialPower = 10  # W

scalingPower = numberOfElements / (wedgeAngle / 360)


#==============================================================================*
# Functions

def fuchs_hansen_solution(t, a0, b, initialPower):
    """
    Fuchs-Hansen analytical solution for power pulse
    """
    c = np.sqrt(a0**2 + 2 * b * initialPower)
    A = (a0 + c) / (c - a0)
    exp_term = np.exp(-c * t)
    return 2 * c**2 * A * exp_term / (b * (A * exp_term + 1)**2)


def calculate_metrics(
        t_sim, P_sim, promptGenerationTime, beta_total, reactivityInsertion, gamma, initialPowerScaled,
        wedgeAngle
    ):
    """
    Calculate FWHM, peak power, and energy metrics
    """
    time_array = np.array(t_sim)
    power_array = np.array(P_sim)
    mask_transient = (time_array >= tStartReactitivity) & (time_array <= tEndPlot)
    tSimulatedReduced = time_array[mask_transient] # - tStartReactitivity
    powerSimulatedReduced = power_array[mask_transient]

    # Fuchs-Hansen parameters
    rho_prime = reactivityInsertion - beta_total
    a0 = rho_prime / promptGenerationTime
    b = gamma / promptGenerationTime

    # Analytical solution (scaled to MW)
    P_analytical = fuchs_hansen_solution(tSimulatedReduced - tStartReactitivity, a0, b, initialPowerScaled) * 1e-6
    P_analytical_max = (a0**2 / (2 * b)) * 1e-6
    FWHM_analytical = 3.525 / a0
    E_analytical_max = a0 / b * 1e-6

    # Simulation metrics
    P_sim_max = powerSimulatedReduced.max()
    idx_peak = np.argmax(powerSimulatedReduced)
    above_half_max = powerSimulatedReduced >= P_sim_max / 2
    crossing_indices = np.where(above_half_max)[0]
    if len(crossing_indices) > 1:
        FWHM_sim = tSimulatedReduced[crossing_indices[-1]] - tSimulatedReduced[crossing_indices[0]]
    else:
        FWHM_sim = None
    # Energy from start to peak
    E_sim = scipy.integrate.trapezoid(powerSimulatedReduced[:idx_peak+1], tSimulatedReduced[:idx_peak+1])

    # Align peaks
    t_peak_sim = tSimulatedReduced[idx_peak]
    t_peak_analytical = tSimulatedReduced[np.argmax(P_analytical)]
    time_shift = t_peak_analytical - t_peak_sim
    tAnalyticalShifted = tSimulatedReduced - time_shift

    # Calculate relative errors
    FWHM_rel_error = (FWHM_sim - FWHM_analytical) / FWHM_analytical * 100 if FWHM_sim else None
    E_rel_error = (E_sim - E_analytical_max) / E_analytical_max * 100
    Pmax_rel_error = (P_sim_max - P_analytical_max) / P_analytical_max * 100

    return {
        't_simulated': tSimulatedReduced,
        'P_simulated': powerSimulatedReduced,
        't_analytical': tAnalyticalShifted,
        'P_analytical': P_analytical,
        'Pmax_sim': P_sim_max,
        'Pmax_analytical': P_analytical_max,
        'FWHM_sim': FWHM_sim,
        'FWHM_analytical': FWHM_analytical,
        'E_sim': E_sim,
        'E_analytical': E_analytical_max,
        'FWHM_rel_error': FWHM_rel_error,
        'E_rel_error': E_rel_error,
        'Pmax_rel_error': Pmax_rel_error
    }


def createMesh(region: str):
    newMesh = mesh.BlockMesh(region=region)
    zone0 = newMesh.create_wedge("zone0", innerRadius, outerRadius, 0, fuelHeight, wedgeAngle, nr=nr, nz=nz)

    wedgeFaces = [
        ("front", zone0.frontFace()),
        ("back", zone0.backFace()),
    ]

    for facename, subface in wedgeFaces:
        face = ffn.mesh.Face(facename, boundaryType="wedge")
        face.add_sub_face(subface)
        newMesh.add_boundary(face)

    outer= ffn.mesh.Face("outer", boundaryType="wall")
    outer.add_sub_face(zone0.rightFace())
    newMesh.add_boundary(outer)

    inner = ffn.mesh.Face("inner", boundaryType="wall")
    inner.add_sub_face(zone0.leftFace())
    newMesh.add_boundary(inner)

    top = ffn.mesh.Face("top", boundaryType="wall")
    top.add_sub_face(zone0.topFace())
    newMesh.add_boundary(top)

    bottom = ffn.mesh.Face("bottom", boundaryType="wall")
    bottom.add_sub_face(zone0.bottomFace())
    newMesh.add_boundary(bottom)

    return(newMesh)


#==============================================================================*
# Meshes

nMesh = createMesh(region="neutroRegion")
tmMesh = createMesh(region="thermoMechanicalRegion")


#==============================================================================*
# Fields

timeFolder0 = ffn.TimeFolder(time=0)

# Neutronics
defaultFlux = ffn.fields.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.fields.Dimension(default='flux')
defaultFlux.internalField = 1e21
defaultFlux.set_boundary_condition('(front|back)', bc.Wedge())
defaultFlux.set_boundary_condition('inner', bc.ZeroGradient())
defaultFlux.set_boundary_condition('outer', bc.FixedValue(0))
defaultFlux.set_boundary_condition('top', bc.FixedValue(0))
defaultFlux.set_boundary_condition('bottom', bc.FixedValue(0))

# Thermomechanics
T_coolant = 293  # Coolant temperature [K]
htc = 70  # Heat transfer coefficient [W/m²/K]

Tmech = ffn.fields.Field("T", region=tmMesh.region)
Tmech.dimensions = ffn.fields.Dimension(default="T")
Tmech.internalField = T_coolant
Tmech.set_boundary_condition('(front|back)', bc.Wedge())
Tmech.set_boundary_condition('inner', bc.ZeroGradient())
Tmech.set_boundary_condition('outer', bc.CustomPatch(
    parameters={
        'type': 'externalWallHeatFluxTemperature',
        'mode': 'coefficient',
        'Ta': f'constant {T_coolant}',
        'h': f'constant {htc}',
        'kappaMethod': 'lookup',
        'kappa': 'k'
    },
    value=T_coolant
))
Tmech.set_boundary_condition('bottom', bc.ZeroGradient())
Tmech.set_boundary_condition('top', bc.ZeroGradient())

timeFolder0.append(defaultFlux)
timeFolder0.append(Tmech)


#==============================================================================*
# Neutronics solver

initialPowerScaled = initialPower / scalingPower

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    solver="diffusionNeutronics",
    mesh=nMesh,
    region=nMesh.region,
    isMeshDeformation=False,
    power=initialPowerScaled,
    keff=1,
)

# neutronicsSolver.timeFolder = timeFolder0

neutronicsSolver.fvSchemes.ddtSchemes['default'] = 'steadyState'
neutronicsSolver.fvSchemes.laplacianSchemes["div(facePhi_,angularFlux_)"] = "Gauss upwind"

neutronicsSolver.fvSolution.solvers["\"flux.*\""] = ffn.numerics.fvSolutionSolver(
    solver="PCG",
    preconditioner="DIC",
    tolerance=1e-6,
    relTol=1e-4
)
neutronicsSolver.fvSolution.solvers["\"precStar.*\""] = neutronicsSolver.fvSolution.solvers["\"prec.*\""]

nuclearData = neutronicsSolver.nuclearData
nuclearData.axialOrientation = ffn.common.Vector(1, 1, 1)

refState = ffn.nuclearData.NuclearDataState(
     name="reference",
     zones=[
         ffn.nuclearData.NuclearDataZone(
            "zone0",
            fuelFraction=1,
            removalXS=[45.1082],
            nuFissionXS=[2554.13],
            powerXS=[1],
            scatteringMatrixP0=[[0.1]],
            discFactor=[1],
            chiPrompt=[1],
            chiDelayed=[1],
            inverseVelocity=[7.1e-08],
            diffusionCoefficient=[0.1275],
            integralFlux=[1],
            decayConstant=[
                0.0124, 0.0305, 0.1115,
                0.3010, 1.1380, 3.0100
            ],
            delayedFraction=[
                0.00023097, 0.00153278, 0.00137180,
                0.00276451, 0.00080489,  0.00029396
            ],
        )
    ]
)

nuclearData.add_state(refState)

#==============================================================================*
# Thermomechanics solver

tmSolver = ffn.solvers.Offbeat(
    region=tmMesh.region,
    solver="extendedThermoMechanics",
    mesh=tmMesh,
    thermalSolver=ffn.offbeat_lib.thermal_solver.SolidConduction(
        heatFluxSummary=True,
        calculateEnthalpy=True
    ),
    mechanicsSolver=ffn.offbeat_lib.mechanics_solver.MechanicsSubSolver(),
    couplingOptions=ffn.offbeat_lib.ThermoMechanicsCouplingOptions(
        correctTFromTH=False,
        correctDispForNeutro=False
    ),
    isSetFvSchemesToDefault=False,
    isSetFvSolutionToDefault=False
)

tmSolver.globalOptions.pinDirection = [0, 0, 1]

fuelMat = ffn.offbeat_lib.materials.Constant(
    name="zone0",
    density=6120,
    heatCapacity=341.6,
    conductivity=8,
    emissivity=0,
    YoungModulus=1e9,
    PoissonRatio=0.3,
    thermalExpansion=8.93e-06,
    Tref=T_coolant
)

tmSolver.add_material(fuelMat)
tmSolver.add_relaxation_on_field('T', 0.8)

tmSolver.fvSchemes.d2dt2Schemes['default'] = "backward"
tmSolver.fvSchemes.ddtSchemes['default'] = "backward"
tmSolver.fvSchemes.gradSchemes['default'] = "Gauss linear"
tmSolver.fvSchemes.divSchemes['default'] = "Gauss linear"
tmSolver.fvSchemes.laplacianSchemes['default'] = "Gauss linear uncorrected"
tmSolver.fvSchemes.snGradSchemes['default'] = "uncorrected"
tmSolver.fvSchemes.interpolationSchemes['default'] = "linear"
tmSolver.fvSchemes.fluxRequired['default'] = "true"

solidSolution = ffn.numerics.fvSolution()
solidSolution.append('T', ffn.numerics.fvSolutionSolver(
    solver='PCG',
    preconditioner='DIC',
    tolerance=1e-10,
    relTol=0.0001,
    minIter=2,
    maxIter=500
))
tmSolver.fvSolution = solidSolution


#==============================================================================*
# Solvers

solvers = ffn.solvers.Solvers([neutronicsSolver, tmSolver])


#==============================================================================*
# Coupling

coupling = ffn.coupling.Coupling(solvers=solvers)

coupling.add_field_transfer(neutronicsSolver, tmSolver, "powerDensity", "Q")
coupling.add_field_transfer(tmSolver, neutronicsSolver, "T", "TFuel")


#==============================================================================*
# Settings

model = ffn.Case(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0],
    caseFolder='steadyState'
)

settings = model.settings

settings.application = 'GeN-Foam'
settings.startFrom = "latestTime"
settings.endTime = 1
settings.deltaT = 1e-3
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 1
settings.runTimeModifiable = True
settings.adjustTimeStep = True

# Add idxField to visualize cellZones, not used in calculations
idxField = neutronicsSolver.create_zone_field(nMesh.cellZones)
timeFolder0.append(idxField)

#==============================================================================*
# Transient configuration (outside if __name__ for import)

tStartReactitivity = 2  # Time when reactivity insertion starts
tEndPlot = tStartReactitivity + 0.3  # Upper time limit for plotting

pointKineticsData = ffn.nuclearData.PointKineticsData(
    promptGenerationTime=4e-05,
    decayConstants=[0.0124, 0.0305, 0.1115, 0.3010, 1.1380, 3.0100],
    delayedFractions=[0.00023097, 0.00153278, 0.00137180, 0.00276451, 0.00080489, 0.00029396],
    feedbackCoeffTFuel=-1e-04
)

beta_total = sum(pointKineticsData.delayedFractions)
reactivityInsertion = 2 * beta_total  # Total reactivity insertion
promptGenerationTime = pointKineticsData.promptGenerationTime
volumeFuel = np.pi * (outerRadius**2 - innerRadius**2) * fuelHeight
massFuel = fuelMat.density.value * volumeFuel * numberOfElements
gamma = np.abs(pointKineticsData.feedbackCoeffTFuel / (fuelMat.heatCapacity.value * massFuel))

pointKineticsData.externalReactivityTimeProfile = ffn.TimeProfile(
    'table',
    startTime=tStartReactitivity,
    table=[
        (0.000, 0.00000),
        (0.025, 0.00120),
        (0.050, 0.00520),
        (0.070, 0.01020),
        (0.080, reactivityInsertion)
    ]
)

#==============================================================================*
# Run

if __name__ == "__main__":
    print(model)

    # Export to OpenFOAM
    model.export_to_openfoam()

    coupling.plot_coupling_graph()
    coupling.plot_solving_graph()
    coupling.plot_solving_flowchart()

    ffn.run(model, is_preprocessing=True)

    model.plot_mesh(region=nMesh, show_edges=True)

    # Post-processing
    lastTime = model.get_time_steps()[-1]
    print(f"keff = {model.keff(time=lastTime)}")

    model.plot_mesh(
        region=nMesh,
        time=lastTime,
        fieldName='TFuel',
        cmap='RdBu_r',
        unit='K'
    )
    model.plot_mesh(
        region=nMesh,
        time=lastTime,
        fieldName='flux0',
        cmap='Blues_r',
        unit='neutron/m2/s'
    )


#==============================================================================*
# Restart for transient

def run_transient(model):
    newFolderName = 'transient'

    ffn.duplicateFolder(model.caseFolder, newFolderName)

    model.caseFolder = newFolderName

    model.settings.endTime = 45
    model.settings.maxPowerVariation = 0.01
    model.settings.maxDeltaT = 3e-3

    neutronicsSolver = [solver for solver in model.solvers if isinstance(solver, ffn.solvers.NeutronicsSolver)][0]

    neutronicsSolver.nuclearData = pointKineticsData
    neutronicsSolver.solver = "pointKinetics"
    neutronicsSolver.eigenvalueNeutronics = False
    neutronicsSolver.fastNeutrons = True

    model.export_to_openfoam()

    ffn.run(model)

    return(model)


def getMaxPowerRelErr(model):
    """
    Used for integration testing
    """
    res = model.get_parameters_from_point_kinetics()

    time = res['time']
    power = res['totalPower']

    powerScaled = [p * 1e-6 * scalingPower for p in power]

    metrics = calculate_metrics(
        time,
        powerScaled,
        promptGenerationTime, beta_total, reactivityInsertion, gamma, initialPowerScaled,
        wedgeAngle
    )

    return(metrics['Pmax_rel_error'])


if __name__ == "__main__":
    model = run_transient(model)

    # Extract transient data

    res = model.get_parameters_from_point_kinetics()

    time = res['time']
    power = res['totalPower']
    TFuel = res["TFuel"]

    powerScaled = [p * 1e-6 * scalingPower for p in power]

    metrics = calculate_metrics(
        time,
        powerScaled,
        promptGenerationTime, beta_total, reactivityInsertion, gamma, initialPowerScaled,
        wedgeAngle
    )

    # Print relative errors only
    print(tabulate.tabulate(
        [
            ['FWHM [s]', metrics['FWHM_sim'], metrics['FWHM_analytical'], f"{metrics['FWHM_rel_error']:.2f}%"],
            ['Energy [MJ]', metrics['E_sim'], metrics['E_analytical'], f"{metrics['E_rel_error']:.2f}%"],
            ['Power max [MW]', metrics['Pmax_sim'], metrics['Pmax_analytical'], f"{metrics['Pmax_rel_error']:.2f}%"],
        ],
        headers=['', 'Simulated', 'Analytical', 'Rel. err [%]'],
        tablefmt='pipe'
    ))


    # Plot 1a: Reactivity components (zoomed in)
    #-------------------------------------------

    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    ax.plot(time, res["totalReactivity"], label="Total")
    ax.plot(time, res["extReactivity"], label="External")
    ax.plot(time, res["dopplerReactivity"], label="Doppler")
    ax.plot(time, res["TFuelReactivity"], label="TFuel")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Reactivity [pcm]")
    ax.set_xlim(tStartReactitivity, tEndPlot)
    ax.legend()
    fig.tight_layout()
    fig.savefig("fig_results_reactivity_zoom.png")
    plt.close()


    # Plot 1b: Reactivity components (full time range)
    #-------------------------------------------------

    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    ax.plot(time, res["totalReactivity"], label="Total")
    ax.plot(time, res["extReactivity"], label="External")
    ax.plot(time, res["dopplerReactivity"], label="Doppler")
    ax.plot(time, res["TFuelReactivity"], label="TFuel")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Reactivity [pcm]")
    ax.set_xlim(time[0], time[-1])
    ax.legend()
    fig.tight_layout()
    fig.savefig("fig_results_reactivity_full.png")
    plt.close()


    # Plot 2: Power comparison (whole reactor)
    #-----------------------------------------

    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    ax.plot(
        time,
        powerScaled,
        color='black',
        label='Simulation'
    )
    ax.plot(
        metrics['t_analytical'],
        metrics['P_analytical'],
        color='black',
        ls='--',
        label='Fuchs-Hansen Analytical'
    )
    ax.set_xlabel("Time after insertion [s]")
    ax.set_ylabel("Power [MW]")
    ax.set_xlim(tStartReactitivity, tEndPlot)
    ax.set_ylim(0)
    ax.legend(
        bbox_to_anchor=(0, 1.02, 1, 0.2),
        loc="lower left",
        mode="expand",
        ncol=2
    )
    fig.tight_layout()
    fig.savefig("fig_results_power_comparison_zoom.png")
    plt.close()


    # Plot 3: Fuel temperature evolution during transient
    #----------------------------------------------------

    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    ax.plot(time, TFuel, color='black')
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Fuel Temperature [K]")
    ax.set_xlim(time[0], time[-1])
    fig.tight_layout()
    fig.savefig("fig_results_fuel_temperature.png")
    plt.close()


    # Plot 4: Fuel temperature with power overlay (dual axis)
    #--------------------------------------------------------

    fig, axTfuel = plt.subplots(figsize=(5, 4), dpi=200)
    axTfuel.plot(
        time,
        TFuel,
        color='tab:red',
        ls='--',
        label='Fuel Temperature'
    )
    axTfuel.set_xlabel("Time [s]")
    axTfuel.set_ylabel("Fuel Temperature [K]", color='tab:red')
    axTfuel.tick_params(axis='y', labelcolor='tab:red')
    axTfuel.set_xlim(tStartReactitivity, tEndPlot)

    axPower = axTfuel.twinx()
    axPower.plot(
        time,
        powerScaled,
        color='black',
        label='Power'
    )
    axPower.set_ylabel("Power [MW]")
    axPower.set_ylim(0)

    fig.tight_layout()
    fig.savefig("fig_results_fuel_temperature_power_zoom.png")
    plt.close()


#==============================================================================*
