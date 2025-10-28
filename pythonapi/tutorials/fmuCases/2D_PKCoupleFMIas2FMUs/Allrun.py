"""

"""
#==============================================================================*
# Imports

from copy import copy
import pandas as pd
import matplotlib.pyplot as plt
import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


def generateSteadyStateModel():
    #==============================================================================*
    # Mesh

    nMesh = mesh.BlockMesh(region="neutroRegion")

    zone0 = nMesh.createCube("zone0", 0, 0, 0, 1, 0.1, 1, nx=25, nz=25)

    inlet = mesh.Face("inlet")
    inlet.addSubFace(zone0.bottomFace())

    outlet = mesh.Face("outlet")
    outlet.addSubFace(zone0.topFace())

    fixedWalls = mesh.Face("fixedWalls", boundaryType="wall")
    fixedWalls.addSubFace(zone0.leftFace())
    fixedWalls.addSubFace(zone0.rightFace())

    frontAndBack = mesh.Face("frontAndBack", boundaryType="empty")
    frontAndBack.addSubFace(zone0.frontFace())
    frontAndBack.addSubFace(zone0.backFace())

    nMesh.addBoundary(inlet)
    nMesh.addBoundary(outlet)
    nMesh.addBoundary(fixedWalls)
    nMesh.addBoundary(frontAndBack)


    thMesh = copy(nMesh)
    thMesh.region = "fluidRegion"


    #==============================================================================*
    # Time folder

    timeFolder0 = ffn.TimeFolder(0)

    defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
    defaultFlux.dimensions = ffn.Dimension(default='neutronFlux')
    defaultFlux.internalField = 1
    defaultFlux.set_boundary_condition(fixedWalls, bc.FixedValue(0))
    defaultFlux.set_boundary_condition(inlet, bc.FixedValue(0))
    defaultFlux.set_boundary_condition(outlet, bc.FixedValue(0))
    defaultFlux.set_boundary_condition(frontAndBack, bc.Empty())

    T = ffn.Field("T", region=thMesh.region)
    T.dimensions = ffn.Dimension(default='T')
    T.internalField = 600
    T.set_boundary_condition(fixedWalls, bc.ZeroGradient())
    T.set_boundary_condition(inlet, bc.FixedValue(600))
    T.set_boundary_condition(outlet, bc.ZeroGradient())
    T.set_boundary_condition(frontAndBack, bc.Empty())

    U = ffn.Field("U", region=thMesh.region)
    U.dimensions = ffn.Dimension(default='U')
    U.internalField = ffn.Vector(0, 0, 1)
    U.set_boundary_condition(fixedWalls, bc.Slip())
    U.set_boundary_condition(inlet, bc.FixedValue(U.internalField))
    U.set_boundary_condition(outlet, bc.ZeroGradient())
    U.set_boundary_condition(frontAndBack, bc.Empty())

    p = ffn.Field("p", region=thMesh.region)
    p.dimensions = ffn.Dimension(default='p')
    p.internalField = 1e5
    p.set_boundary_condition(fixedWalls, bc.Calculated(p.internalField))
    p.set_boundary_condition(inlet, bc.Calculated(p.internalField))
    p.set_boundary_condition(outlet, bc.Calculated(p.internalField))
    p.set_boundary_condition(frontAndBack, bc.Empty())

    p_rgh = ffn.Field("p_rgh", region=thMesh.region)
    p_rgh.dimensions = ffn.Dimension(default='p')
    p_rgh.internalField = 1e5
    p_rgh.set_boundary_condition(fixedWalls, bc.FixedFluxPressure(p_rgh.internalField))
    p_rgh.set_boundary_condition(inlet, bc.FixedFluxPressure(p_rgh.internalField))
    p_rgh.set_boundary_condition(outlet, bc.FixedValue(p_rgh.internalField))
    p_rgh.set_boundary_condition(frontAndBack, bc.Empty())

    timeFolder0.append(defaultFlux)
    timeFolder0.append(T)
    timeFolder0.append(U)
    timeFolder0.append(p)
    timeFolder0.append(p_rgh)


    #==============================================================================*
    # Thermal-hydraulics solver

    thSolver = ffn.ThermalHydraulicsSolver(
        "fluidRegion", "onePhase",
        mesh=thMesh,
        removeBaffles=True,
        isSetFvSolutionToDefault=False
    )

    thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumConst()

    thSolver.turbulenceProperties.simulationType = 'laminar'

    core = ffn.StructureProperty(['zone0'], volumeFraction=0.5, Dh=0.01)

    core.powerModel = ffn.NuclearFuelPin(
        fuelInnerRadius=0.0012,
        fuelOuterRadius=0.004715,
        cladInnerRadius=0.004865,
        cladOuterRadius=0.005365,
        fuelMeshSize=16,
        cladMeshSize=4,
        fuelRho=10480,
        fuelCp=250,
        fuelK=3,
        cladRho=7500,
        cladCp=500,
        cladK=20,
        fuelT=600,
        cladT=600,
        gapH=3000
        # gapHPowerDensity=[
        #     (   0              , 528.4782467867  ),
        #     (   795443280.75916, 3170.8694807204 )
        # ]
    )

    core.add_passive_structure(
        volumetricArea=2,
        rho=7700,
        Cp=500,
        T=600
    )

    thSolver.add_structure_property(core)

    thSolver.add_drag_model(
        ffn.ReynoldsPower(0.687, -0.25, zones=['zone0'])
    )
    thSolver.add_heat_transfer_model(
        ffn.NusseltReynoldsPrandtlPower(4.82, 0.0185, 0.827, 0.827, zones=['zone0'])
    )


    thSolution = ffn.fvSolution()

    thSolution.append('"p_rgh.*"', ffn.fvSolutionSolver(
        solver='GAMG', smoother='DIC', tolerance=1e-8, relTol=0
    ))
    thSolution.append('"e.*"', ffn.fvSolutionSolver(
        solver='smoothSolver',
        smoother='symGaussSeidel',
        tolerance=1e-8, relTol=0, minIter=0
    ))
    thSolution.append('"h.*"', ffn.fvSolutionSolver(
        solver='smoothSolver',
        smoother='symGaussSeidel',
        tolerance=1e-8, relTol=0, minIter=0
    ))
    thSolution.append('".*"', ffn.fvSolutionSolver(
        solver='PBiCGStab',
        preconditioner='diagonal',
        tolerance=1e-6, relTol=0.001
    ))

    thSolver.fvSolution = thSolution

    thSolver.pimpleOptions.nOuterCorrectors = 3


    thSolver.fvSchemes.divSchemes['default'] = 'none'
    thSolver.fvSchemes.divSchemes['div(phi,alpha)'] = 'Gauss vanLeer'
    thSolver.fvSchemes.divSchemes['div(phir,alpha)'] = 'Gauss vanLeer'
    thSolver.fvSchemes.divSchemes['div\(phi.*,U.*\)'] = 'Gauss upwind'
    thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = 'Gauss upwind'
    thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = 'Gauss linear'
    thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = 'Gauss upwind'
    thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,k.*\)'] = 'Gauss upwind'
    thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,epsilon.*\)'] = 'Gauss upwind'
    thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,(h|e).*\)'] = 'Gauss upwind'

    thSolver.fvSchemes.laplacianSchemes['default'] = 'Gauss linear uncorrected'

    thSolver.fvSchemes.snGradSchemes['default'] = 'uncorrected'


    #==============================================================================*
    # Neutronics solver

    neutronicsSolver = ffn.NeutronicsSolver(
        "neutroRegion",
        "diffusionNeutronics",
        mesh=nMesh, power=10000000
    )

    neutronicsSolver.nuclearData.add_state(
        ffn.NuclearDataState(
            name='reference',
            zones=[
                ffn.NuclearDataZone(
                    name='zone0',
                    fuelFraction=1,
                    inverseVelocity=[7.1e-08, 0.0045],
                    diffusionCoefficient=[0.1, 0.04],
                    nuFissionXS=[1, 1.75],
                    powerXS=[1, 1],
                    scatteringMatrixP0=[[10, 5], [0, 10]],
                    removalXS=[5, 4],
                    chiPrompt=[1, 0],
                    chiDelayed=[0.5, 0.5],
                    delayedFraction=[7.2315e-05, 0.000609661, 0.000471181, 0.00118907, 0.000445487, 9.58515e-05],
                    decayConstant=[0.0125371, 0.0300828, 0.109879, 0.325484, 1.3036, 9.51817],
                    discFactor=[1, 1],
                    integralFlux=[1, 1]
                )
            ]
        )
    )

    neutronicsSolver.fvSchemes.divSchemes['div(facePhi_,angularFlux_)'] = "Gauss upwind"

    neutronicsSolver.neutronTransportOptions.integralPredictor = True
    neutronicsSolver.neutronTransportOptions.aitkenAcceleration = True
    neutronicsSolver.neutronTransportOptions.maxNeutronIterations = 50

    #==============================================================================*
    # Solvers

    solvers = ffn.Solvers([neutronicsSolver, thSolver])

    #==============================================================================*
    # Coupling

    coupling = ffn.Coupling(solvers)

    coupling.add_field_transfer(neutronicsSolver, thSolver, 'powerDensity', 'powerDensityNeutronics')
    coupling.add_field_transfer(neutronicsSolver, thSolver, 'secondaryPowerDensity', 'powerDensityNeutronicsToLiquid')
    coupling.add_field_transfer(thSolver, neutronicsSolver, 'T', 'TCool')
    coupling.add_field_transfer(thSolver, neutronicsSolver, 'thermo:rho', 'rhoCool')
    coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.fuelAvForNeutronics', 'TFuel')
    coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.cladAvForNeutronics', 'TClad')
    coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.passiveStructure', 'TStructMech')


    #==============================================================================*
    # FMU
    # The FMU needs to be linked even during steady-state, this improves stability
    # of the results especially if PID controllers are present

    # Coupling interface
    externalCouplingDict = ffn.ExternalCouplingDict()

    fieldIntegralToFMU = ffn.FieldIntegralToFMU(
        name="fieldIntegralToFMU",
        nameFMU="gfPower",
        fieldName="powerDensity",
        region="neutroRegion",
        cellZone="zone0",
        initValue=10e6
    )
    externalCouplingDict.append(fieldIntegralToFMU)



    #==============================================================================*
    # Model

    model = ffn.Model(
        solvers=solvers,
        coupling=coupling,
        timeFolders=[timeFolder0],
        externalCouplingDict=externalCouplingDict,
        caseFolder='steadyState'
    )

    settings = model.settings

    settings.application = "GeN-Foam"
    settings.startTime = 0
    settings.endTime = 100
    settings.deltaT = 1
    settings.startFrom = 'latestTime'
    settings.writeControl = 'adjustableRunTime'
    settings.writeInterval = settings.endTime
    settings.adjustTimeStep = True
    settings.runTimeModifiable = True
    settings.maxDeltaT = 1
    settings.maxCo = 1
    settings.maxPowerVariation = 0.01

    coupling.plot_coupling_graph()
    coupling.plot_solving_graph()
    coupling.plot_solving_flowchart()

    print(model)

    model.export_to_openfoam()

    return(model, neutronicsSolver)


def generateTransientModel(model: ffn.Model, neutronicsSolver: ffn.Solver):
    try:
        ffn.duplicateFolder(model.caseFolder, "transient")
    except:
        pass

    model.caseFolder = "transient"
    model.settings.endTime = 600
    model.settings.deltaT = 1e-3
    # model.settings.writeControl = 'timeStep'
    model.settings.writeInterval = 5
    model.settings.solveFMI = True

    # Change nuclear data
    pointKineticsData = ffn.PointKineticsData(
        fastNeutrons=True,
        promptGenerationTime=1e-06,
        delayedFractions=[
            7.2315e-05,
            0.000609661,
            0.000471181,
            0.00118907,
            0.000445487,
            9.58515e-05,
        ],
        decayConstants=[
            0.0125371,
            0.0300828,
            0.109879,
            0.325484,
            1.3036,
            9.51817,
        ],
        feedbackCoeffTFuel=-3e-06
    )

    pointKineticsData.controlRodReactivityMap = [
        (  0.1, -0.01 ),
        (  0.0,  0.0  ),
        ( -0.1,  0.01 ),
    ]

    pointKineticsData.externalReactivityTimeProfile = ffn.TimeProfile(
        type='fmi',
        nameFromFMU='gfExtReact',
        initialValue=0
    )

    neutronicsSolver.nuclearData = pointKineticsData
    neutronicsSolver.solver = "pointKinetics"
    neutronicsSolver.eigenvalueNeutronics = False


    # Export and overwrite the case
    model.export_to_openfoam()

    print(model)

    return(model)


def runCoupledFMU():
    MasterRunner = ffn.FMPyMasterRunnerParallel(
    # MasterRunner = ffn.FMPyMasterRunner(
        fmuNames={
            "controller": "ExternalReactivityController.fmu",
            "core1": "CoreFMU.fmu",
            "core2": "CoreFMU.fmu",
            "core3": "CoreFMU.fmu",
            "core4": "CoreFMU.fmu",
            "core5": "CoreFMU.fmu",
            "core6": "CoreFMU.fmu",
            "core7": "CoreFMU.fmu",
            "core8": "CoreFMU.fmu",
        },
        connections=[
            ("controller", "externalReactivity", "core1", "gfExtReact"),
            ("controller", "externalReactivity", "core2", "gfExtReact"),
            ("controller", "externalReactivity", "core3", "gfExtReact"),
            ("controller", "externalReactivity", "core4", "gfExtReact"),
            ("controller", "externalReactivity", "core5", "gfExtReact"),
            ("controller", "externalReactivity", "core6", "gfExtReact"),
            ("controller", "externalReactivity", "core7", "gfExtReact"),
            ("controller", "externalReactivity", "core8", "gfExtReact"),
            ("core1", "gfPower", "controller", "power")
        ],
        startTime=100,
        endTime=450,
        stepSize=0.001,
        minStepSize=0.001,
        maxStepSize=1,
        minNImplicitSteps=0,
        maxNImplicitSteps=2,
        adjustStepSize=True,
        writeInterval=1,
        parametersToRecord={
            "core1": ["gfPower", "gfExtReact"],
            "core2": ["gfPower", "gfExtReact"],
            "core3": ["gfPower", "gfExtReact"],
            "core4": ["gfPower", "gfExtReact"],
            "core5": ["gfPower", "gfExtReact"],
            "core6": ["gfPower", "gfExtReact"],
            "core7": ["gfPower", "gfExtReact"],
            "core8": ["gfPower", "gfExtReact"],
            "controller": ["power", "externalReactivity", "pid.u_s"]
        },
        initialParametersBeforeFMUInit={
            "core1": {"port": 8001, 'outputPath': 'CoreFMU1'},
            "core2": {"port": 8002, 'outputPath': 'CoreFMU2'},
            "core3": {"port": 8003, 'outputPath': 'CoreFMU3'},
            "core4": {"port": 8004, 'outputPath': 'CoreFMU4'},
            "core5": {"port": 8005, 'outputPath': 'CoreFMU5'},
            "core6": {"port": 8006, 'outputPath': 'CoreFMU6'},
            "core7": {"port": 8007, 'outputPath': 'CoreFMU7'},
            "core8": {"port": 8008, 'outputPath': 'CoreFMU8'},
        }
    )
    MasterRunner.simulate()



#==============================================================================*
# Main

if __name__ == "__main__":

    # Steady-state calculations
    model, neutronicsSolver = generateSteadyStateModel()

    ffn.run(model=model, is_preprocessing=True)

    # Transient
    model = generateTransientModel(model=model, neutronicsSolver=neutronicsSolver)

    # Export GeN-Foam as an FMU
    ffn.generateCaseAsFMU(model.caseFolder, "CoreFMU")

    runCoupledFMU()



    # Post processing
    model.caseFolder = 'CoreFMU1'

    model.plot_residuals(
        parameters=['fluxStar0'],
        title="Neutronics"
    )
    model.plot_residuals(
        parameters=['p_rgh', 'h'],
        title="Thermal-hydraulics"
    )

    dataCsvCtrl = pd.read_csv("ExternalReactivityController_0.csv")

    res = model.get_parameters_from_point_kinetics()


    fig, ax = plt.subplots()
    # ax.plot(res['time'], [e * 1e-6 for e in res['totalPower']], label="GeN-Foam PK")
    ax.plot(dataCsvCtrl['time'], [e * 1e-6 for e in dataCsvCtrl['pid.u_s']], label="Command")
    for i in range(8):
        dataCsvCore = pd.read_csv(f"CoreFMU_{i+1}.csv")
        ax.plot(dataCsvCore['time'], [e * 1e-6 for e in dataCsvCore['gfPower']], label=f"GeN-Foam integralFMI {i+1}")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Power [MW]")
    ax.legend()
    fig.tight_layout()
    fig.savefig("fig_results_power.png")
    plt.close()


    fig, ax = plt.subplots()
    ax.plot(res['time'], res["totalReactivity"], label="Total")
    ax.plot(res['time'], res["extReactivity"], label="External")
    ax.plot(res['time'], res["dopplerReactivity"], label="Doppler")
    ax.plot(res['time'], res["TFuelReactivity"], label="TFuel")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Reactivity [pcm]")
    ax.legend()
    fig.tight_layout()
    fig.savefig("fig_results_reactivity.png")
    plt.close()

#==============================================================================*
