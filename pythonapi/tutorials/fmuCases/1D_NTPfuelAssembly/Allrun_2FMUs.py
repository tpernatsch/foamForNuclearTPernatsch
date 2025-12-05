#==============================================================================*

import foamForNuclear as ffn

from model import *


#==============================================================================*

# ffn.allclean()

model.caseFolder = "steadyState_2FMUs"

ffn.duplicateFolder('rootCase', model.caseFolder)

print(model)

model.export_to_openfoam()


if isSimpleMesh:
    ffn.copyFolder(f"./nuclearDataSimple", f"{model.caseFolder}/constant/{nMesh.region}/nuclearData")
else:
    ffn.copyFolder(f"./nuclearData", f"{model.caseFolder}/constant/{nMesh.region}")

ffn.run_preprocessing(model=model)

model.plot_mesh(
    region=nMesh,
    show_edges=True
)

# Run
ffn.run(model)


#==============================================================================*
#==============================================================================*
# Restart for transient

oldFolder = model.caseFolder
model.caseFolder = "transient_2FMUs"

ffn.duplicateFolder(oldFolder, model.caseFolder)

timeFolder1 = ffn.TimeFolder(settings.endTime)

timeFolder1.append(alphat)
timeFolder1.append(epsilon)
timeFolder1.append(k)
timeFolder1.append(nut)
timeFolder1.append(p)
timeFolder1.append(p_rgh)
timeFolder1.append(T)
# timeFolder1.append(Tmatrix)
timeFolder1.append(U)

pointKineticsData = ffn.PointKineticsData(
    fastNeutrons=True,
    promptGenerationTime=3.571485e-05,
    delayedFractions=[
        2.285127e-04,  # (4.402223e-08)
        1.180026e-03,  # (2.272815e-07)
        1.126848e-03,  # (2.170125e-07)
        2.527636e-03,  # (4.866769e-07)
        1.037503e-03,  # (1.996539e-07)
        4.345650e-04,  # (8.363008e-08)
    ],
    decayConstants=[
        1.333623e-02,  # (2.566366e-06)
        3.273731e-02,  # (6.297438e-06)
        1.207837e-01,  # (2.322766e-05)
        3.028201e-01,  # (5.820689e-05)
        8.496638e-01,  # (1.631340e-04)
        2.853584e+00,  # (5.479324e-04)
    ],
    feedbackCoeffDoppler=-1.474e-04, # (+/- 1.771e-04)
    feedbackCoeffTFuel=-8.796e-07, # [1/K] (+/- 3.874e-08)
    feedbackCoeffRhoCool=6.484e-03, # [1/(kg/m3)] (+/- 2.202e-03)
    fuelFeedbackZones=["fuelElement"],
    coolFeedbackZones=["fuelElement"],
    structFeedbackZones=["fuelElement"],
)

pointKineticsData.externalReactivityTimeProfile = ffn.TimeProfile(
    type='fmi',
    nameFromFMU='gfExtReact_in',
    initialValue=0
)
neutronicsSolver.nuclearData = pointKineticsData
neutronicsSolver.solver = "pointKinetics"
neutronicsSolver.eigenvalueNeutronics = False

settings.endTime = 1000
settings.deltaT = 1e-9
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 10
# Not possible to do adaptative time step with FMU !? yes you can but be
# careful, can be unstable
settings.adjustTimeStep = True
settings.runTimeModifiable = True
settings.maxCo = 0.5 # 0.1
settings.maxDeltaT = 1 # 100e-5
settings.maxPowerVariation = 0.001
settings.solveFMI = True

model.add_time_folder(timeFolder1)

model.coupling.plot_coupling_graph()
model.coupling.plot_solving_graph()
model.coupling.plot_solving_flowchart()

model.export_to_openfoam()

ffn.run_preprocessing(model=model)

print(model)

ffn.generateCaseAsFMU(model.caseFolder, 'CoreFMU')


# MasterRunner = ffn.FMPyMasterRunnerParallel(
MasterRunner = ffn.FMPyMasterRunner(
    fmuNames={
        'core': 'CoreFMU.fmu',
        'turbomachinery': 'modelica/NTPTurbomachine.fmu'
    },
    connections=[
        # Modelica to GeN-Foam
        ('turbomachinery', "Tturbine_out",                  'core', "gfT_in"),
        ('turbomachinery', "mFlowTurbine_out",              'core', "gfMassFlowInlet_in"),
        ('turbomachinery', "pNozzleChamber_out",            'core', "gfpNozzleChamber_in"),
        ('turbomachinery', "controlSystem.extReactivity",   'core', "gfExtReact_in"),
        # GeN-Foam to Modelica
        ('core', "gfTnozzleChamber_out",    'turbomachinery', "Tnozzle_in"),
        ('core', "gfpInlet_out",            'turbomachinery', "pTurbine_in"),
        ('core', "gfMassFlowOutlet_out",    'turbomachinery', "mFlowCoreOutlet_in"),
        ('core', "gfPower_out",             'turbomachinery', "corePower_in")
    ],
    outputFilenames={
        'core': 'CoreFMU/CoreFMU.csv',
        'turbomachinery': 'CoreFMU/Turbomachine.csv',
    },
    startTime=10,
    endTime=settings.endTime,
    stepSize=settings.deltaT,
    minStepSize=settings.deltaT/10,
    maxStepSize=1,
    minNImplicitSteps=0,
    maxNImplicitSteps=1, # Higher than 1 causes problem on integrator block!? Wrong propellant mass depletion
    adjustStepSize=True,
    writeInterval=0.1,
    relaxationFactorTimeStep=0.01,
    parametersToRecord={
        'turbomachinery': [
            'pTurbine_in',
            'SinkP1.in_p0',
            'pNozzleChamber_out',
            'sensTnozzleWall_in.T',
            'sensTnozzleWall_out.T',
            'Tturbine_out',
            'Tnozzle_in',
            'sourceMassFlow.in_T',
            'NozzleChamber.T[1]',
            'NozzleChamber.T[11]',
            'sensW.w',
            'corePower_in',
            'mFlowTurbine_out',
            'mFlowCoreOutlet_in',
            'Turbine1.omega',
            'valveLin.cmd',
            'controlSystem.PI_massFlow.u_s',
            'controlSystem.PI_massFlow.u_m',
            'controlSystem.extReactivity',
            'controlSystem.controlDrumReactivity.angle',
            'controlSystem.PI_Power.u_s',
            'controlSystem.PI_Power.u_m',
            'controlSystem.corePower.y',
            'controlSystem.PI_Power.y',
            'controlSystem.PI_DoublingTime.y',
            'controlSystem.PI_DoublingTime.u_m',
            # 'controlSystem.minCDangle.y',
            'controlSystem.PI_Thrust.y',
            'controlSystem.PI_Thrust.u_s',
            'controlSystem.PI_Thrust.u_m',
            'controlSystem.PI_DeltaV.u_s',
            'controlSystem.PI_DeltaV.u_m',
            'controlSystem.thrustCalculator.thrust',
            'controlSystem.thrustCalculator.exhaustVelocity',
            'controlSystem.doublingTimeCalculator.doublingTime',
            'controlSystem.thrustCmd.y',
            'controlSystem.deltaVCalculator.deltaV',
            'controlSystem.deltaVCalculator.endPropellant.activePort',
            'controlSystem.deltaVCalculator.totalMassPropellant.y',
            'controlSystem.deltaVCalculator.totalMassPropellant.u1',
            'controlSystem.deltaVCalculator.totalMassPropellant.u2',
            # 'controlSystem.feedbackThrust.u1',
            'firstOrder.T',    # Tturbine_out
            'firstOrder1.T',     # mFlowTurbine_out
            'firstOrder2.T',     # pTurbine_in
            'firstOrder5.T',     # mFlowCoreOutlet_in
            'firstOrder3.T',     # Tnozzle_in
            'firstOrder4.T',     # pNozzleChamber_out
        ]
    },
    initialParametersBeforeFMUInit={
        'turbomachinery': {
            'firstOrder.T': 0.1,    # 0.5, Tturbine_out
            'firstOrder1.T': 0.1,     # 1, mFlowTurbine_out
            'firstOrder2.T': 0.1,     # 1, pTurbine_in
            'firstOrder5.T': 0.1,     # 1, mFlowCoreOutlet_in
            'firstOrder3.T': 0.1,     # 1, Tnozzle_in
            'firstOrder4.T': 0.1,     # 1, pNozzleChamber_out
        }
    }
)
MasterRunner.simulate()


#==============================================================================*
# Post-process

tPlot = 340

try:
    model.plot_mesh(
        region=nMesh,
        time=tPlot,
        fieldName='powerDensity',
        cmap='RdBu_r',
        unit='W/m3',
        show_edges=True,
        lighting=False
    )
    for field, unit in [('T', 'K'), ('p', 'Pa'), ('U', 'm/s')]:
        model.plot_mesh(
            region=thMesh,
            time=tPlot,
            fieldName=field,
            cmap='RdBu_r',
            unit=unit,
            show_edges=True,
            lighting=False
        )
        model.plot_slice(
            region=thMesh,
            time=tPlot,
            fieldName=field,
            cmap='RdBu_r',
            unit=unit,
            show_edges=True,
            lighting=False
        )
except:
    pass

model.plot_animation(
    region=thMesh,
    fieldName="T",
    cmap='RdBu_r',
    show_edges=True,
    unit='K',
    lighting=False,
    fps=10
)
model.plot_animation(
    region=thMesh,
    fieldName="U",
    cmap='RdBu_r',
    show_edges=True,
    unit='m/s',
    lighting=False,
    fps=10
)
model.plot_animation(
    region=nMesh,
    fieldName="powerDensity",
    cmap='RdBu_r',
    show_edges=True,
    unit='W/m3',
    lighting=False,
    fps=10
)

#==============================================================================*
