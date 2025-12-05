#==============================================================================*

import foamForNuclear as ffn

from model import *


#==============================================================================*

# FMU simulator
prefix = "model.root.system1." if False else ""
FMUSimulator = ffn.FMUSimulator(
    name="FMUSimulator",
    pyClassName="Turbomachine",
    pyFileName="Turbomachine",
    mapping=[
        # Modelica to GeN-Foam
        ("to", "gfT_in", prefix + "Tturbine_out"),
        ("to", "gfMassFlowInlet_in", prefix + "mFlowTurbine_out"),
        ("to", "gfpNozzleChamber_in", prefix + "pNozzleChamber_out"),
        ("to", "gfExtReact_in", prefix + "controlSystem.extReactivity"),
        # Modelica from GeN-Foam
        ("from", "gfTnozzleChamber_out", prefix + "Tnozzle_in"),
        ("from", "gfpInlet_out", prefix + "pTurbine_in"),
        ("from", "gfMassFlowOutlet_out", prefix + "mFlowCoreOutlet_in"),
        ("from", "gfPower_out", prefix + "corePower_in")
    ]
)
FMUSimulator.plot_coupling_graph()

# ffn.allclean()

model.caseFolder = "steadyState_masterGenfoam_PID"

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
model.caseFolder = "transient_masterGenfoam_PID"

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
settings.deltaT = 1e-6
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 10
# Not possible to do adaptative time step with FMU !? yes you can but be
# careful, can be unstable
settings.adjustTimeStep = True
settings.runTimeModifiable = True
settings.maxCo = 0.5 # 0.1
settings.maxDeltaT = 1 # 100e-5
settings.maxPowerVariation = 0.001


model.add_time_folder(timeFolder1)

# model.solvers.append(thSolver)
# model.coupling.append(thSolver)

# model.coupling.add_field_transfer(neutronicsSolver, thSolver, "powerDensity", "powerDensityNeutronics")
# model.coupling.add_field_transfer(neutronicsSolver, thSolver, "secondaryPowerDensity", "powerDensityNeutronicsToLiquid")

# model.coupling.add_field_transfer(thSolver, neutronicsSolver, "T", "TCool")
# model.coupling.add_field_transfer(thSolver, neutronicsSolver, "thermo:rho", "rhoCool")
# model.coupling.add_field_transfer(thSolver, neutronicsSolver, "T.fuelAvForNeutronics", "TFuel")
# model.coupling.add_field_transfer(thSolver, neutronicsSolver, "T.cladAvForNeutronics", "TClad")
# model.coupling.add_field_transfer(thSolver, neutronicsSolver, "T.passiveStructure", "TStructMech")

model.add_function_object(FMUSimulator)

model.coupling.plot_coupling_graph()
model.coupling.plot_solving_graph()
model.coupling.plot_solving_flowchart()

model.export_to_openfoam()

print(model)

ffn.run(model, is_preprocessing=True)


#==============================================================================*
# Post-process

tPlot = 400

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
