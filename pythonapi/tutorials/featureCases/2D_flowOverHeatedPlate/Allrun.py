"""

"""
#==============================================================================*
# %% 0) Imports

import foamForNuclear as ffn
from foamForNuclear import (
    mesh, functions, fields, boundaryConditions as bc, solvers, thermo, offbeat_lib
)

from foamForNuclear.offbeat_lib import materials
from mesh import build_mesh

# %% 1) Mesh
fluid_mesh, solid_mesh = build_mesh()

# %% 2) Case and solvers definition
case = ffn.Case('case/')

fluid = case.add_solver(solvers.thermal_hydraulics.OnePhase(mesh=fluid_mesh))
solid = case.add_solver(solvers.offbeat.Offbeat(mesh=solid_mesh))

# %% 3) Fields
T_fluid = fields.Temperature(internalField=300)
U_fluid = fields.U(internalField=[1, 0, 0])
p_fluid = fields.p_rgh(internalField=0.0)
T_solid = fields.Temperature(internalField=310)

fluid.fields += [T_fluid, U_fluid, p_fluid]
solid.fields += [T_solid]

T_fluid.boundaryField = {
    "interface": bc.Mixed(
        value=T_fluid.internalField,
        refValue=T_fluid.internalField,
        refGradient=0,
        valueFraction=0
    ),
    "inlet": bc.FixedValue(value=T_fluid.internalField),
}

U_fluid.boundaryField = {
    "interface|bottom": bc.FixedValue(value=[0, 0, 0]),
    "inlet": bc.FixedValue(value=U_fluid.internalField),
}

p_fluid.boundaryField = {
    "outlet": bc.FixedValue(value=0.0),
}

T_solid.boundaryField = {
    "top": bc.Mixed(
        value=T_solid.internalField,
        refValue=T_solid.internalField,
        refGradient=0,
        valueFraction=0
    ),
    "bottom": bc.FixedValue(value=T_solid.internalField),
}

# %% 4) Solvers settings

fluid.g = [0, 0, 0]
fluid.turbulenceProperties.simulationType = "laminar"
fluid.thermophysicalProperties = thermo.SodiumBoussinesq(
    molWeight = 1,
    rho0 = 1,
    T0 = 303.0,
    beta = 0,
    Cp = 2000,
    Hf = 0,
    mu = 2.0e-4,
    Pr = 0.01,
)

solid.thermalSolver = offbeat_lib.thermal_solver.SolidConduction()
solid.mechanicsSolver = offbeat_lib.mechanics_solver.Constant()
solid.materials += [materials.Constant(
        name="solid",
        density=1,
        heatCapacity=100,
        conductivity=100,
        emissivity=0,
        YoungModulus=2e+11,
        PoissonRatio=0.3,
        thermalExpansion=1e-5,
        Tref=0
)]

# %% 5) Numerics

# Fluid outer corrector
fluid.pimpleOptions.nOuterCorrectors = 1
fluid.pimpleOptions.nCorrectors = 3
fluid.pimpleOptions.nNonOrthogonalCorrectors = 1
fluid.pimpleOptions.solveEnergy = True
fluid.pimpleOptions.solveFluidMechanics = True
fluid.pimpleOptions.pMin = -1e6
fluid.pimpleOptions.pRefCell = 0
fluid.pimpleOptions.pRefValue = 0

# TODO: I don't like this two commands
fluid.pimpleOptions.add_residual_control_on_field(
    fieldName="U",
    relTol=1e-6,
    tolerance=1e-6,
)

fluid.add_relaxation_on_equation('"U|UFinal"', 0.9)

# Solid outer corrector
solid.stressAnalysis.nCorrectors = 5
solid.stressAnalysis.maxOuterIter = 1
solid.stressAnalysis.relT = 1e-6

# %% 6) Solver coupling
# TODO: I am no sure I understand how this works, and if I like it. 
# Also I am not sure if coupling belongs here and not at the beginning.

chtLoop = ffn.coupling.CHTLoop(
    region="Level_1",
    solvers=[fluid, solid],
    maxResidual=1e-6,
    maxIterations=100,
    fluidRegionName=fluid.region,
    solidRegionName=solid.region,
    fluidPatches=["interface"],
    solidPatches=["top"],
    useHTC=False
)

case.coupling = ffn.coupling.Coupling(solvers=[chtLoop])


# %% 7) Case settings

case.settings.application = 'GeN-Foam'
case.settings.endTime = 8
case.settings.deltaT = 0.01
case.settings.writeControl = 'runTime'
case.settings.writeInterval = 0.5
case.settings.runTimeModifiable = True

# %% 8) Function objects

case.add_function_object(ffn.functions.Probes(
    name="probes",
    fields=["T"],
    enabled=True,
    writeControl="writeTime",
    writeInterval=1,
    region=fluid.region,
    probeLocations=[
        (0.01, 0, 0),
        (0.06, 0, 0),
        (0.11, 0, 0),
        (0.16, 0, 0),
        (0.21, 0, 0),
        (0.26, 0, 0),
        (0.31, 0, 0),
        (0.36, 0, 0),
        (0.41, 0, 0),
        (0.46, 0, 0),
        (0.51, 0, 0),
        (0.56, 0, 0),
        (0.61, 0, 0),
        (0.66, 0, 0),
        (0.71, 0, 0),
        (0.76, 0, 0),
        (0.81, 0, 0),
        (0.86, 0, 0),
        (0.91, 0, 0),
        (0.96, 0, 0),
    ]
))

# %% 9) Preprocessing

print(case)
ffn.allclean()

# Export to OpenFOAM
case.export_to_openfoam()

# coupling.plot_solving_flowchart()
# coupling.plot_solving_graph()

ffn.run_preprocessing(model)

# model.plot_mesh(region=[fluid_mesh, solid_mesh], show_edges=True, normal="z")
# model.plot_mesh(region=fluid_mesh, show_edges=True, normal="z")
# model.plot_mesh(region=solid_mesh, show_edges=True, normal="z")


# %% 10) Run

# ffn.run(model, is_preprocessing=False)


# #==============================================================================*
# # Post-processing

# model.plot_residuals(
#     parameters=["h", "p_fluid"],
#     title="Thermal-hydraulics"
# )

# for fieldName, unit in [('T', 'K'), ("U", "m/s")]:
#     model.plot_slice(
#         region=fluid_mesh,
#         time=settings.endTime,
#         fieldName=fieldName,
#         cmap='RdBu_r',
#         unit=unit,
#         normal="z"
#     )
#     model.plot_animation(
#         region=fluid_mesh,
#         fieldName=fieldName,
#         cmap="RdBu_r",
#         unit=unit,
#         fps=2,
#         normal="z"
#     )

# model.plot_slice(
#     region=[solid_mesh, fluid_mesh],
#     time=settings.endTime,
#     fieldName='T',
#     cmap='RdBu_r',
#     unit='K',
#     normal="z"
# )
# model.plot_animation(
#     region=solid_mesh,
#     fieldName='T',
#     cmap="RdBu_r",
#     unit='K',
#     fps=2,
#     normal="z"
# )


# data, locations = probesFunctionObject.read_from_case(startTime=0, fieldName="T")

# xRange = [loc.x for loc in locations]
# lastTime = model.settings.endTime

# # Normalize temperatures
# T_norm = (data[lastTime] - 300) / 10

# # Get expected values for comparison
# # numerical = pd.read_csv('benchmark/expected_numerical.csv', header=None, names=['x', 'y'])
# # analytical = pd.read_csv('benchmark/expected_analytical.csv', header=None, names=['x', 'y'])

# # Plot
# plt.figure(figsize=(8, 4))
# plt.plot(xRange, T_norm, marker='o', label='foamForNuclear')
# # plt.plot(numerical['x'], numerical['y'], label='Expected Numerical', linestyle='--')
# # plt.plot(analytical['x'], analytical['y'], label='Expected Analytical', linestyle='-.')
# plt.xlabel('Probe Index')
# plt.ylabel(r'Normalized Temperature $\frac{T - 300}{T_s - 300}$')
# plt.title(f'Normalized Temperature at t = {lastTime}s')
# plt.xlim((0, 1))
# plt.grid(True)
# plt.tight_layout()
# plt.legend()
# plt.savefig("fig_results_temperatureComparison.png")


#==============================================================================*
