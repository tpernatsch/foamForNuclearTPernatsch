
#==============================================================================*
# %% 0) Imports

import foamForNuclear as ffn
from foamForNuclear import (
    mesh, functions, fields, boundaryConditions as bc, solvers, thermo, offbeat_lib
)

from foamForNuclear.offbeat_lib import materials
from mesh import build_mesh


def build_case():

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
        Cp = 500,
        Hf = 0,
        mu = 1.0e-4,
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

    # Fluid schemes
    fluid.fvSchemes.ddtSchemes['default'] = "Euler"
    fluid.fvSchemes.gradSchemes['default'] = "Gauss linear"
    fluid.fvSchemes.divSchemes['default'] = "none"
    fluid.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = "Gauss linear"
    fluid.fvSchemes.divSchemes['div(alphaRhoPhi,he)'] = "Gauss linear"
    fluid.fvSchemes.divSchemes['div((nuEff*dev2(T(grad(U)))))'] = "Gauss linear"
    fluid.fvSchemes.divSchemes['div((nuEff*dev(T(grad(U)))))'] = "Gauss linear"
    fluid.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = "Gauss linear"
    fluid.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = "Gauss linear"
    fluid.fvSchemes.laplacianSchemes['default'] = "Gauss linear uncorrected"
    fluid.fvSchemes.snGradSchemes['default'] = "uncorrected"
    fluid.fvSchemes.fluxRequired['default'] = "no"
    fluid.fvSchemes.fluxRequired['p_rgh'] = ""
    fluid.fvSchemes.interpolationSchemes['default'] = "linear"

    # Fluid fvSolution
    fluidSolution = ffn.numerics.fvSolution()
    fluidSolution.append('"p_rgh|p_rghFinal"', ffn.numerics.fvSolutionSolver(
        solver='GAMG',
        smoother='GaussSeidel',
        tolerance=1e-8,
        relTol=1e-3,
    ))
    fluidSolution.append('"U|UFinal"', ffn.numerics.fvSolutionSolver(
        solver='PBiCG',
        preconditioner="DILU",
        tolerance=1e-8,
        relTol=1e-3,
        minIter=1
    ))
    fluidSolution.append('"h|hFinal"', ffn.numerics.fvSolutionSolver(
        solver='PBiCG',
        preconditioner="DILU",
        tolerance=1e-12,
        relTol=1e-3,
        minIter=1
    ))

    fluid.fvSolution = fluidSolution

    # Solid schemes
    solid.fvSchemes.d2dt2Schemes['default'] = "backward"
    solid.fvSchemes.ddtSchemes['default'] = "backward"
    solid.fvSchemes.gradSchemes['default'] = "Gauss linear"
    solid.fvSchemes.divSchemes['default'] = "Gauss linear"
    solid.fvSchemes.laplacianSchemes['default'] = "Gauss linear uncorrected"
    solid.fvSchemes.snGradSchemes['default'] = "uncorrected"
    solid.fvSchemes.interpolationSchemes['default'] = "linear"
    solid.fvSchemes.fluxRequired['default'] = "true"

    # Solid fvSolution
    solidSolution = ffn.numerics.fvSolution()
    solidSolution.append('T', ffn.numerics.fvSolutionSolver(
        solver='PCG',
        preconditioner='DIC',
        tolerance=1e-10,
        relTol=0.01,
        minIter=1,
        maxIter=100
    ))

    solid.fvSolution = solidSolution

    solid.stressAnalysis.nCorrectors = 5
    solid.stressAnalysis.maxOuterIter = 1
    solid.stressAnalysis.relT = 1e-6    

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
        fieldName="p_rgh",
        relTol=1e-6,
        tolerance=1e-6,
    )

    fluid.add_relaxation_on_equation('"U|UFinal"', 0.7)
    fluid.add_relaxation_on_field('"p_rgh|p_rghFinal"', 0.3)

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
        maxResidual=5e-5,
        maxIterations=50,
        fluidRegionName=fluid.region,
        solidRegionName=solid.region,
        fluidPatches=["interface"],
        solidPatches=["top"],
        useHTC=False
    )

    case.coupling = ffn.coupling.Coupling(solvers=[chtLoop])

    # %% 7) Case settings

    case.settings.application = 'GeN-Foam'
    case.settings.endTime = 100
    case.settings.deltaT = 0.1
    case.settings.writeControl = 'runTime'
    case.settings.writeInterval = 10
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

    # %% 9) Return
    return case, fluid_mesh, solid_mesh