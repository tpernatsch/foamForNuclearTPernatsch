# case.py
import foamForNuclear as ffn
from foamForNuclear import (
    mesh, functions, fields, boundaryConditions as bc,)


def build_case(region: str = "fluidRegion"):
    # %% 1) Mesh
    th_mesh = mesh.BlockMesh(region=region)

    low_block = th_mesh.create_cube("low", 0, 0, 0, 0.1, 0.01, 0.5, nz=10)
    mid_block = th_mesh.extrude_top([low_block], "mid", dz=1, nz=20)
    top_block = th_mesh.extrude_top([mid_block], "top", dz=0.5, nz=10)

    outlet = mesh.Face("outlet", boundaryType="wall")
    outlet.add_sub_face(top_block.topFace())

    inlet = mesh.Face("inlet", boundaryType="wall")
    inlet.add_sub_face(low_block.bottomFace())

    walls = mesh.Face("walls", boundaryType="empty")
    for block in [low_block, mid_block, top_block]:
        walls.add_sub_face(block.leftFace())
        walls.add_sub_face(block.rightFace())
        walls.add_sub_face(block.frontFace())
        walls.add_sub_face(block.backFace())

    th_mesh.add_boundary(outlet)
    th_mesh.add_boundary(inlet)
    th_mesh.add_boundary(walls)

    # %% 2) Fields
    t_liquid = fields.Tliquid(internalField=652.15, region=region)
    t_liquid.boundaryField = {"inlet": bc.FixedValue(value=t_liquid.internalField)}

    t_structure = fields.Tstructure(internalField=652.15, region=region)
    t_structure.boundaryField = {"inlet": bc.FixedValue(value=t_structure.internalField)}

    t_vapour = fields.Tvapour(internalField=652.15, region=region)

    u_liquid = fields.Uliquid(internalField=[0, 0, 1], region=region)
    u_vapour = fields.Uvapour(internalField=[0, 0, 1], region=region)

    alpha_liquid = fields.alphaLiquid(internalField=1, region=region)
    alpha_liquid.boundaryField = {
        "inlet": bc.Calculated(alpha_liquid.internalField),
        "outlet": bc.Calculated(alpha_liquid.internalField)}

    alpha_vapour = fields.alphaVapour(internalField=0, region=region)

    p = fields.p(internalField=1e5, region=region)
    p.boundaryField = {
        "inlet": bc.Calculated(p.internalField),
        "outlet": bc.Calculated(p.internalField)}

    p_rgh = fields.p_rgh(internalField=1e5, region=region)
    p_rgh.boundaryField = {
        "inlet": bc.FixedValue(value=1.4e5),
        "outlet": bc.FixedValue(value=p_rgh.internalField)}

    # %% 3) Solver (two-phase)
    solver = ffn.solvers.thermal_hydraulics.TwoPhaseThermalHydraulicsSolver(
        region=region,
        mesh=th_mesh
    )

    solver.residualKd = 10

    # Fluid properties (liquid)
    solver.fluid1 = ffn.solvers.thermal_hydraulics.Fluid(
        name="liquid",
        thermophysicalProperties=ffn.thermo.SodiumPolynomial(),
        stateOfMatter="liquid",
        dispersedDiameterModel=ffn.porous_medium.dispersed_diameter.Constant(value=0.01),
    )
    solver.fluid1.turbulenceProperties.simulationType = "laminar"

    # Fluid properties (vapour)
    solver.fluid2 = ffn.solvers.thermal_hydraulics.Fluid(
        name="vapour",
        thermophysicalProperties=ffn.thermo.SodiumVapourPerfectGas(),
        stateOfMatter="gas",
        thermoResidualAlpha=0.1,
        dispersedDiameterModel=ffn.porous_medium.dispersed_diameter.Constant(value=0.01),
    )
    solver.fluid2.turbulenceProperties.simulationType = "laminar"

    # Structure properties
    lowTopStruct = ffn.porous_medium.Structure(zones=["low", "top"], volumeFraction=0.5, Dh=0.01)

    midStruct = ffn.porous_medium.Structure(zones=["mid"], volumeFraction=0.5, Dh=0.005)
    midStruct.powerModel = ffn.porous_medium.power_models.FixedPower(
        volumetricArea=100,
        T=652.15,
        Cp=500,
        rho=7700,
        powerDensity=4.5e8,
    )

    solver.structures += [lowTopStruct, midStruct]

    # Regime maps
    regimeMap = ffn.porous_medium.regime_map.OneParameter(
        name="slugMist",
        parameter="normalized.alpha.vapour",
        regimeBounds=dict(
            slug = [0, 0.85],
            mist = [0.95, 1],
        ),
    )
    solver.regimeMapModels += [regimeMap]

    # Fluid-structure models
    solver.fluid1_structure.dragModels.append(
        ffn.porous_medium.drag.ReynoldsPower(coeff=2, exp=-0.125, zones=["low", "mid", "top"])
    )

    solver.fluid1_structure.heatTransferModels.append(
        ffn.porous_medium.heat_transfer.SuperpositionNucleateBoiling(
            zones=["mid"],
            forcedConvection=ffn.porous_medium.heat_transfer.NusseltReynoldsPrandtlPower(
                const=7.48467, coeff=0.02994, expRe=0.77, expPr=0.77
            ),
            poolBoiling=ffn.porous_medium.heat_transfer.Shah(useExplicitHeatFlux=False),
            flowEnhancementFactor=ffn.porous_medium.heat_transfer.enhancement.CobraTf(),
            suppressionFactor=ffn.porous_medium.heat_transfer.suppression.CobraTf(),
        )
    )

    solver.fluid1_structure.pairGeometryModel = ffn.porous_medium.pair_geometry.PairGeometryModel(
        contactPartitionModel=ffn.porous_medium.pair_geometry.contact_partition.ByRegime(
            regimeMap="slugMist",
            regimes=[
                dict(name="slug", type="constant", value= 1.0),
                dict(name="mist", type="constant", value= 0.1),
            ],
        ),
    )

    # Fluid-fluid models
    solver.fluid_fluid.dragModels.append(ffn.porous_medium.drag.SchillerNaumann())

    solver.fluid_fluid.heatTransferModels.append(
        ffn.porous_medium.heat_transfer.NusseltReynoldsPrandtlPower(
            const=10, coeff=0, expRe=0, expPr=0, phaseName="liquid"
        )
    )
    solver.fluid_fluid.heatTransferModels.append(
        ffn.porous_medium.heat_transfer.NusseltReynoldsPrandtlPower(
            const=10, coeff=0, expRe=0, expPr=0, phaseName="vapour"
        )
    )

    solver.fluid_fluid.twoPhaseDragMultiplierModel = (
        ffn.porous_medium.two_phase_drag_multiplier.Kaiser88(multiplierFluid="liquid")
    )

    solver.fluid_fluid.pairGeometryModel = ffn.porous_medium.pair_geometry.PairGeometryModel(
        dispersionModel=ffn.porous_medium.pair_geometry.dispersion.ByRegime(
            regimeMap="slugMist",
            regimes=[
                dict(name="slug", type="constant", dispersedPhase="vapour"),
                dict(name="mist", type="constant", dispersedPhase="liquid"),
            ],
        ),
        interfacialAreaDensityModel=ffn.porous_medium.pair_geometry.interfacial_area_density.Spherical(),
    )

    solver.fluid_fluid.phaseChangeModel = ffn.porous_medium.phase_change.HeatDriven(
        mode="conductionLimited",
        correctLatentHeat=False,
        latentHeatModel=ffn.porous_medium.phase_change.latent_heat.FinkLeibowitz(adjust=True),
        saturationModel=ffn.porous_medium.phase_change.saturation.BrowningPotter(),
    )

    # %% 3b) Numerics
    solver.fvSchemes.laplacianSchemes["default"] = "Gauss linear uncorrected"
    solver.fvSchemes.divSchemes["default"] = "none"
    solver.fvSchemes.divSchemes["div(phi,alpha)"] = "Gauss vanLeer"
    solver.fvSchemes.divSchemes["div(phir,alpha)"] = "Gauss vanLeer"
    solver.fvSchemes.divSchemes[r"div\(phi.*,U.*\)"] = "Gauss upwind"
    solver.fvSchemes.divSchemes["div(alphaRhoPhi,U)"] = "Gauss upwind"
    solver.fvSchemes.divSchemes["div(alphaRhoPhiNu,U)"] = "Gauss linear"
    solver.fvSchemes.divSchemes["div(alphaRhoPhi,K)"] = "Gauss upwind"
    solver.fvSchemes.divSchemes[r"div\(alphaRhoPhi.*,(h|e).*\)"] = "Gauss upwind"

    th_solution = ffn.numerics.fvSolution()
    th_solution.append(
        '"p_rgh.*"',
        ffn.numerics.fvSolutionSolver(
            solver="GAMG",
            smoother="DIC",
            tolerance=1e-6,
            relTol=0,
        ),
    )
    smoothSolver = ffn.numerics.fvSolutionSolver(
        solver="smoothSolver",
        smoother="symGaussSeidel",
        tolerance=1e-6,
        relTol=0,
        minIter=1,
    )
    th_solution.append('"e.*"', smoothSolver)
    th_solution.append('"h.*"', smoothSolver)
    th_solution.append(
        "alpha",
        ffn.numerics.fvSolutionSolver(
            solver="MULES",
            adjustSubCycles=True,
            alphaMaxCo=0.25,
        ),
    )
    th_solution.append(
        '".*"',
        ffn.numerics.fvSolutionSolver(
            solver="PBiCGStab",
            preconditioner="diagonal",
            tolerance=1e-6,
            relTol=0.001,
        ),
    )
    solver.fvSolution = th_solution

    solver.pimpleOptions.nCorrectors = 12
    solver.pimpleOptions.nOuterCorrectors = 24
    solver.pimpleOptions.partialEliminationMode = "implicit"
    solver.pimpleOptions.momentumMode = "faceCentered"
    solver.pimpleOptions.oscillationLimiterFraction = 0
    solver.pimpleOptions.correctUntilConvergence = True
    solver.pimpleOptions.massTransferSafetyFactor = 0.1
    solver.pimpleOptions.maxTInterfaceDdt = 10000
    solver.pimpleOptions.enthalpyStabilizationMode = "source"

    solver.pimpleOptions.add_residual_control_on_field(
        fieldName="p_rgh",
        tolerance=1e-5,
        relTol=0,
        useFirstPISOInitialResidual=True,
    )
    solver.pimpleOptions.add_residual_control_on_field(
        fieldName="h.liquid",
        tolerance=1e-5,
        relTol=0,
    )

    solver.add_relaxation_on_equation('".*"', 1)
    solver.add_relaxation_on_field('"alpha.*"', 0.5)
    solver.add_relaxation_on_field('"dmdt.*"', 0.125)
    solver.add_relaxation_on_field("T.interface", 0.125)
    solver.add_relaxation_on_field("T.interfaceFinal", 0.125)
    solver.add_relaxation_on_field("h.vapour", 0.0125)
    solver.add_relaxation_on_field("h.vapourFinal", 0.0125)

    # %% 4) Case
    case = ffn.case.Case(caseFolder="case")
    case.add_solver(solver)
    case.fields  = [
        t_liquid, t_structure, t_vapour,
        u_liquid, u_vapour,
        alpha_liquid, alpha_vapour,
        p, p_rgh
    ]

    # %% 5) Settings (controlDict)
    settings: ffn.control.ControlDict = case.settings
    
    settings.application = "GeN-Foam"
    settings.endTime = 12
    settings.deltaT = 0.001
    settings.writeControl = "adjustableRunTime"
    settings.writeInterval = 0.1
    settings.writePrecision = 6
    settings.timePrecision = 6
    settings.runTimeModifiable = True
    settings.adjustTimeStep = True
    settings.maxDeltaT = 0.01
    settings.maxCo = 2
    settings.maxCoTwoPhase = 0.25
    settings.writeContinuityErrors = True

    return case, th_mesh


if __name__ == "__main__":
    case, th_mesh = build_case()
