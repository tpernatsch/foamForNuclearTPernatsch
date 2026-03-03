"""
1D Slab Reactor with Constant External Source - Neutronics Diffusion
"""
#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import numpy as np

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# Nuclear data parameters

H = 0.8  # Slab height [m]
S0 = 3.2e3  # External source [1/(m^3*s)]
D = 0.1275  # Diffusion coefficient [m]
sigma_a = 5.725  # Absorption cross section [m^-1]
nu_sigma_f = 5.6675  # Nu-fission cross section [m^-1]

# Calculate k_inf
k_inf = nu_sigma_f / sigma_a
L_squared = D / sigma_a  # Diffusion area L² = D/Σₐ
B_squared = (np.pi / H)**2  # Geometric buckling B² = (π/H)²
k_eff = k_inf / (1 + L_squared * B_squared)

print(f"\n{'='*60}")
print(f"  k_inf = {k_inf:.6f}")
print(f"  k_eff = {k_eff:.6f}")
print(f"{'='*60}\n")


#==============================================================================*
# Mesh

def createMesh(nz: int):
    nMesh = mesh.BlockMesh(region='neutroRegion')

    zone0 = nMesh.create_cube('zone0', -0.05, -0.05, 0, 0.05, 0.05, H, 1, 1, nz)

    walls = ffn.mesh.Face("walls", boundaryType="wall")
    walls.add_sub_face(zone0.frontFace())
    walls.add_sub_face(zone0.backFace())
    walls.add_sub_face(zone0.leftFace())
    walls.add_sub_face(zone0.rightFace())

    top = ffn.mesh.Face("top")
    top.add_sub_face(zone0.topFace())

    bottom = ffn.mesh.Face("bottom")
    bottom.add_sub_face(zone0.bottomFace())

    nMesh.add_boundary(walls)
    nMesh.add_boundary(top)
    nMesh.add_boundary(bottom)

    return(nMesh)

nz = 101

nMesh = createMesh(nz=nz)


#==============================================================================*
# Fields

timeFolder0 = ffn.timeFolder.TimeFolder(0)

defaultFlux = ffn.fields.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.fields.Dimension(default='flux')
defaultFlux.internalField = 0.5e3
defaultFlux.set_boundary_condition("walls", bc.ZeroGradient())
defaultFlux.set_boundary_condition("top", bc.FixedValue(0))
defaultFlux.set_boundary_condition("bottom", bc.FixedValue(0))

timeFolder0.append(defaultFlux)

defaultExternalSourceFlux = ffn.fields.Field("defaultExternalSourceFlux", region=nMesh.region)
defaultExternalSourceFlux.dimensions = ffn.fields.Dimension(length='-3', time='-1')
defaultExternalSourceFlux.internalField = S0
defaultExternalSourceFlux.set_boundary_condition("walls", bc.FixedValue(0))
defaultExternalSourceFlux.set_boundary_condition("top", bc.FixedValue(0))
defaultExternalSourceFlux.set_boundary_condition("bottom", bc.FixedValue(0))

timeFolder0.append(defaultExternalSourceFlux)


#==============================================================================*
# Solvers

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    region=nMesh.region,
    solver="diffusionNeutronics",
    eigenvalueNeutronics=False,
    externalSourceNeutronics=True,
    isSetFvSolutionToDefault=True,
    keff=1, # keff must be set to 1 for external source problems
)

fluxSolution: ffn.numerics.fvSolutionSolver = neutronicsSolver.fvSolution.solvers['"flux.*"']
fluxSolution.relTol = 1e-10

print(f"List of require fields: {neutronicsSolver.get_required_fields()}")

# External source modulation time profile (table)
source_table = ffn.timeProfile.TimeProfile( 'table',
    startTime=0,
    table=[
        (0, 0.01),
        (1, 0.01),
        (2, 1),
        (4, 1),
        (5, 0)
    ]
)

neutronicsSolver.externalSource.externalSourceModulationTimeProfile = source_table

neutronicsSolver.mesh = nMesh

neutronicsSolver.fvSchemes.ddtSchemes['default'] = 'Euler'

nuclearData = neutronicsSolver.nuclearData

refState = ffn.nuclearData.NuclearDataState(
    name="reference",
    zones=[
        ffn.nuclearData.NuclearDataZone(
            "zone0",
            fuelFraction=1,
            removalXS=[sigma_a],
            nuFissionXS=[nu_sigma_f],
            powerXS=[1],
            scatteringMatrixP0=[[0.1]],
            discFactor=[1],
            chiPrompt=[1],
            chiDelayed=[1],
            inverseVelocity=[7.1e-8],
            diffusionCoefficient=[D],
            integralFlux=[1],
            decayConstant=[
                0.0125371, 0.0300828, 0.109879, 0.125484, 1.3036, 9.51817
            ],
            delayedFraction=[
                0, 0, 0, 0, 0, 0
            ],
        )
    ]
)

nuclearData.add_state(refState)


#==============================================================================*
# Settings

model = ffn.case.Case()

model.solvers.append(neutronicsSolver)
model.timeFolders = [timeFolder0]

settings = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 6
settings.deltaT = 1e-9
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 0.5
settings.runTimeModifiable = True
settings.maxPowerVariation = 0.05
settings.adjustTimeStep = True

print(model)


#==============================================================================*

# Define analytical solutions for different k_inf regimes
def get_analytical_solution(k_inf_value, L_sq, H_val, S0_val, D_val):

    #Returns the appropriate analytical flux function based on k_inf value.

    if k_inf_value==1:
        # Critical case: k_inf = 1
        def flux_critical(x):
            """Φ(x) = S₀/(2D) * (H²/4 - x²)"""
            return (S0_val / (2 * D_val)) * (H_val**2 / 4 - x**2)

        return flux_critical

    elif k_inf_value > 1.0:
        # Supercritical case: k_inf > 1
        alpha_sq = (k_inf_value - 1) / L_sq
        alpha_val = np.sqrt(alpha_sq)

        def flux_supercritical(x):
            """Φ(x) = S₀/(Dα²) * [cos(αx)/cos(αH/2) - 1]"""
            return (S0_val / (D_val * alpha_sq)) * (np.cos(alpha_val * x) / np.cos(alpha_val * H_val / 2) - 1)

        return flux_supercritical

    else:
        # Subcritical case: k_inf < 1
        alpha_sq = (1 - k_inf_value) / L_sq
        alpha_val = np.sqrt(alpha_sq)

        def flux_subcritical(x):
            """Φ(x) = S₀/(Dα²) * [1 - cosh(αx)/cosh(αH/2)]"""
            return (S0_val / (D_val * alpha_sq)) * (1 - np.cosh(alpha_val * x) / np.cosh(alpha_val * H_val / 2))

        return flux_subcritical


#==============================================================================*
# Run and post-processing

if __name__ == '__main__':
    #==========================================================================*
    # Clean, export and run

    ffn.allclean()
    model.export_to_openfoam()

    ffn.run(model, is_preprocessing=True)


    #==========================================================================*
    # Post-processing

    model.plot_mesh(region=nMesh)

    model.plot_slice(
        region=nMesh,
        time=source_table.table[-2][0],
        fieldName="flux0",
        show_edges=False,
        cmap="Blues_r",
        unit="n/m2/s"
    )

    model.plot_mesh(
        region=nMesh,
        time=source_table.table[-2][0],
        fieldName="flux0",
        cmap="Blues_r",
        show_edges=False
    )

    points, flux0 = model.sample_over_line(
        region=nMesh,
        time=source_table.table[-2][0],
        fieldName="flux0",
        point1=(0, 0, H/nz * 0.5),
        point2=(0, 0, H/nz * (nz-0.5)),
    )

    z = points[:, 2]

    # Transform coordinates: simulation uses [0, H], analytical solution uses [-H/2, H/2]
    z_centered = z - H/2

    # Get appropriate analytical solution
    func_flux_analytical = get_analytical_solution(k_inf, L_squared, H, S0, D)
    flux_analytical = func_flux_analytical(z_centered)

    # Normalize both solutions for comparison
    flux_sim_norm = flux0 / np.max(flux0)
    flux_ana_norm = flux_analytical / np.max(flux_analytical)

    # Calculate relative error
    rel_error = (flux0 - flux_analytical) / flux_analytical * 100

    # Plotting
    fig, (ax_flux, ax_error) = plt.subplots(nrows=2, figsize=(5, 6), dpi=200)

    # Plot 1: Flux comparison
    ax_flux.plot(z, flux0, 'o-', label="GeN-Foam", markersize=3, alpha=0.7)
    ax_flux.plot(z, flux_analytical, 'k--', label=f"Analytical (k_inf={k_inf:.4f})", linewidth=2)
    ax_flux.set_ylabel("Neutron flux [n/(m²s)]")
    ax_flux.set_title(f"Flux Distribution (k_inf = {k_inf:.4f})")
    ax_flux.legend()
    ax_flux.grid(True, alpha=0.3)
    ax_flux.set_xlim(0, H)

    # Plot 2: Relative error
    ax_error.plot(z, rel_error, 'r-', linewidth=2)
    ax_error.set_xlabel("Axial position [m]")
    ax_error.set_ylabel("Relative error [%]")
    ax_error.set_title("Relative Error")
    ax_error.grid(True, alpha=0.3)
    ax_error.set_xlim(0, H)

    fig.tight_layout()
    fig.savefig("fig_results_fluxDistribution.png")
    plt.close()

    # Normalized flux comparison with error
    rel_error_norm = (flux_sim_norm - flux_ana_norm) / flux_ana_norm * 100

    fig_norm, (ax_norm, ax_error_norm) = plt.subplots(nrows=2, figsize=(5, 6), dpi=200)

    # Plot 1: Normalized flux comparison
    ax_norm.plot(z, flux_sim_norm, 'o-', label="GeN-Foam", markersize=3, alpha=0.7)
    ax_norm.plot(z, flux_ana_norm, 'k--', label=f"Analytical (k_inf={k_inf:.4f})", linewidth=2)
    ax_norm.set_ylabel("Normalized neutron flux [a.u.]")
    ax_norm.set_title(f"Normalized Flux Distribution (k_inf = {k_inf:.4f})")
    ax_norm.set_xlim(0, H)
    ax_norm.grid(True, alpha=0.3)
    ax_norm.legend()

    # Plot 2: Relative error
    ax_error_norm.plot(z, rel_error_norm, 'r-', linewidth=2)
    ax_error_norm.set_xlabel("Axial position [m]")
    ax_error_norm.set_ylabel("Relative error [%]")
    ax_error_norm.set_title("Relative Error - Normalized Flux")
    ax_error_norm.set_xlim(0, H)
    ax_error_norm.grid(True, alpha=0.3)

    fig_norm.tight_layout()
    fig_norm.savefig("fig_results_normFluxDistribution.png")
    plt.close()

    # Time evolution

    model.plot_animation(
        region=nMesh,
        fieldName="flux0",
        fps=1,
        cmap="Blues_r",
        limits=[0, 2000]
    )

    # Power Plot
    power_profile = model.get_parameters_from_log(parameters=[("power", "power", 1)])
    fig, axPower = plt.subplots(figsize=(6, 4), dpi=200)

    # Plot power on primary axis
    axPower.plot(power_profile["time"], power_profile["power"], 'o-', label="Power", color='C0', markersize=4)
    axPower.set_ylabel("Total Power [W]", color='C0')
    axPower.tick_params(axis='y', labelcolor='C0')
    axPower.set_xlabel("Time [s]")
    axPower.set_xlim((0, power_profile["time"][-1] + settings.writeInterval))
    axPower.set_ylim(0)

    # Create secondary axis for source
    axSource = axPower.twinx()

    # Extract source modulation times and values from source_table
    source_times = [t for t, v in source_table.table]
    source_values = [v * S0 for t, v in source_table.table]

    # Plot source with linear interpolation between points on secondary axis
    axSource.plot(source_times, source_values, label="Source", color='C1', marker='o')
    axSource.set_ylabel("Source [1/(m³s)]", color='C1')
    axSource.tick_params(axis='y', labelcolor='C1')
    axSource.set_ylim(0)

    # Combined legend
    lines1, labels1 = axPower.get_legend_handles_labels()
    lines2, labels2 = axSource.get_legend_handles_labels()
    axPower.legend(lines1 + lines2, labels1 + labels2, loc='best')

    fig.tight_layout()
    fig.savefig("fig_results_powerProfile.png")
    plt.close()


#==============================================================================*
