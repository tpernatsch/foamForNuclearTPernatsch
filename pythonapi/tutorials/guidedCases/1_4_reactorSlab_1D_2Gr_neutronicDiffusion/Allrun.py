"""

"""
#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import numpy as np
import tabulate

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh

#==============================================================================*

ffn.allclean()

#==============================================================================*
# Nuclear data parameters - PWR two-group cross sections

# Group 0 = fast, Group 1 = thermal
nu_sigma_f = np.array([0.85, 14.20472])  # nu*Sigma_fission [m⁻¹]
removal_xs = np.array([3.62, 12.1])  # Removal cross-section [m⁻¹]
# kappa: energy released per fission (200 MeV → Joules)
kappa_J = 200.0e6 * 1.60218e-19 # J/fission
# Assume average nu (neutrons per fission)
nu = 2.4
sigma_fission = nu_sigma_f / nu
power_xs = kappa_J * sigma_fission

# Scattering matrix: scatter[from][to] - upper triangular (no upscattering)
# Self-scattering terms are small in this example
scattering_matrix = [[0.5, 2.41], [0, 0.5]]  # Scattering matrix P0 [m⁻¹]

diffusion_coefficient = [0.01267, 0.00354]  # Diffusion coefficient [m]
chiPrompt = [1, 0]  # Fission spectrum (all fast)
chiDelayed = [1, 0]  # Fission spectrum (all fast)


#==============================================================================*
# Geometry parameters

fuelLength = 1.5
lx = 0.1
ly = 0.1

target_power = 1  # W (from neutronicsSolver power setting)
cross_section_area = lx * ly  # m^2 (from mesh dimensions)

#==============================================================================*
# Mesh

def createMesh(nz: int):
    nMesh = mesh.BlockMesh(region='neutroRegion')

    zone0 = nMesh.create_cube('zone0', -lx/2, -ly/2, 0, lx/2, ly/2, fuelLength, 1, 1, nz)

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
defaultFlux.internalField = 1e16
defaultFlux.set_boundary_condition("walls", bc.ZeroGradient())
defaultFlux.set_boundary_condition("top", bc.FixedValue(0))
defaultFlux.set_boundary_condition("bottom", bc.FixedValue(0))

timeFolder0.append(defaultFlux)


#==============================================================================*
# Solvers

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    region="neutroRegion",
    solver="diffusionNeutronics",
    power=target_power,
    mesh=nMesh,
    isSetFvSolutionToDefault=False
)

print(f"List of require fields: {neutronicsSolver.get_required_fields()}")

neutronicsSolver.fvSchemes.ddtSchemes['default'] = 'steadyState'

neutronicsSolution = ffn.numerics.fvSolution()
neutronicsSolution.append('".*"', ffn.numerics.fvSolutionSolver(
    solver='PCG',
    preconditioner='DIC',
    tolerance=1e-10,
    relTol=1e-5
))

neutronicsSolver.fvSolution = neutronicsSolution

refState = ffn.nuclearData.NuclearDataState(
    name="reference",
    zones=[
        ffn.nuclearData.NuclearDataZone(
            "zone0",
            fuelFraction=0.4,
            removalXS=removal_xs,
            nuFissionXS=nu_sigma_f,
            powerXS=power_xs,
            scatteringMatrixP0=scattering_matrix,
            discFactor=[1, 1],
            chiPrompt=chiPrompt,
            chiDelayed=chiDelayed,
            inverseVelocity=[7.1e-08, 7.1e-08],
            diffusionCoefficient=diffusion_coefficient,
            integralFlux=[1, 1],
            decayConstant=[
                0.0125371, 0.0300828, 0.109879,
                0.325484, 1.3036, 9.51817
            ],
            delayedFraction=[
                7.2315e-05, 0.000609661, 0.000471181,
                0.00118907, 0.000445487, 9.58515e-05
            ],
        )
    ]
)

neutronicsSolver.nuclearData.add_state(refState)


#==============================================================================*
# Settings

model = ffn.case.Case(timeFolders=[timeFolder0])

model.solvers.append(neutronicsSolver)

settings = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 1
settings.deltaT = 1e-9
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 1
settings.runTimeModifiable = True
settings.adjustTimeStep = True

print(model)


#==============================================================================*
# Run and Post-processing (only when run directly)

if __name__ == '__main__':
    # Export to OpenFOAM
    model.export_to_openfoam()

    # Run
    ffn.run(model, is_preprocessing=True)


    # Post-processing
    model.plot_mesh(region=nMesh)

    for group in range(2):
        model.plot_slice(
            region=nMesh,
            time=settings.endTime,
            fieldName=f"flux{group}",
            cmap='Blues_r',
            unit="n/m2/s"
        )
        model.plot_mesh(
            region=nMesh,
            time=settings.endTime,
            fieldName=f"flux{group}",
            cmap="Blues_r",
            unit="n/m2/s"
        )

    points, flux0 = model.sample_over_line(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux0",
        point1=(0, 0, fuelLength/nz * 0.5),
        point2=(0, 0, fuelLength/nz * (nz-0.5)),
    )

    points, flux1 = model.sample_over_line(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux1",
        point1=(0, 0, fuelLength/nz * 0.5),
        point2=(0, 0, fuelLength/nz * (nz-0.5)),
    )

    z = points[:, 2]

    # Analytical solution calculation
    B_squared = (np.pi / fuelLength)**2

    # Calculate analytical k-effective using two-group formula
    # k = νΣ_f^1/(Σ_a^1 + Σ_s^(1→2) + D^1 B^2) + [Σ_s^(1→2)/(Σ_a^1 + Σ_s^(1→2) + D^1 B^2)] * [νΣ_f^2/(Σ_a^2 + D^2 B^2)]
    denominator_fast = removal_xs[0] + diffusion_coefficient[0] * B_squared
    denominator_thermal = removal_xs[1] + diffusion_coefficient[1] * B_squared

    k_analytical = (
        nu_sigma_f[0] / denominator_fast
        + (scattering_matrix[0][1] / denominator_fast) * (nu_sigma_f[1] / denominator_thermal)
    )

    # Get simulated k-effective
    k_simulated = model.keff()

    # Calculate relative error in k-effective
    k_rel_error = (k_simulated - k_analytical) / k_analytical * 100

    print(f"K-effective Comparison:")
    print(tabulate.tabulate(
        [
            ["Analytical", f"{k_analytical:.6f}"],
            ["Simulated", f"{k_simulated:.6f}"],
            ["Rel. Error [%]", f"{k_rel_error:.4f}"],
        ],
        headers=["", "keff"],
        tablefmt='pipe'
    ))

    # Flux ratio φ₀/φ₁ (fast/thermal) from thermal group balance
    flux_ratio_fast_to_thermal = (
        diffusion_coefficient[1] * B_squared + removal_xs[1]
    ) / scattering_matrix[0][1]

    # Amplitude calculation
    amplitude_1 = target_power * np.pi / (
        2 * fuelLength * cross_section_area * kappa_J * (
            sigma_fission[0] * flux_ratio_fast_to_thermal + sigma_fission[1]
        )
    )
    amplitude_0 = amplitude_1 * flux_ratio_fast_to_thermal

    # Analytical flux distributions
    flux_analytical_shape = np.sin(z * np.pi / fuelLength)
    flux_analytical_0 = amplitude_0 * flux_analytical_shape
    flux_analytical_1 = amplitude_1 * flux_analytical_shape

    # Relative errors for flux
    rel_error_0 = (flux0 - flux_analytical_0) / flux_analytical_0 * 100
    rel_error_1 = (flux1 - flux_analytical_1) / flux_analytical_1 * 100


    # Plot 1: Flux and relative error comparison
    #===========================================

    fig1, axes = plt.subplots(2, 1, figsize=(5, 6), dpi=200)

    ax_flux = axes[0]
    ax_flux.plot(z, flux0, 'o-', label="GeN-Foam g=0 (fast)", markersize=3, alpha=0.7)
    ax_flux.plot(z, flux_analytical_0, 'k--', label="Analytical g=0", linewidth=2)
    ax_flux.plot(z, flux1, 's-', label="GeN-Foam g=1 (thermal)", markersize=3, alpha=0.7)
    ax_flux.plot(z, flux_analytical_1, 'r--', label="Analytical g=1", linewidth=2)
    ax_flux.set_ylabel("Neutron flux [n/(m$^2$ s)]")
    ax_flux.set_title("Flux Distribution")
    ax_flux.legend()
    ax_flux.set_xlim(0, fuelLength)

    ax_error = axes[1]
    ax_error.plot(z, rel_error_0, label="g=0 (fast)")
    ax_error.plot(z, rel_error_1, label="g=1 (thermal)", linestyle='--')
    ax_error.set_xlabel("Axial position [m]")
    ax_error.set_ylabel("Relative error [%]")
    ax_error.set_title("Relative Error")
    ax_error.legend()
    ax_error.set_xlim(0, fuelLength)

    fig1.tight_layout()
    fig1.savefig("fig_results_fluxDistribution.png")
    plt.close()


    # Plot 2: Normalized flux comparison
    #===================================

    fig2, axes = plt.subplots(2, 1, figsize=(5, 6), dpi=200)

    flux0_norm = flux0 / np.max(flux0)
    flux1_norm = flux1 / np.max(flux1)
    flux_analytical_norm = flux_analytical_shape / np.max(flux_analytical_shape)

    ax_norm = axes[0]
    ax_norm.plot(z, flux0_norm, 'o-', label="GeN-Foam g=0 (fast)", markersize=3, alpha=0.7)
    ax_norm.plot(z, flux1_norm, 's-', label="GeN-Foam g=1 (thermal)", markersize=3, alpha=0.7)
    ax_norm.plot(z, flux_analytical_norm, 'k--', label="Analytical", linewidth=2)
    ax_norm.set_ylabel("Normalized neutron flux [a.u.]")
    ax_norm.set_title("Normalized Flux Distribution")
    ax_norm.legend()
    ax_norm.set_xlim(0, fuelLength)

    rel_error_norm_0 = (flux0_norm - flux_analytical_norm) / flux_analytical_norm * 100
    rel_error_norm_1 = (flux1_norm - flux_analytical_norm) / flux_analytical_norm * 100

    ax_error_norm = axes[1]
    ax_error_norm.plot(z, rel_error_norm_0, label="g=0 (fast)")
    ax_error_norm.plot(z, rel_error_norm_1, label="g=1 (thermal)", linestyle='--')
    ax_error_norm.set_xlabel("Axial position [m]")
    ax_error_norm.set_ylabel("Relative error [%]")
    ax_error_norm.set_title("Normalized Flux Relative Error")
    ax_error_norm.legend()
    ax_error_norm.set_xlim(0, fuelLength)

    fig2.tight_layout()
    fig2.savefig("fig_results_normFluxDistribution.png")
    plt.close()


#==============================================================================*
