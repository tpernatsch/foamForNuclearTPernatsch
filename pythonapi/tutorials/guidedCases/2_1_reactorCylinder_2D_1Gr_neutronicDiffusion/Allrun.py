"""
2D Cylindrical Reactor (R-Z geometry) - Neutronics Diffusion
"""
#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import j0, jn_zeros

import fluidfoam  # type: ignore
import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# Parameters

rOut = 0.6  # Reactor radius [m]
fuelHeight = 1.3  # Reactor height [m]
nr = 11  # Number of radial cells
nz = 51  # Number of axial cells
wedgeAngle = 1  # Wedge angle [degrees]
threshold = 0.05  # Threshold for error calculation (fraction of max flux)

# Nuclear data
lambda_0 = jn_zeros(0, 1)[0]  # First zero of J₀: λ₀ ≈ 2.4048
nuSigmaf = 7.80108  # m⁻¹
Sigmar = 5.0082  # m⁻¹
D = 0.1275  # m
keff_anal = nuSigmaf / (Sigmar + D * ((lambda_0/rOut)**2 + (np.pi/fuelHeight)**2))


#==============================================================================*
# Analytical solution

def analytical_flux(r, z):
    """Analytical flux: φ(r,z) = J₀(λ₀·r/R) · sin(π·z/H)"""
    return j0(lambda_0 * r / rOut) * np.sin(np.pi * z / fuelHeight)


#==============================================================================*
# Mesh

def createMesh(nr: int, nz: int):
    nMesh = mesh.BlockMesh(region="neutroRegion")

    zone0 = nMesh.create_wedge(
        "zone0", 0, rOut, 0, fuelHeight, wedgeAngle=wedgeAngle, nr=nr, nz=nz
    )

    wedgeFaces = [
        ("front", zone0.frontFace()),
        ("back", zone0.backFace()),
    ]

    for facename, subface in wedgeFaces:
        face = ffn.mesh.Face(facename, boundaryType="wedge")
        face.add_sub_face(subface)
        nMesh.add_boundary(face)

    outer = ffn.mesh.Face("outer", boundaryType="wall")
    outer.add_sub_face(zone0.rightFace())
    nMesh.add_boundary(outer)

    inner = ffn.mesh.Face("inner", boundaryType="empty")
    inner.add_sub_face(zone0.leftFace())
    nMesh.add_boundary(inner)

    top = ffn.mesh.Face("top", boundaryType="wall")
    top.add_sub_face(zone0.topFace())
    nMesh.add_boundary(top)

    bottom = ffn.mesh.Face("bottom", boundaryType="wall")
    bottom.add_sub_face(zone0.bottomFace())
    nMesh.add_boundary(bottom)

    return(nMesh, wedgeFaces)

nMesh, wedgeFaces = createMesh(nr=nr, nz=nz)


#==============================================================================*
# Fields

timeFolder0 = ffn.timeFolder.TimeFolder(0)

defaultFlux = ffn.fields.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.fields.Dimension(default='flux')
defaultFlux.internalField = 1e4
for facename, _ in wedgeFaces:
    defaultFlux.set_boundary_condition(facename, bc.Wedge())
defaultFlux.set_boundary_condition("top", bc.FixedValue(0))
defaultFlux.set_boundary_condition("bottom", bc.FixedValue(0))
defaultFlux.set_boundary_condition("outer", bc.FixedValue(0))
defaultFlux.set_boundary_condition("inner", bc.Empty())

timeFolder0.append(defaultFlux)

#==============================================================================*
# Solvers

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    region=nMesh.region,
    solver="diffusionNeutronics",
    isSetFvSolutionToDefault=True,
)

print(f"List of require fields: {neutronicsSolver.get_required_fields()}")

neutronicsSolver.mesh = nMesh

neutronicsSolver.fvSchemes.ddtSchemes['default'] = 'steadyState'

nuclearData = neutronicsSolver.nuclearData

refState = ffn.nuclearData.NuclearDataState(
    name="reference",
    zones=[
        ffn.nuclearData.NuclearDataZone(
            "zone0",
            fuelFraction=0.4,
            removalXS=[Sigmar],
            nuFissionXS=[nuSigmaf],
            powerXS=[1],
            scatteringMatrixP0=[[0.1]],
            discFactor=[1],
            chiPrompt=[1],
            chiDelayed=[1],
            inverseVelocity=[7.1e-08],
            diffusionCoefficient=[D],
            integralFlux=[1],
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

nuclearData.add_state(refState)

#==============================================================================*
# Settings

model = ffn.case.Case()

model.solvers.append(neutronicsSolver)
model.timeFolders = [timeFolder0]

settings = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 1
settings.deltaT = 1e-9
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 1
settings.runTimeModifiable = True
settings.adjustTimeStep = True


#==============================================================================*
# Run

if __name__ == "__main__":
    print(model)

    # Export to OpenFOAM
    ffn.allclean()
    model.export_to_openfoam()

    ffn.run(model, is_preprocessing=True)


    # Post-processing
    #----------------

    print(f"Analytical k_eff = {keff_anal:.6f}")
    print(f"keff = {model.keff()}")

    model.plot_mesh(region=nMesh)

    model.plot_mesh(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux0",
        cmap="Blues_r",
        show_edges=False,
        normal='y'
    )

    model.plot_residuals(
        parameters=['fluxStar0'],
        title="Neutronics"
    )

    # Read mesh and flux data
    X, Y, Z = fluidfoam.readmesh(".", region=nMesh.region, structured=False)
    flux0 = fluidfoam.readscalar(".", time_name="1", name="flux0", region=nMesh.region, structured=False)
    R = np.sqrt(X**2 + Y**2)  # Radial coordinate

    # Normalize solutions and compute error
    flux0_norm = flux0 / np.max(flux0)
    flux_analytical_norm = analytical_flux(R, Z) / np.max(analytical_flux(R, Z))

    # Apply threshold mask and calculate relative error
    mask = flux_analytical_norm > threshold
    relative_error = np.full_like(flux0_norm, np.nan)
    relative_error[mask] = ((flux0_norm[mask] / flux_analytical_norm[mask]) - 1) * 100

    # Reshape data to 2D grid
    sort_idx = np.lexsort((Z, R))
    flux0_grid = flux0_norm[sort_idx].reshape(nr, nz).T
    flux_analytical_grid = flux_analytical_norm[sort_idx].reshape(nr, nz).T
    relative_error_grid = relative_error[sort_idx].reshape(nr, nz).T
    mask_grid = mask[sort_idx].reshape(nr, nz).T
    R_grid = R[sort_idx].reshape(nr, nz).T
    Z_grid = Z[sort_idx].reshape(nr, nz).T
    extent = [0, rOut, 0, fuelHeight]


    # Plot 1: Numerical flux (GeN-Foam)
    #----------------------------------

    fig1, ax1 = plt.subplots(figsize=(5, 4), dpi=200)
    im1 = ax1.imshow(
        flux0_grid,
        cmap='Blues_r',
        origin='lower',
        extent=extent,
        aspect='auto',
        vmin=0,
        vmax=1
    )
    ax1.set_xlabel('Radial position [m]')
    ax1.set_ylabel('Axial position [m]')
    ax1.set_title('GeN-Foam normalized flux')
    plt.colorbar(im1, ax=ax1, label='Normalized flux [a.u.]')
    fig1.tight_layout()
    fig1.savefig("fig_results_flux2D_GenFoam.png")
    plt.close()


    # Plot 2: Analytical flux
    #------------------------

    fig2, ax2 = plt.subplots(figsize=(5, 4), dpi=200)
    im2 = ax2.imshow(
        flux_analytical_grid,
        cmap='Blues_r',
        origin='lower',
        extent=extent,
        aspect='auto',
        vmin=0,
        vmax=1
    )
    ax2.set_xlabel('Radial position [m]')
    ax2.set_ylabel('Axial position [m]')
    ax2.set_title('Analytical normalized flux')
    plt.colorbar(im2, ax=ax2, label='Normalized flux [a.u.]')
    fig2.tight_layout()
    fig2.savefig("fig_results_flux2D_analytical.png")
    plt.close()


    # Plot 3: Relative error
    #-----------------------

    fig3, ax3 = plt.subplots(figsize=(5, 4), dpi=200)
    im3 = ax3.imshow(
        relative_error_grid,
        cmap='viridis',
        origin='lower',
        extent=extent,
        aspect='auto'
    )
    ax3.set_xlabel('Radial position [m]')
    ax3.set_ylabel('Axial position [m]')
    ax3.set_title(f'Relative error (flux > {threshold*100:.0f}% max)')
    plt.colorbar(im3, ax=ax3, label='Relative error [%]')
    fig3.tight_layout()
    fig3.savefig("fig_results_flux2D_error.png")
    plt.close()


    # Plot 4: Threshold mask visualization
    #-------------------------------------

    fig4, ax4 = plt.subplots(figsize=(5, 4), dpi=200)
    im4 = ax4.imshow(
        mask_grid.astype(float),
        cmap='gray',
        origin='lower',
        extent=extent,
        aspect='auto',
        vmin=0,
        vmax=1
    )
    ax4.set_xlabel('Radial position [m]')
    ax4.set_ylabel('Axial position [m]')
    ax4.set_title(f'Region with flux > {threshold*100:.0f}% max')
    plt.colorbar(im4, ax=ax4, label='Included in error [-]', ticks=[0, 1])
    fig4.tight_layout()
    fig4.savefig("fig_results_flux2D_mask.png")
    plt.close()

    # Extract 1D profiles at middle coordinates
    mid_z_idx, mid_r_idx = nz // 2, nr // 2

    # Radial profile at mid-height (extract from actual mesh coordinates)
    r_profile = R_grid[mid_z_idx, :]
    flux_r = flux0_grid[mid_z_idx, :]
    flux_r_anal = flux_analytical_grid[mid_z_idx, :]
    mask_r = flux_r_anal > threshold
    rel_error_r = np.full_like(flux_r, np.nan)
    rel_error_r[mask_r] = ((flux_r[mask_r] / flux_r_anal[mask_r]) - 1) * 100

    # Axial profile at mid-radius (extract from actual mesh coordinates)
    z_profile = Z_grid[:, mid_r_idx]
    flux_z = flux0_grid[:, mid_r_idx]
    flux_z_anal = flux_analytical_grid[:, mid_r_idx]
    mask_z = flux_z_anal > threshold
    rel_error_z = np.full_like(flux_z, np.nan)
    rel_error_z[mask_z] = (flux_z[mask_z] - flux_z_anal[mask_z]) / flux_z_anal[mask_z] * 100


    # Plot 5: Radial profile at mid-height
    #-------------------------------------

    fig5, (ax5a, ax5b) = plt.subplots(nrows=2, figsize=(5, 6), dpi=200, sharex=True)
    ax5a.plot(r_profile, flux_r, 'o-', label='GeN-Foam', markersize=4)
    ax5a.plot(r_profile, flux_r_anal, 'k--', label='Analytical', linewidth=2)
    ax5a.set_xlim(0, rOut)
    ax5a.set_ylim(bottom=0)
    ax5a.set_ylabel('Normalized flux [a.u.]')
    ax5a.set_title(f'Radial profile at Z = {Z_grid[mid_z_idx, 0]:.3f} m')
    ax5a.legend()
    ax5a.grid(True, alpha=0.3)

    ax5b.plot(r_profile, rel_error_r, 'r-', linewidth=2)
    ax5b.set_xlim(0, rOut)
    ax5b.set_xlabel('Radial position [m]')
    ax5b.set_ylabel('Relative error [%]')
    ax5b.grid(True, alpha=0.3)
    fig5.tight_layout()
    fig5.savefig("fig_results_radial_profile.png")
    plt.close()


    # Plot 6: Axial profile at mid-radius
    #------------------------------------

    fig6, (ax6a, ax6b) = plt.subplots(nrows=2, figsize=(5, 6), dpi=200, sharex=True)
    ax6a.plot(z_profile, flux_z, 'o-', label='GeN-Foam', markersize=4)
    ax6a.plot(z_profile, flux_z_anal, 'k--', label='Analytical', linewidth=2)
    ax6a.set_xlim(0, fuelHeight)
    ax6a.set_ylabel('Normalized flux [a.u.]')
    ax6a.set_title(f'Axial profile at R = {R_grid[0, mid_r_idx]:.3f} m')
    ax6a.legend()
    ax6a.grid(True, alpha=0.3)

    ax6b.plot(z_profile, rel_error_z, 'r-', linewidth=2)
    ax6b.set_xlim(0, fuelHeight)
    ax6b.set_xlabel('Axial position [m]')
    ax6b.set_ylabel('Relative error [%]')
    ax6b.grid(True, alpha=0.3)
    fig6.tight_layout()
    fig6.savefig("fig_results_axial_profile.png")
    plt.close()


#==============================================================================*
