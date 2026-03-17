#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import numpy as np
import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# Mesh

def createMesh(nz: int):
    nMesh = mesh.BlockMesh(region='neutroRegion')

    zone0 = nMesh.create_cube('zone0', -widthX/2, -widthY/2, 0, widthX/2, widthY/2, fuelLength, 1, 1, nz)

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

fuelLength = 3/2
widthX = widthY = 0.1

nMesh = createMesh(nz=nz)


#==============================================================================*

# Parameters (k_inf = nuSigmaf/removalXS = 1)
D = 0.1275  # Diffusion coefficient [1/m]
removalXS = 12.725  # Removal cross-section [1/m]
nuSigmaf = 10  # nu*Sigma_f [1/m]
externalSourceCurrent = 3e14  # External source [1/(m^2*s)]
fluxGradient = externalSourceCurrent / D


#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

defaultFlux = ffn.fields.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.fields.Dimension(default='flux')
defaultFlux.internalField = 1e21

defaultFlux.set_boundary_condition("walls",bc.ZeroGradient())
defaultFlux.set_boundary_condition(
    "bottom",
    bc.CustomPatch(parameters={
        "type": "fixedGradient",
        "gradient": f"uniform {fluxGradient}"
    })
)
defaultFlux.set_boundary_condition("top", bc.FixedValue(0))

timeFolder0.append(defaultFlux)


#==============================================================================*
# Solvers

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    region=nMesh.region,
    solver="diffusionNeutronics",
    eigenvalueNeutronics=False
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
            fuelFraction=1,
            removalXS=[removalXS],
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

model = ffn.Case()
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

print(model)

# Export to OpenFOAM
model.export_to_openfoam()

# Analytical solution function
def flux_analytical(z, nuSigmaf, removalXS, D, a, J0):
    B2 = (nuSigmaf - removalXS) / D
    if np.isclose(B2, 0.0, atol=1e-12):  # Critical case
        return (J0 / D) * (a - z)
    B = np.sqrt(B2) if B2 > 0 else 1j * np.sqrt(-B2)
    flux = (J0 / (D * B * np.cos(B * a))) * np.sin(B * (a - z))
    return np.real_if_close(flux)

# Calculate k_inf
k_inf = nuSigmaf / removalXS


#==============================================================================*
# Main

if __name__ == '__main__':
    ffn.run(model, is_preprocessing=True)

    #==========================================================================*
    # Post-processing
    model.plot_mesh(region=nMesh)
    model.plot_boundary(region=nMesh, boundaryName="walls")

    model.plot_slice(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux0",
        show_edges=False,
        cmap='Blues_r',
        unit="n/m2/s"
    )

    model.plot_mesh(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux0",
        cmap="Blues_r",
        show_edges=False
    )

    points, flux0 = model.sample_over_line(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux0",
        point1=(0, 0, fuelLength/nz * 0.5),
        point2=(0, 0, fuelLength/nz * (nz-0.5))
    )
    z = points[:,2]

    print(f"k_inf = {k_inf:.6f}")

    fluxTh = lambda z_: flux_analytical(z_, nuSigmaf, removalXS, D, fuelLength, externalSourceCurrent)

    relError = [(flux0_ / fluxTh(z_) - 1) * 100 for z_, flux0_ in zip(z, flux0)]

    # Plot absolute flux comparison
    fig, (axFlux, axError) = plt.subplots(nrows=2, figsize=(5, 6), dpi=200)
    axFlux.plot(z, flux0, label="GeN-Foam")
    axFlux.plot(z, fluxTh(z), label="Analytic", ls="--")
    axFlux.set_ylabel("Neutron flux [a.u]")
    axFlux.set_title("Flux Distribution ($k_{inf}$ = "+f"{k_inf:.4f})")
    axFlux.legend()
    axError.plot(z, relError)
    axError.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    axError.set_ylabel("Relative error [%]")
    for ax in [axFlux, axError]:
        ax.set_xlabel("Axial position [m]")
        ax.set_xlim((0, fuelLength))
    fig.tight_layout()
    fig.savefig("fig_results_fluxDistribution.png")
    plt.close()

    # Plot normalized flux comparison
    flux0 = np.array(flux0)
    fluxTh_vals = fluxTh(z)
    normFlux_sim = flux0 / np.max(flux0)
    normFlux_ana = fluxTh_vals / np.max(fluxTh_vals)
    relError_norm = (normFlux_sim / normFlux_ana - 1) * 100

    fig, (axFlux, axError) = plt.subplots(nrows=2, figsize=(5, 6), dpi=200)
    axFlux.plot(z, normFlux_sim, label="GeN-Foam")
    axFlux.plot(z, normFlux_ana, label="Analytic", ls="--")
    axFlux.set_ylabel("Neutron flux [a.u]")
    axFlux.set_title("Normalized Flux Distribution ($k_{inf}$ = "+f"{k_inf:.4f})")
    axFlux.legend()
    axError.plot(z, relError_norm)
    axError.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    axError.set_ylabel("Relative error [%]")
    for ax in [axFlux, axError]:
        ax.set_xlabel("Axial position [m]")
        ax.set_xlim((0, fuelLength))
    fig.tight_layout()
    fig.savefig("fig_results_normFluxDistribution.png")
    plt.close()


#==============================================================================*
