"""

"""
#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import numpy as np

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh
import fluidfoam  # type: ignore


#==============================================================================*
# Mesh

def createMesh(dr: float, nr: int):
    opening = 2 # degrees

    nMesh = mesh.BlockMesh(region='neutroRegion')

    zone0 = nMesh.create_sphere_1D(
        name='sph1',
        innerRadius=0, outerRadius=dr,
        opening=opening,
        nr=nr
    )

    wedgeFaces = [
        ("front",  zone0.frontFace()),
        ("back",   zone0.backFace()),
        ("top",    zone0.topFace()),
        ("bottom", zone0.bottomFace()),
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

    return(nMesh, wedgeFaces)


dr = 1

nMesh, wedgeFaces = createMesh(dr=dr, nr=101)


#==============================================================================*
# Fields

timeFolder0 = ffn.TimeFolder(0)

defaultFlux = ffn.fields.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.fields.Dimension(default='flux')
defaultFlux.internalField = 1
for facename, _ in wedgeFaces:
    defaultFlux.set_boundary_condition(facename, bc.Wedge())
defaultFlux.set_boundary_condition("outer", bc.FixedValue(0))
defaultFlux.set_boundary_condition("inner", bc.Empty())

timeFolder0.append(defaultFlux)


#==============================================================================*
# Solvers

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    region=nMesh.region,
    solver="diffusionNeutronics",
    mesh=nMesh,
    isSetFvSolutionToDefault=False,
)

neutronicsSolution = ffn.numerics.fvSolution()
neutronicsSolution.append('".*"', ffn.numerics.fvSolutionSolver(
    solver='PCG',
    preconditioner='DIC',
    tolerance=1e-9,
    relTol=1e-3
))

neutronicsSolver.fvSolution = neutronicsSolution

print(f"List of require fields: {neutronicsSolver.get_required_fields()}")

refState = ffn.nuclearData.NuclearDataState(
    name="reference",
    zones=[
        ffn.nuclearData.NuclearDataZone(
            "sph1",
            fuelFraction=0.4,
            removalXS=[4.612608],
            nuFissionXS=[4.8101],
            powerXS=[1],
            scatteringMatrixP0=[[0.1]],
            discFactor=[1],
            chiPrompt=[1],
            chiDelayed=[1],
            inverseVelocity=[7.1e-08],
            diffusionCoefficient=[0.02],
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

neutronicsSolver.nuclearData.add_state(refState)


#==============================================================================*
# Settings

model = ffn.Case(timeFolders=[timeFolder0])

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
    print(f"keff = {model.keff()}")

    model.plot_mesh(region=nMesh)
    for facename, _ in wedgeFaces:
        model.plot_boundary(region=nMesh, boundaryName=facename)

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

    X, Y, Z = fluidfoam.readmesh('.', region=nMesh.region, structured=False)
    flux0 = fluidfoam.readscalar(
        '.',
        time_name=str(settings.endTime),
        name='flux0',
        region=nMesh.region,
        structured=False
    )

    # Calculate radial coordinate from Cartesian coordinates
    r_mesh = np.sqrt(X**2 + Y**2 + Z**2)

    # Normalize the flux
    normFlux = [flux0_ / max(flux0) for flux0_ in flux0]

    # Analytic solution
    fluxTh = lambda r_: np.sin(r_ * np.pi/dr)/r_/np.pi*dr

    # Analytical flux at mesh points
    analytical_flux = fluxTh(r_mesh)

    # Threshold for relative error calculation (avoid division by near-zero values)
    threshold = 0.05  # Only plot relative error where flux > 5% of maximum
    max_analytical_flux = np.max(analytical_flux)
    mask = analytical_flux > threshold * max_analytical_flux

    # Compute the relative error between GeN-Foam and the analytic solution
    relError = np.array([(normFlux_ / analytical_flux_ - 1) * 100
                          for normFlux_, analytical_flux_ in zip(normFlux, analytical_flux)])

    # Apply mask to data for plotting
    r_masked = r_mesh[mask]
    relError_masked = relError[mask]


    fig, axes = plt.subplots(2, 1, figsize=(5, 6))
    axFlux, axError = axes.flatten()

    axFlux.plot(r_mesh, normFlux, label="GeN-Foam")
    axFlux.plot(r_mesh, analytical_flux, label="Analytic", ls="--")
    axFlux.set_ylabel("Neutron flux [a.u]")
    axFlux.legend()

    axError.plot(r_masked, relError_masked)
    axError.set_ylabel("Relative error [%]")

    for ax in axes:
        ax.set_xlabel("Radial position [m]")
        ax.set_xlim((0, dr))

    fig.tight_layout()
    fig.savefig("fig_results_fluxDistribution.png")
    plt.close()


#==============================================================================*
