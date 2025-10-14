"""

"""
#==============================================================================*
# Imports

import matplotlib.pyplot as plt
import numpy as np

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh

#==============================================================================*

ffn.allclean()

#==============================================================================*
# Mesh

nMesh = mesh.BlockMesh(region='neutroRegion')

fuelLength = 1.5
nz = 150

zone0 = nMesh.createCube('zone0', -0.05, -0.05, 0, 0.05, 0.05, fuelLength, 1, 1, nz)

walls = ffn.Face("walls", boundaryType="wall")
walls.addSubFace(zone0.frontFace())
walls.addSubFace(zone0.backFace())
walls.addSubFace(zone0.leftFace())
walls.addSubFace(zone0.rightFace())

top = ffn.Face("top")
top.addSubFace(zone0.topFace())

bottom = ffn.Face("bottom")
bottom.addSubFace(zone0.bottomFace())

nMesh.addBoundary(walls)
nMesh.addBoundary(top)
nMesh.addBoundary(bottom)


#==============================================================================*
# Fields

timeFolder0 = ffn.TimeFolder(0)

defaultFlux = ffn.Field("defaultFlux", region="neutroRegion")
defaultFlux.dimensions = ffn.Dimension(default='flux')
defaultFlux.internalField = 1e21
defaultFlux.set_boundary_condition("walls", bc.ZeroGradient())
defaultFlux.set_boundary_condition("top", bc.FixedValue(0))
defaultFlux.set_boundary_condition("bottom", bc.FixedValue(0))

timeFolder0.append(defaultFlux)


#==============================================================================*
# Solvers

neutronicsSolver = ffn.NeutronicsSolver(
    region="neutroRegion",
    solver="diffusionNeutronics",
)

print(f"List of require fields: {neutronicsSolver.get_required_fields()}")

neutronicsSolver.mesh = nMesh

neutronicsSolver.fvSchemes.ddtSchemes['default'] = 'steadyState'

nuclearData = neutronicsSolver.nuclearData

refState = ffn.NuclearDataState(
    name="reference",
    zones=[
        ffn.NuclearDataZone(
            "zone0",
            fuelFraction=0.4,
            removalXS=[5.0082],
            nuFissionXS=[5.5675],
            powerXS=[1],
            scatteringMatrixP0=[[0.1]],
            discFactor=[1],
            chiPrompt=[1],
            chiDelayed=[1],
            inverseVelocity=[7.1e-08],
            diffusionCoefficient=[0.1275],
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

model = ffn.Model()

model.solvers.append(neutronicsSolver)
model.timeFolders = [timeFolder0]

settings: ffn.ControlDict = model.settings

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


#==============================================================================*
# Run

ffn.run(model=model, is_preprocessing=True)


#==============================================================================*
# Post-processing

print(f"keff = {model.keff()}")

model.plot_mesh(region="neutroRegion")
model.plot_boundary(region="neutroRegion", boundaryName=walls.name)

model.plot_slice(
    region="neutroRegion",
    time=settings.endTime,
    fieldName="flux0",
    show_edges=False,
    cmap='Blues',
    unit="n/m2/s"
)

model.plot_mesh(
    region="neutroRegion",
    time=settings.endTime,
    fieldName="flux0",
    cmap="Blues",
    show_edges=False
)

points, flux0 = model.sample_over_line(
    region="neutroRegion",
    time=settings.endTime,
    fieldName="flux0",
    point1=(0, 0, fuelLength/nz * 0.5),
    point2=(0, 0, fuelLength/nz * (nz-0.5)),
)

z = points[:,2]

# Normalize the flux
normFlux = [flux0_ / max(flux0) for flux0_ in flux0]

# Analytic solution
fluxTh = lambda z_: np.sin(z_ * np.pi/fuelLength)

# Compute the relative error between GeN-Foam and the analytic solution
relError = [(normFlux_ / fluxTh(z_) - 1) * 100 for z_, normFlux_ in zip(z, normFlux)]


fig, axes = plt.subplots(2, 1, figsize=(5, 6))
axFlux, axError = axes.flatten()

axFlux.plot(z, normFlux, label="GeN-Foam")
axFlux.plot(z, fluxTh(z), label="Analytic", ls="--")
axFlux.set_ylabel("Neutron flux [a.u]")
axFlux.legend()

axError.plot(z, relError)
axError.set_ylabel("Relative error [%]")

for ax in axes:
    ax.set_xlabel("Axial position [m]")
    ax.set_xlim((0, fuelLength))

fig.tight_layout()
fig.savefig("fig_results_fluxDistribution.png")
plt.close()


#==============================================================================*
