"""

"""
#==============================================================================*
# Imports

import os
import sys

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# Mesh

nMesh = mesh.BlockMesh(region='neutroRegion')

side = 0.1
lCore = 1.1969
lRefl = 0.5

reflBot = nMesh.create_cube('reflector', 0, 0, -lCore/2-lRefl, side, side, -lCore/2, 1, 1, 150)
core = nMesh.extrude_top(reflBot, 'core', lCore, 151)
reflTop = nMesh.extrude_top(core, 'reflector', lRefl, 150)

walls = ffn.mesh.Face("walls", boundaryType="wall")
for zone in [core, reflTop, reflBot]:
    walls.add_sub_face(zone.frontFace())
    walls.add_sub_face(zone.backFace())
    walls.add_sub_face(zone.leftFace())
    walls.add_sub_face(zone.rightFace())

top = ffn.mesh.Face("top")
top.add_sub_face(reflTop.topFace())

bottom = ffn.mesh.Face("bottom")
bottom.add_sub_face(reflBot.bottomFace())

nMesh.add_boundary(walls)
nMesh.add_boundary(top)
nMesh.add_boundary(bottom)


#==============================================================================*
# Fields

timeFolder0 = ffn.TimeFolder(0)

defaultFlux = ffn.fields.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.fields.Dimension(default='flux')
defaultFlux.internalField = 1e21
defaultFlux.set_boundary_condition("walls", bc.ZeroGradient())
defaultFlux.set_boundary_condition("top", bc.FixedValue(0))
defaultFlux.set_boundary_condition("bottom", bc.FixedValue(0))

timeFolder0.append(defaultFlux)


#==============================================================================*
# Solvers

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    region=nMesh.region,
    solver="diffusionNeutronics",
    mesh=nMesh
)

neutronicsSolver.fvSchemes.ddtSchemes['default'] = 'steadyState'

nuclearData = neutronicsSolver.nuclearData

refState = ffn.nuclearData.NuclearDataState(
    name="reference",
    zones=[
        ffn.nuclearData.NuclearDataZone(
            "core",
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
            decayConstant=[0.0125371, 0.0300828, 0.109879, 0.325484, 1.3036, 9.51817],
            delayedFraction=[7.2315e-05, 0.000609661, 0.000471181, 0.00118907, 0.000445487, 9.58515e-05],
        ),
        ffn.nuclearData.NuclearDataZone(
            "reflector",
            fuelFraction=0,
            removalXS=[3],
            nuFissionXS=[0],
            powerXS=[0],
            scatteringMatrixP0=[[0.1]],
            discFactor=[1],
            chiPrompt=[0],
            chiDelayed=[0],
            inverseVelocity=[7.1e-08],
            diffusionCoefficient=[0.2],
            integralFlux=[1],
            decayConstant=[0.0125371, 0.0300828, 0.109879, 0.325484, 1.3036, 9.51817],
            delayedFraction=[7.2315e-05, 0.000609661, 0.000471181, 0.00118907, 0.000445487, 9.58515e-05],
        )
    ]
)

nuclearData.add_state(refState)


idxField = neutronicsSolver.create_zone_field(nMesh.cellZones)
timeFolder0.append(idxField)


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


#==============================================================================*
# Run

ffn.run(model, is_preprocessing=True)

model.plot_mesh(nMesh, fieldName=idxField.name)


#==============================================================================*
# Post-processing

keff = model.keff()

print(f"keff = {keff}")


# Only executed when >>> python3 Allrun.py
if __name__ == "__main__":

    model.plot_mesh(region=nMesh.region)
    model.plot_boundary(region=nMesh.region, boundaryName=walls.name)

    model.plot_slice(
        region=nMesh.region,
        time=settings.endTime,
        fieldName="flux0",
        show_edges=False,
        cmap='Blues_r',
        unit="n/m2/s"
    )

    model.plot_mesh(
        region=nMesh.region,
        time=settings.endTime,
        fieldName="flux0",
        cmap="Blues_r",
        show_edges=False
    )


#==============================================================================*
