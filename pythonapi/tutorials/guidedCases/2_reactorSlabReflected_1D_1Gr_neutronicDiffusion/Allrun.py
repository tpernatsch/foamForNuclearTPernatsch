"""

"""
#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh

#==============================================================================*

ffn.allclean()

#==============================================================================*
# Mesh

nMesh = mesh.BlockMesh(region='neutroRegion')

side = 0.1
lCore = 1.197
lRefl = 0.5

reflBot = nMesh.createCube('reflector', 0, 0, -lCore/2-lRefl, side, side, -lCore/2, 1, 1, 150)
core = nMesh.extrudeTop([reflBot], 'core', lCore, 151)
reflTop = nMesh.extrudeTop([core], 'reflector', lRefl, 150)

walls = ffn.Face("walls", boundaryType="wall")
for zone in [core, reflTop, reflBot]:
    walls.addSubFace(zone.frontFace())
    walls.addSubFace(zone.backFace())
    walls.addSubFace(zone.leftFace())
    walls.addSubFace(zone.rightFace())

top = ffn.Face("top")
top.addSubFace(reflTop.topFace())

bottom = ffn.Face("bottom")
bottom.addSubFace(reflBot.bottomFace())

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
    solver="diffusionNeutronics"
)
neutronicsSolver.mesh = nMesh

neutronicsSolver.fvSchemes.ddtSchemes['default'] = 'steadyState'

nuclearData = neutronicsSolver.nuclearData

refState = ffn.NuclearDataState(
    name="reference",
    zones=[
        ffn.NuclearDataZone(
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
        ffn.NuclearDataZone(
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

model.plot_mesh(nMesh, fieldName=idxField.name)


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


#==============================================================================*
