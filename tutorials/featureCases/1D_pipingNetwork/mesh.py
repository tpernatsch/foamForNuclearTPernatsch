import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import foamForNuclear.boundaryConditions as bc
import math
import matplotlib.pyplot as plt

# === Secondary mesh creation ===

refinementFactor = 12

secondaryMesh = mesh.BlockMesh(region='fluid')

# === Constants ===

pipe1Diameter = 1
pipe2Diameter = 0.8
pipe3Diameter = 1.1
pipe4Diameter = 0.9



# === Create main pipe ===

pipe1 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe1",
    originPosition=ffn.Vector(0, 0, 0),
    direction=ffn.Vector(1, 0, 0),
    length=10,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=10 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True
)


pipe2 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe2",
    originPosition=pipe1,
    direction=ffn.Vector(0, 0, 1),
    length=10,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=10 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
    isCustomNames=True,
    outletCustomName="pipe1topipe2_outletAMI"
)
pipe3 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe3",
    originPosition=pipe2,
    direction=ffn.Vector(1, 0, 0),
    length=10,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=10 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
)
pipe4 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe4",
    originPosition=pipe3,
    direction=ffn.Vector(0, 0, 1),
    length=5,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=5 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
    isCustomNames=True,
    outletCustomName="pipe3topipe4_outletAMI"
)
pipe5 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe5",
    originPosition=pipe4,
    direction=ffn.Vector(1, 0, 0),
    length=20,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=20 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
)



pipe6 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe6",
    originPosition=pipe3,
    direction=ffn.Vector(0, 0, -1),
    length=5,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=5 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
    isCustomNames=True,
    outletCustomName="pipe3topipe6_outletAMI"
)
pipe7 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe7",
    originPosition=pipe6,
    direction=ffn.Vector(1, 0, 0),
    length=20,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=20 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
)


pipe8 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe8",   
    originPosition=pipe5,
    direction=ffn.Vector(0, 0, -1),
    length=5,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=5 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
)

pipe9 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe9",   
    originPosition=pipe7,
    direction=ffn.Vector(0, 0, 1),
    length=5,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=5 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
)

pipe10 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe10",   
    originPosition=pipe8,
    direction=ffn.Vector(1, 0, 0),
    length=10,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=10 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
    isCustomNames=True,
    inletCustomName="pipe10frompipe8_inletAMI"
)

secondaryMesh.connect_pipes_custom_names(
    pipe9,
    pipe10,
    1,
    customOutletName="pipe9_outletAMI",
    customInletName="pipe10frompipe9_inletAMI"
)

pipe11 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe11",   
    originPosition=pipe10,
    direction=ffn.Vector(0, 0, -1),
    length=10,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=10 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True
)

pipe12 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe12",   
    originPosition=pipe1,
    direction=ffn.Vector(0, 0, -1),
    length=10,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=10 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
    isCustomNames=True,
    outletCustomName="pipe1topipe12_outletAMI"
)

pipe13 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe13",   
    originPosition=pipe12,
    direction=ffn.Vector(1, 0, 0),
    length=44,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=44 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True
)

pipe14 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe14",   
    originPosition=pipe13,
    direction=ffn.Vector(0, 0, 1),
    length=10,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=10 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True
)

pipe15 = secondaryMesh.add_pipe_1D_from_direction(
    name="pipe15",   
    originPosition=pipe14,
    direction=ffn.Vector(1, 0, 0),
    length=10,
    equivalentHydraulicDiameter=pipe1Diameter,
    n=10 * refinementFactor,
    elbowRadius=1,
    isAddBoundaryConditions=True,
    isCustomNames=True,
    inletCustomName="pipe15frompipe14_inletAMI"
)

secondaryMesh.connect_pipes_custom_names(
    pipe11,
    pipe15,
    1,
    customOutletName="pipe11_outletAMI",
    customInletName="pipe15frompipe11_inletAMI"
)

inlet = ffn.Face("inlet", boundaryType="patch")
outlet = ffn.Face("outlet", boundaryType="patch")

inlet.add_sub_face(pipe1.bottomFace())
outlet.add_sub_face(pipe15.topFace())

secondaryMesh.faces.append(inlet)
secondaryMesh.faces.append(outlet)

secondaryMesh.export_to_openfoam()
