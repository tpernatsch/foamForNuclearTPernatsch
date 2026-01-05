import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


thMesh = mesh.BlockMesh(region='fluidRegion')

if (True):
    outerRadius = 2

    frontBlock, rightBlock, backBlock, leftBlock = thMesh.create_ring_along_z(
        name="core",
        innerRadius=1,
        outerRadius=outerRadius,
        lowZ=0,
        highZ=1,
        nr=2,
        nt=30,
        nz=4,
    )

    # Extrude out concentric
    for block, facename, dr in [
        (frontBlock, "front", 0.5),
        (rightBlock, "right", 1),
        (backBlock, "back", 1.5),
        (leftBlock, "left", 2)
    ]:
        thMesh.extrude_normal_ring_section(
            targetBlock=block,
            facename=facename,
            name=f"extrude_out_{facename}",
            r=outerRadius + dr,
            xCenter=0,
            yCenter=0,
            nr=4
        )

    # Extrude in with eccentricity
    for block, facename in [
        (frontBlock, "back"),
        (rightBlock, "left"),
        (backBlock, "front"),
        (leftBlock, "right")
    ]:
        thMesh.extrude_normal_ring_section(
            targetBlock=block,
            facename=facename,
            name=f"extrude_in_{facename}",
            r=0.25,
            xCenter=-0.5,
            yCenter=0.25,
            nr=3
        )


standaloneFaces = thMesh.get_standalone_faces_as_list_blocks()

print(len(standaloneFaces))

for block, facename, barycenter in standaloneFaces:
    if ("extrude_out_" not in block.name or facename in ['top', 'bottom']):
        continue

    normal = block.get_face_normal(faceName=facename)

    print(f"{block.name:18} {facename:6} {barycenter} {normal} {normal.dot(barycenter)}")

    if (normal.dot(barycenter) <= 0):
        continue

    thMesh.extrude_normal_ring_section(
        targetBlock=block,
        facename=facename,
        name=f"extrude_expand_{facename}",
        r=5,
        xCenter=0,
        yCenter=0,
        nr=3
    )

thMesh.isMergeCoincidentPoints = True


solver = ffn.ThermalHydraulicsSolver(region=thMesh.region, mesh=thMesh, solver='onePhase')

model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=thMesh, show_edges=True)
for normal in ['x', 'y', 'z']:
    model.plot_mesh(region=thMesh, show_edges=True, normal=normal)

# model.plot_boundary(region=thMesh, boundaryName='defaultFaces', show_edges=True)
# model.plot_boundary(region=thMesh, boundaryName='topFuel', show_edges=True)
# model.plot_boundary(region=thMesh, boundaryName='bottomFuel', show_edges=True)
# model.plot_boundary(region=thMesh, boundaryName='outerClad', show_edges=True)
