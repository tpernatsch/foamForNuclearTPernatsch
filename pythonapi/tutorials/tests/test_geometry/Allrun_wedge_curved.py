import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh

inch = 2.54

nozzleInletRadius  = 45.39209016393443 # 38.250/2.0 * inch #
nozzleOutletRadius = 30.370/2.0 * inch
nozzleInletCoreBotSupport  =  62.694 * inch - 53.5 * inch
nozzleMinCoreBotSupport    =  82.250 * inch - 53.5 * inch
nozzleOutletCoreBotSupport = 115.112 * inch - 53.5 * inch
bot = nozzleOutletCoreBotSupport + 1500.0
maxRadius = 3*nozzleOutletRadius


nMesh = mesh.BlockMeshWedge(appertureAngle=20, region='neutroMesh', scale=0.01)

block = nMesh.create_wedge_conical(
    "block",
    innerRadiusBottom=0,
    outerRadiusBottom=nozzleInletRadius,
    innerRadiusTop=0,
    outerRadiusTop=nozzleOutletRadius,
    lowZ=0,
    highZ=nozzleOutletCoreBotSupport,
    nr=10,
    nz=200
)


points = [
    # (45.39209016393443, 0),
    (45.39209016393443, 23.14354985915492),
    (45.12663934426229, 25.347697464788727),
    (44.59573770491804, 27.000808169014075),
    (43.53393442622951, 28.929437323943656),
    (42.472131147540985, 30.85806647887323),
    (16.192500000000003, 59.236466901408434),
    (14.86524590163934, 61.165096056338),
    (13.537991803278686, 63.36924366197181),
    (13.00709016393442, 65.29787281690139),
    (12.476188524590157, 67.22650197183097),
    (12.210737704918024, 69.15513112676054),
    (12.210737704918024, 71.63479718309857),
    (12.476188524590157, 74.66550014084504),
    (13.272540983606554, 77.42068464788731),
    (14.86524590163934, 80.17586915492956),
    (16.192500000000003, 82.93105366197182),
    (17.51975409836065, 85.96175661971827),
    (18.58155737704918, 88.4414226760563),
    (19.908811475409838, 91.74764408450702),
    (20.97061475409836, 94.22731014084503),
    (22.29786885245901, 97.53353154929574),
    (23.625122950819677, 100.56423450704224),
    (24.686926229508195, 103.5949374647887),
    (26.545081967213108, 108.27875112676053),
    (27.87233606557377, 112.41152788732391),
    (29.46504098360656, 116.5443046478873),
    (30.79229508196721, 121.22811830985913),
    (32.119549180327866, 125.63641352112673),
    (33.181352459016395, 129.49367183098587),
    (33.977704918032785, 133.350930140845),
    (35.039508196721314, 137.20818845070417),
    (35.57040983606557, 140.51440985915488),
    (36.36676229508196, 145.47374197183095),
    (37.16311475409836, 149.8820371830985),
    # (nozzleOutletRadius, nozzleOutletCoreBotSupport),
]

nMesh.add_right_face_edge_polyline(
    block=block,
    rzCoords=points
)
# atmos1 = nMesh.extrude_top([block], "atmosphere", dz=bot, nz=1500)
# atmos2 = nMesh.extrude_right([atmos1], "atmosphere", dr=maxRadius, nr=30)

solver = ffn.solvers.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.case.Case()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
model.plot_mesh(region=nMesh, show_edges=True, normal='y')
model.plot_mesh(region=nMesh, show_edges=True, normal='z')
