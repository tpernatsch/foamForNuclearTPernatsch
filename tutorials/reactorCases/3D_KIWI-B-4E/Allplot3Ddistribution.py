"""
Data taken from LA-3185-MS, Chapiter X, Figure 1 to 8, at time 20800 s where
the core power was at 905 MW and mass flow at 68.7 lb/s (page 110)
"""

#==============================================================================*
# Imports

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colormaps
from matplotlib.colors import Normalize, TwoSlopeNorm
import scipy as sp
from scipy.ndimage import gaussian_filter
import sys
from fluidfoam import readmesh, readscalar


#==============================================================================*
# Usage

if (len(sys.argv) != 3):
    print(f"\n    Usage: python3 {sys.argv[0]} path/to/case timestep\n")
    sys.exit(1)


#==============================================================================*
# Units and constants

inch = 2.54e-2 # m

coreHeight = 52 * inch


#==============================================================================*
# Extract data from GeN-Foam

caseName = sys.argv[1] # steadyStateTH
region = 'fluidRegion'
timeStep = sys.argv[2] # 6

x, y, z = readmesh(caseName, region=region, structured=False)

Tmatrix = readscalar(caseName, timeStep, 'Tmatrix.lumpedNuclearStructure', region=region, structured=False)

r = [np.sqrt(x_**2+y_**2) for x_, y_ in zip(x, y)]

rmin, rmax, zmin, zmax = 0, 17*inch, 0, coreHeight

rMesh = np.linspace(rmin, rmax, 12)
zMesh = np.linspace(zmin, zmax, 12)

TmatrixMesh = np.zeros((len(zMesh)-1, len(rMesh)-1))
TmatrixMeshCount = np.zeros((len(zMesh)-1, len(rMesh)-1))

# Radial average
for k, Tmatrix_ in enumerate(Tmatrix[:len(r)]):
    if (Tmatrix_ > 100):
        for i, ri in enumerate(rMesh[:-1]):
            if (ri <= r[k] and r[k] <= rMesh[i+1]):
                for j, zj in enumerate(zMesh[:-1]):
                    if (zj-coreHeight/2 <= z[k] and z[k] <= zMesh[j+1]-coreHeight/2):
                        TmatrixMesh[j][i] += Tmatrix_
                        TmatrixMeshCount[j][i] += 1
                        break

for i in range(len(TmatrixMesh)):
    for j in range(len(TmatrixMesh[i])):
        if (TmatrixMeshCount[i][j] > 0):
            TmatrixMesh[i][j] = TmatrixMesh[i][j]/TmatrixMeshCount[i][j]


#==============================================================================*
# Extract data from report LA-3185-MS, Chapiter X, Figure 1 to 8

# Station 2 (2 inch above core inlet)
r2 = [
    0.0420495196240371, 0.0936050570545152, 0.144659927848857,
    0.195524044143424, 0.247291303726467, 0.297430757484358, 0.348268110189592,
    0.397492999254643
]
T2 = [
    414.017455464, 425.090430304428, 380.855906596721, 394.881885280788,
    377.930274624716, 346.974266952597, 322.406433014241, 305.958043690525
]
z2 = [2*inch]*len(r2)

# Station 8 (8 inch above core inlet)
r8 = [
    0.009592224983394, 0.0404547593288014, 0.0915866598068164, 0.142361055251745,
    0.192316625987591, 0.242858235697397, 0.293850358439483, 0.344765417410809,
    0.396425780482331
]
T8 = [
    1005.75249361382, 1000.18077548079, 988.554167958181, 944.82381863114,
    932.948445607337, 888.466055878403, 849.751634141868, 832.858282135301,
    804.927755981671
]
z8 = [8*inch]*len(r8)

# Station 14 (14 inch above core inlet)
r14 = [
    0.0438059474002053, 0.0938430171994903, 0.144748497780786,
    0.195066283938294, 0.245617586358393, 0.297307014392678, 0.348053505178166,
    0.399474637932025
]
T14 = [
    1305.52568053193, 1293.74575130708, 1271.6510099114, 1232.97049661335,
    1194.03641780308, 1182.99519218781, 1171.71139752386, 1189.33031951837
]
z14 = [14*inch]*len(r14)

# Station 20 (20 inch above core inlet)
r20 = [
    0.0105466529575585, 0.0410136469178339, 0.0914452577607398, 0.141113628533103,
    0.192349815131315, 0.242035135263957, 0.293955868474272, 0.344921211213324,
    0.395096333267512

]
T20 = [
    1638.8813164557, 1583.28394936709, 1577.57529113924, 1537.86288607595,
    1510.06420253165, 1466.87696202532, 1489.21518987342, 1421.70410126582,
    1420.95949367089
]
z20 = [20*inch]*len(r20)

# Station 26 (26 inch above core inlet)
r26 = [
    0.0423059723756073, 0.0938795427218828, 0.144335199145148,
    0.195230086898673, 0.246084693131405, 0.297067874723436, 0.347903560847917,
    0.398442425059082
]
T26 = [
    1865.67745905529, 1838.89435953247, 1849.78512832735, 1910.17466587602,
    1759.91060814296, 1776.98489494552, 1755.05039094377, 1721.28077827514
]
z26 = [26*inch]*len(r26)

# Station 30 (30 inch above core inlet)
r30 = [
    0.00989858622363774, 0.0397817999963739, 0.0605345150392525,
    0.090629105098499, 0.110502825336774, 0.141479019784408, 0.161454513790262,
    0.191596293822938, 0.211832528741884, 0.242076300007454, 0.262783999730057,
    0.293034729885717, 0.313419711080356, 0.345017203933916, 0.365464597674058,
    0.395089897582801, 0.416053553981559
]
T30 = [
    2055.46573425517, 2094.20555785746, 2050.06850882653, 2066.29006710704,
    2060.2377880665, 2038.09657308104, 2021.20216629298, 2004.61849867345,
    1987.72493126726, 1932.49796896058, 1976.49047643337, 2031.63417205848,
    1998.89441794436, 1882.7872564318, 1871.17659064345, 1881.83640464866,
    1870.78342417521
]
z30 = [30*inch]*len(r30)

# Station 38 (38 inch above core inlet)
r38 = [
    0.0433771178404587, 0.0946157493882506, 0.145360874056926, 0.195838182208934,
    0.246776232453445, 0.297454199484962, 0.348461442446545, 0.39868334806936
]
T38 = [
    2199.21252833578, 2089.29346027192, 2166.20741438902, 2166.77988907938,
    2160.89615633976, 2116.8429406918, 2084.07677680103, 1973.09256647576
]
z38 = [38*inch]*len(r38)

# Station 44 (44 inch above core inlet)
r44 = [
    0.0438154647667456, 0.0937768266169, 0.144183678540273, 0.194728142580504,
    0.245890694942216, 0.296324369566333, 0.347336481563006, 0.398626150202155
]
T44 = [
    2028.29107650729, 1984.27808176998, 1956.13238533593, 1955.89083881406,
    2111.58415206724, 2016.96114671777, 1967.20386359802, 1972.7077432955
]
z44 = [44*inch]*len(r44)

# Surface interpolation
rAllData, zAllData, TAllData = r2+r8+r14+r20+r26+r30+r38+r44, z2+z8+z14+z20+z26+z30+z38+z44, T2+T8+T14+T20+T26+T30+T38+T44
splineExp = sp.interpolate.Rbf(rAllData, zAllData, TAllData, function='thin-plate')


rInterp = [(ri+ro)/2 for ri, ro in zip(rMesh[:-1], rMesh[1:])]
zInterp = np.flip([(zi+zo)/2 for zi, zo in zip(zMesh[:-1], zMesh[1:])])

rInterp, zInterp = np.meshgrid(rInterp, zInterp)

TInterpExp = splineExp(rInterp, zInterp)

# Compute relative difference
diffT = np.matrix(TmatrixMesh) - np.matrix(TInterpExp)


#==============================================================================*
# Plot

# Create figures
fig, (ax, axDiff) = plt.subplots(ncols=2, figsize=(9, 5))

# Colormap
cmap = colormaps['inferno']
normalizer = Normalize(0, np.matrix(TmatrixMesh).max())
im = cm.ScalarMappable(norm=normalizer, cmap=cmap)

# Plot GeN-Foam data
ax.imshow(
    TmatrixMesh,
    interpolation='bilinear',
    cmap=cmap,
    extent=[rmin, rmax, zmin, zmax],
    label="GeN-Foam",
)
contour = ax.contour(
    rInterp,
    zInterp,
    gaussian_filter(TmatrixMesh, 0.7),
    cmap=colormaps['binary'],
    levels=[400, 600, 800, 1000, 1200, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500]
)
ax.clabel(contour, inline=True)

# Plot experimental data
ax.scatter(
    rAllData,
    zAllData,
    c=TAllData,
    label="Experimental data",
    cmap=cmap,
    norm=normalizer,
    edgecolor='black',
    linewidth=1
)

# Plot difference
imDiff = axDiff.imshow(
    diffT,
    interpolation='gaussian',
    cmap=colormaps['RdBu_r'],
    extent=[rmin, rmax, zmin, zmax],
    norm=TwoSlopeNorm(vcenter=0)
)
contourDiff = axDiff.contour(
    rInterp, zInterp, gaussian_filter(diffT, 0.7),
    colors="black"
)
contourDiff.set_linestyle('solid')
axDiff.clabel(contourDiff, inline=True)

# Legend
ax.legend(
    fontsize=10,
    bbox_to_anchor=(0, 1.01, 1, 0.2),
    loc="lower left",
    mode='expand',
    ncol=2
)

# Labels
ax.set_xlabel('Radius [m]')
ax.set_ylabel('Axial [m]')
axDiff.set_xlabel('Radius [m]')
axDiff.set_ylabel('Axial [m]')
fig.colorbar(im, ax=ax, label="Fuel temperature [K]")
fig.colorbar(imDiff, ax=axDiff, label="Temperature difference (GeN-Foam - Exp.) [K]")

# Aspect
ax.set_aspect(0.5)
axDiff.set_aspect(0.5)
fig.tight_layout()

# Save the figure
fig.savefig("fig_results_2dMapFuelTemperatureRZ.png")


#==============================================================================*

fig3D = plt.figure()
ax3D = fig3D.add_subplot(projection='3d')

ax3D.scatter(rAllData, zAllData, TAllData, label="Experimental data")

surf = ax3D.plot_surface(rInterp, zInterp, TmatrixMesh, cmap=cmap)

ax3D.set_xlabel("Radius [m]")
ax3D.set_ylabel("Axial [m]")
ax3D.set_zlabel("Temperature [K]")

fig3D.colorbar(surf)
fig3D.savefig("fig_results_3dTfuel.png")


fig3Derr = plt.figure()
ax3Derr = fig3Derr.add_subplot(projection='3d')

surfErr = ax3Derr.plot_surface(rInterp, zInterp, diffT, cmap=cmap)

ax3Derr.set_xlabel("Radius [m]")
ax3Derr.set_ylabel("Axial [m]")
ax3Derr.set_zlabel("Temperature difference [K]")

fig3Derr.colorbar(surfErr)
plt.close()


#==============================================================================*
