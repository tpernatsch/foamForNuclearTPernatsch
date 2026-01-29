import numpy as np
import matplotlib.pyplot as plt
import json

with open("./rod3Ddict.json","r") as f:
    inputDict = json.load(f)

unit_str = inputDict["unit"]
if (unit_str=="mm"):
    c2m = 0.001 
elif (unit_str=="cm"): 
    c2m = 0.01
elif (unit_str=="dm"): 
    c2m = 0.1
elif (unit_str=="m"): 
    c2m = 1
else:
    print("ERROR in input: unit can be 'm', 'dm', 'cm' or 'mm'")
    quit()

rif = inputDict["fuelInnerRadius"][0]
rof = inputDict["fuelOuterRadius"][0]
ric = inputDict["cladInnerRadius"][0]
roc = inputDict["cladOuterRadius"][0]
hei = inputDict["rodHeight"][0]
nRadF = inputDict["nRadialCellsFuel"][0]
nRadC = inputDict["nRadialCellsClad"][0]
nTheta = inputDict["nAzimuthalCellsPerQuarter"][0]
nZeta = inputDict["nZetaCells"][0]
c = inputDict["center"]

# Compute vertices
v=[
    [ c[0]+rif, c[1], c[2] ],
    [ c[0]+rof, c[1], c[2] ],
    [ c[0], c[1]+rof, c[2] ],
    [ c[0], c[1]+rif, c[2] ],
    [ c[0]+rif, c[1], c[2]+hei ],
    [ c[0]+rof, c[1], c[2]+hei ],
    [ c[0], c[1]+rof, c[2]+hei ],
    [ c[0], c[1]+rif, c[2]+hei ],
    
    [ c[0]-rof, c[1], c[2] ],
    [ c[0]-rif, c[1], c[2] ],
    [ c[0]-rof, c[1], c[2]+hei ],
    [ c[0]-rif, c[1], c[2]+hei ],

    [ c[0], c[1]-rof, c[2] ],
    [ c[0], c[1]-rif, c[2] ],
    [ c[0], c[1]-rof, c[2]+hei ],
    [ c[0], c[1]-rif, c[2]+hei ],

    # Clad
    [ c[0]+ric, c[1], c[2] ],
    [ c[0]+roc, c[1], c[2] ],
    [ c[0], c[1]+roc, c[2] ],
    [ c[0], c[1]+ric, c[2] ],
    [ c[0]+ric, c[1], c[2]+hei ],
    [ c[0]+roc, c[1], c[2]+hei ],
    [ c[0], c[1]+roc, c[2]+hei ],
    [ c[0], c[1]+ric, c[2]+hei ],
    
    [ c[0]-roc, c[1], c[2] ],
    [ c[0]-ric, c[1], c[2] ],
    [ c[0]-roc, c[1], c[2]+hei ],
    [ c[0]-ric, c[1], c[2]+hei ],

    [ c[0], c[1]-roc, c[2] ],
    [ c[0], c[1]-ric, c[2] ],
    [ c[0], c[1]-roc, c[2]+hei ],
    [ c[0], c[1]-ric, c[2]+hei ]
]

# Initialize an empty list to store the rotated points
edge_v = []

# Iterate through each point in v and apply the rotation
for point in v:
    
    # Define the angle in radians
    angle = np.deg2rad(45)

    # Define the rotation matrix for a 3D rotation around the z-axis
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0],
        [0, 0, 1]
    ])
    
    # Convert the point to a NumPy array
    point_array = np.array(point)
    
    # Translate the point by -c to rotate around the origin
    translated_point = point_array - c
    
    # Apply the rotation matrix
    rotated_point = np.dot(rotation_matrix, translated_point)
    
    # Translate the rotated point back to the original position
    final_rotated_point = rotated_point + c
    
    # Append the rotated point to the new list
    edge_v.append(final_rotated_point.tolist())

edges =[
    [1	, 2],
    [2	, 8],
    [8	, 12],
    [12	, 1],
    [5	, 6],
    [6	, 10],
    [10	, 14],
    [14	, 5],
    [0	, 3],
    [3	, 9],
    [9	, 13],
    [13	, 0],
    [4	, 7],
    [7	, 11],
    [11	, 15],
    [15	, 4],

    [17	, 18],
    [18	, 24],
    [24	, 28],
    [28	, 17],
    [21	, 22],
    [22	, 26],
    [26	, 30],
    [30	, 21],
    [16	, 19],
    [19	, 25],
    [25	, 29],
    [29	, 16],
    [20	, 23],
    [23	, 27],
    [27	, 31],
    [31	, 20]
]

edges_str = []
for i,e in enumerate(edges):
    x = edge_v[int(e[0])][0]
    y = edge_v[int(e[0])][1]
    z = edge_v[int(e[0])][2]
    edges_str.append(f"arc {int(e[0])} {int(e[1])} ({x} {y} {z})\n")

head="""
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     4.0                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

blocks=f"""
blocks
(
    hex (0 1 2 3 4 5 6 7)      fuel       ({nRadF} {nTheta} {nZeta}) simpleGrading (1 1 1)  
    hex (3 2 8 9 7 6 10 11)     fuel       ({nRadF} {nTheta} {nZeta}) simpleGrading (1 1 1)  
    hex (9 8 12 13 11 10 14 15)     fuel       ({nRadF} {nTheta} {nZeta}) simpleGrading (1 1 1)  
    hex (1 0 13 12 5 4 15 14)     fuel       ({nRadF} {nTheta} {nZeta}) simpleGrading (1 1 1)  
    hex (16 17 18 19 20 21 22 23)     cladding   ({nRadC} {nTheta} {nZeta}) simpleGrading (1 1 1)  
    hex (19 18 24 25 23 22 26 27)     cladding   ({nRadC} {nTheta} {nZeta}) simpleGrading (1 1 1)  
    hex (25 24 28 29 27 26 30 31)     cladding   ({nRadC} {nTheta} {nZeta}) simpleGrading (1 1 1)  
    hex (17 16 29 28 21 20 31 30)     cladding   ({nRadC} {nTheta} {nZeta}) simpleGrading (1 1 1)  
);
"""

footer="""
// ************************************************************************* //
"""

patches = [
    "fuelInnerSurface",
    "fuelOuterSurface",
    "topFuel",
    "bottomFuel",
    "cladInnerSurface",
    "cladOuterSurface",
    "topClad",
    "bottomClad",
]

pp = [
    [[0,4,7,3],[3,7,11,9],[9,11,15,13],[13,15,4,0]],
    [[1,2,6,5],[2,8,10,6],[8,12,14,10],[12,5,1,14]],
    [[4,5,6,7],[7,6,10,11],[11,10,14,15],[14,5,4,15]],
    [[0,3,2,1],[8,2,3,9],[12,8,9,13],[12,13,0,1]],
    [[16, 20, 23, 19], [19, 23, 27, 25], [25, 27, 31, 29], [29, 31, 20, 16]], 
    [[17, 18, 22, 21], [18, 24, 26, 22], [24, 28, 30, 26], [28, 21, 17, 30]], 
    [[20, 21, 22, 23], [23, 22, 26, 27], [27, 26, 30, 31], [30, 21, 20, 31]], 
    [[16, 19, 18, 17], [24, 18, 19, 25], [28, 24, 25, 29], [28, 29, 16, 17]]
]

ppp=""
for i, p in enumerate(patches):
    ppp += "    " + p + "    \n    {"
    ppp+="""
        type patch;
        faces 
        (
            ({:.0f} {:.0f} {:.0f} {:.0f})
            ({:.0f} {:.0f} {:.0f} {:.0f})
            ({:.0f} {:.0f} {:.0f} {:.0f})
            ({:.0f} {:.0f} {:.0f} {:.0f})
        );""".format( pp[i][0][0],pp[i][0][1],pp[i][0][2],pp[i][0][3],
                pp[i][1][0],pp[i][1][1],pp[i][1][2],pp[i][1][3],
                pp[i][2][0],pp[i][2][1],pp[i][2][2],pp[i][2][3],
                pp[i][3][0],pp[i][3][1],pp[i][3][2],pp[i][3][3])
    ppp+= "\n    }\n"

with open("blockMeshDict", "w") as f:
    f.write(head)

    f.write("\n")
    
    f.write(f"convertToMeters {c2m};")
    
    f.write("\n\n")
    
    f.write("vertices\n(\n")
    for i, vi in enumerate(v):
        f.write(f"    ({vi[0]} {vi[1]} {vi[2]}) // {i}\n")
    f.write(");")
    
    f.write("\n\n")
    
    f.write(blocks)
    
    f.write("\n\n")
    
    f.write("edges\n(\n")
    for i, ei in enumerate(edges_str):
        f.write(f"    {ei}")
    f.write(");\n\n")

    f.write("boundary\n(\n")
    f.write(ppp)
    f.write(");\n\n")

    f.write("mergePatchPairs\n(\n")
    f.write(");\n\n")
    
    f.write(footer)

