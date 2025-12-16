######### This script helps you in the creation of a blockMeshDict for 
######### a 1D or 2D (r-z). The 2D rod can have discrete pellets or a smeared
######### column. The script reads the input from a file named 
######### 'rodDict', which should be placed in the same folder
######### as this script. Read the comment section on the input file for more
######### information about how the rod is modeled.

import math
from collections import defaultdict
import os
import re
import sys

# importing the module 
import ast 

def writeHeader(file):

    file.write("/*--------------------------------*- C++ -*----------------------------------*\\" + "\n") 
    file.write("| ========                 |                                                 |" + "\n") 
    file.write("| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |" + "\n") 
    file.write("|  \    /   O peration     | Version:  5.0                                   |" + "\n") 
    file.write("|   \  /    A nd           | Web:      www.OpenFOAM.org                      |" + "\n") 
    file.write("|    \/     M anipulation  |                                                 |" + "\n")
    file.write("\*---------------------------------------------------------------------------*/" + "\n") 
    file.write("FoamFile" + "\n")
    file.write("{" + "\n")
    file.write("    version     5.0;" + "\n")
    file.write("    format      ascii;" + "\n")
    file.write("    class       dictionary;" + "\n")
    file.write("    object      blockMeshDict;" + "\n")
    file.write("}" + "\n")
    file.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" + "\n")

def writeVertices(wedgeAngle, rInner, rOuter, height, offset, file):
    
    #### Prepare input data ####
    wedgeAngleDegree = float(wedgeAngle)
    wedgeAngleRadiant = wedgeAngleDegree/180*math.pi
    rInner = float(rInner)
    rOuter = float(rOuter)
    height = float(height)
    offset = float(offset)

    # Prepare list of points
    pointList = []

    # Correction for wedge volume
    # Relevant only for large wedge angles
    sinAngle = math.sin(wedgeAngleRadiant/2)
    cosAngle = math.cos(wedgeAngleRadiant/2)
    correction = math.sqrt(sinAngle/(wedgeAngleRadiant/2))

    # Coordinates for vertices
    rInnerX = rInner*cosAngle*correction
    rInnerY = rInner*sinAngle*correction
    rOuterX = rOuter*cosAngle*correction
    rOuterY = rOuter*sinAngle*correction

    # Define block vertices for each block        
    pointList.append([rInnerX, -rInnerY, offset])
    pointList.append([rOuterX, -rOuterY, offset])
    pointList.append([rOuterX,  rOuterY, offset])
    pointList.append([rInnerX,  rInnerY, offset])

    pointList.append([rInnerX, -rInnerY, offset + height])
    pointList.append([rOuterX, -rOuterY, offset + height])
    pointList.append([rOuterX,  rOuterY, offset + height])
    pointList.append([rInnerX,  rInnerY, offset + height])

    # Write fuel vertices
    numberPoints = len(pointList)
    for i in range(numberPoints):
        newLine = "    (" 
        newLine += str(pointList[i][0]) + " "
        newLine += str(pointList[i][1]) + " " 
        newLine += str(pointList[i][2]) + ")" 
        newLine += "\n"
    
        file.write(newLine)

    file.write("\n")

def writePelletVertices(pelletType, wedgeAngle, rInner, rDish, hDish, rLand, hChamfer, rOuter, height, offset, file):
    
    #### Prepare input data ####
    wedgeAngleDegree = float(wedgeAngle)
    wedgeAngleRadiant = wedgeAngleDegree/180*math.pi
    rInner = float(rInner)
    rOuter = float(rOuter)
    height = float(height)
    offset = float(offset)

    # Prepare list of points
    pointList = []

    # Correction for wedge volume
    # Relevant only for large wedge angles
    sinAngle = math.sin(wedgeAngleRadiant/2)
    cosAngle = math.cos(wedgeAngleRadiant/2)
    correction = math.sqrt(sinAngle/(wedgeAngleRadiant/2))

    # Coordinates for vertices
    rInnerX = rInner*cosAngle*correction
    rInnerY = rInner*sinAngle*correction
    rDishX = rDish*cosAngle*correction
    rDishY = rDish*sinAngle*correction
    rLandX = rLand*cosAngle*correction
    rLandY = rLand*sinAngle*correction
    rOuterX = rOuter*cosAngle*correction
    rOuterY = rOuter*sinAngle*correction

    if pelletType == 'dishedChamfered':
        # Define block vertices for each block        
        pointList.append([rInnerX, -rInnerY, offset + hDish])
        pointList.append([rDishX, -rDishY, offset])
        pointList.append([rDishX,  rDishY, offset])
        pointList.append([rInnerX,  rInnerY, offset + hDish])

        pointList.append([rInnerX, -rInnerY, offset - hDish + height])
        pointList.append([rDishX, -rDishY, offset + height])
        pointList.append([rDishX,  rDishY, offset + height])
        pointList.append([rInnerX,  rInnerY, offset - hDish + height])

        pointList.append([rLandX, -rLandY, offset])
        pointList.append([rLandX,  rLandY, offset])

        pointList.append([rLandX, -rLandY, offset + height])
        pointList.append([rLandX,  rLandY, offset + height])

        pointList.append([rOuterX, -rOuterY, offset + hChamfer])
        pointList.append([rOuterX,  rOuterY, offset + hChamfer])

        pointList.append([rOuterX, -rOuterY, offset + height - hChamfer])
        pointList.append([rOuterX,  rOuterY, offset + height - hChamfer])

    elif pelletType == 'dished':
        # Define block vertices for each block        
        pointList.append([rInnerX, -rInnerY, offset + hDish])
        pointList.append([rDishX, -rDishY, offset])
        pointList.append([rDishX,  rDishY, offset])
        pointList.append([rInnerX,  rInnerY, offset + hDish])

        pointList.append([rInnerX, -rInnerY, offset - hDish + height])
        pointList.append([rDishX, -rDishY, offset + height])
        pointList.append([rDishX,  rDishY, offset + height])
        pointList.append([rInnerX,  rInnerY, offset - hDish + height])

        pointList.append([rOuterX, -rOuterY, offset])
        pointList.append([rOuterX,  rOuterY, offset])

        pointList.append([rOuterX, -rOuterY, offset + height])
        pointList.append([rOuterX,  rOuterY, offset + height])

    elif pelletType == 'chamfered':  
        # Define block vertices for each block        
        pointList.append([rInnerX, -rInnerY, offset])
        pointList.append([rLandX, -rLandY, offset])
        pointList.append([rLandX,  rLandY, offset])
        pointList.append([rInnerX,  rInnerY, offset])

        pointList.append([rInnerX, -rInnerY, offset + height])
        pointList.append([rLandX, -rLandY, offset + height])
        pointList.append([rLandX,  rLandY, offset + height])
        pointList.append([rInnerX,  rInnerY, offset + height])

        pointList.append([rOuterX, -rOuterY, offset + hChamfer])
        pointList.append([rOuterX,  rOuterY, offset + hChamfer])

        pointList.append([rOuterX, -rOuterY, offset + height - hChamfer])
        pointList.append([rOuterX,  rOuterY, offset + height - hChamfer])

    elif pelletType == 'flat':  
        # Define block vertices for each block        
        pointList.append([rInnerX, -rInnerY, offset])
        pointList.append([rOuterX, -rOuterY, offset])
        pointList.append([rOuterX,  rOuterY, offset])
        pointList.append([rInnerX,  rInnerY, offset])

        pointList.append([rInnerX, -rInnerY, offset + height])
        pointList.append([rOuterX, -rOuterY, offset + height])
        pointList.append([rOuterX,  rOuterY, offset + height])
        pointList.append([rInnerX,  rInnerY, offset + height])

    # Write fuel vertices
    numberPoints = len(pointList)
    for i in range(numberPoints):
        newLine = "    (" 
        newLine += str(pointList[i][0]) + " "
        newLine += str(pointList[i][1]) + " " 
        newLine += str(pointList[i][2]) + ")" 
        newLine += "\n"
    
        file.write(newLine)

    file.write("\n")

def writeCapVertices(wedgeAngle, rInner, rOuter, height, offset, file):
    
    #### Prepare input data ####
    wedgeAngleDegree = float(wedgeAngle)
    wedgeAngleRadiant = wedgeAngleDegree/180*math.pi
    rInner = float(rInner)
    rOuter = float(rOuter)
    height = float(height)
    offset = float(offset)

    # Prepare list of points
    pointList = []

    # Correction for wedge volume
    # Relevant only for large wedge angles
    sinAngle = math.sin(wedgeAngleRadiant/2)
    cosAngle = math.cos(wedgeAngleRadiant/2)
    correction = math.sqrt(sinAngle/(wedgeAngleRadiant/2))

    # Coordinates for vertices
    rInnerX = rInner*cosAngle*correction
    rInnerY = rInner*sinAngle*correction
    rOuterX = rOuter*cosAngle*correction
    rOuterY = rOuter*sinAngle*correction

    # Define block vertices for each block        
    pointList.append([0, 0, offset])
    pointList.append([rInnerX, -rInnerY, offset])
    pointList.append([rInnerX,  rInnerY, offset])
    pointList.append([0,  0, offset])

    pointList.append([0, 0, offset + height])
    pointList.append([rInnerX, -rInnerY, offset + height])
    pointList.append([rInnerX,  rInnerY, offset + height])
    pointList.append([0,  0, offset + height])

    # pointList.append([rInnerX, -rInnerY, offset])
    pointList.append([rOuterX, -rOuterY, offset])
    pointList.append([rOuterX,  rOuterY, offset])
    # pointList.append([rInnerX,  rInnerY, offset])

    # pointList.append([rInnerX, -rInnerY, offset + height])
    pointList.append([rOuterX, -rOuterY, offset + height])
    pointList.append([rOuterX,  rOuterY, offset + height])
    # pointList.append([rInnerX,  rInnerY, offset + height])

    # Write fuel vertices
    numberPoints = len(pointList)
    for i in range(numberPoints):
        newLine = "    (" 
        newLine += str(pointList[i][0]) + " "
        newLine += str(pointList[i][1]) + " " 
        newLine += str(pointList[i][2]) + ")" 
        newLine += "\n"
    
        file.write(newLine)

    file.write("\n")
    
def writeBlock(startIndex, zoneName, meshX, meshZ, file):

    index = startIndex
    newLine = "    hex ( " 

    for i in range(8):
        newLine += str(index) + " "
        index += 1

    newLine += ") " + str(zoneName) + " ("
    newLine +=  str(meshX) + " 1 " + str(meshZ) 
    newLine += ") simpleGrading (1 1 1)"

    file.write(newLine + "\n")

def writePelletBlock(pelletType, startIndex, zoneName, meshDishX, meshLandX, meshChamferX, meshZ, file):

    index = startIndex
    newLine = "    hex ( " 

    # First block - only block for flat pellet
    for i in range(8):
        newLine += str(index) + " "
        index += 1

    newLine += ") " + str(zoneName) + " ("
    if pelletType == 'flat' or pelletType == 'chamfered':
        newLine +=  str(meshLandX) + " 1 " + str(meshZ) 
    else:
        newLine +=  str(meshDishX) + " 1 " + str(meshZ) 
    newLine += ") simpleGrading (1 1 1)"

    file.write(newLine + "\n")

    if not(pelletType == 'flat'):
        # Write land or chamfer piece depending on pellet type
        index = startIndex
        newLine = "    hex ( " 

        newLine += str(index+1) + " "
        newLine += str(index+8) + " "
        newLine += str(index+9) + " "
        newLine += str(index+2) + " "

        newLine += str(index+5) + " "
        newLine += str(index+10) + " "
        newLine += str(index+11) + " "
        newLine += str(index+6) + " "

        newLine += ") " + str(zoneName) + " ("
        if pelletType == 'chamfered':
            newLine +=  str(meshChamferX) + " 1 " + str(meshZ) 
        else:
            newLine +=  str(meshLandX) + " 1 " + str(meshZ) 
        newLine += ") simpleGrading (1 1 1)"

        file.write(newLine + "\n")

    if pelletType == 'dishedChamfered':
        # Write chamfer piece
        index = startIndex
        newLine = "    hex ( " 

        newLine += str(index+8) + " "
        newLine += str(index+12) + " "
        newLine += str(index+13) + " "
        newLine += str(index+9) + " "

        newLine += str(index+10) + " "
        newLine += str(index+14) + " "
        newLine += str(index+15) + " "
        newLine += str(index+11) + " "

        newLine += ") " + str(zoneName) + " ("
        newLine +=  str(meshChamferX) + " 1 " + str(meshZ) 
        newLine += ") simpleGrading (1 1 1)"

        file.write(newLine + "\n")
        
def writeCapBlock(startIndex, zoneName, meshX, meshXInnerPart, meshZ, file):

    index = startIndex
    newLine = "    hex ( " 

    for i in range(8):
        newLine += str(index) + " "
        index += 1

    newLine += ") " + str(zoneName) + " ("
    newLine +=  str(meshXInnerPart) + " 1 " + str(meshZ) 
    newLine += ") simpleGrading (1 1 1)"

    file.write(newLine + "\n")

    # Write second piece
    index = startIndex
    newLine = "    hex ( " 

    newLine += str(index+1) + " "
    newLine += str(index+8) + " "
    newLine += str(index+9) + " "
    newLine += str(index+2) + " "

    newLine += str(index+5) + " "
    newLine += str(index+10) + " "
    newLine += str(index+11) + " "
    newLine += str(index+6) + " "

    newLine += ") " + str(zoneName) + " ("
    newLine +=  str(meshX) + " 1 " + str(meshZ) 
    newLine += ") simpleGrading (1 1 1)"

    file.write(newLine + "\n")

def writeFuelEdge(startIndex, wedgeAngle, rArcPoint, hArcPoint, height, offset, file):
    
    #### Prepare input data ####
    wedgeAngleDegree = float(wedgeAngle)
    wedgeAngleRadiant = wedgeAngleDegree/180*math.pi
    rArcPoint = float(rArcPoint)
    hArcPoint = float(hArcPoint)
    height = float(height)
    offset = float(offset)
    index = int(startIndex)

    if(hArcPoint > 1e-9):

        # Correction for wedge volume
        # Relevant only for large wedge angles
        sinAngle = math.sin(wedgeAngleRadiant/2)
        cosAngle = math.cos(wedgeAngleRadiant/2)
        correction = math.sqrt(sinAngle/(wedgeAngleRadiant/2))

        # Coordinates for edge arc point
        rArcPointX = rArcPoint*cosAngle*correction
        rArcPointY = rArcPoint*sinAngle*correction

        # Write fuel edges
        newLine = "    arc " + str(index + 0) + " " + str(index + 1)
        newLine += "(" 
        newLine += str(rArcPointX) + " "
        newLine += str(-rArcPointY) + " " 
        newLine += str(hArcPoint + offset) + ")" 
        newLine += "\n"

        file.write(newLine)

        newLine = "    arc " + str(index + 3) + " " + str(index + 2)
        newLine += "(" 
        newLine += str(rArcPointX) + " "
        newLine += str(rArcPointY) + " " 
        newLine += str(hArcPoint + offset) + ")" 
        newLine += "\n"

        file.write(newLine)

        newLine = "    arc " + str(index + 4) + " " + str(index + 5)
        newLine += "(" 
        newLine += str(rArcPointX) + " "
        newLine += str(-rArcPointY) + " " 
        newLine += str(height - hArcPoint + offset) + ")" 
        newLine += "\n"

        file.write(newLine)

        newLine = "    arc " + str(index + 7) + " " + str(index + 6)
        newLine += "(" 
        newLine += str(rArcPointX) + " "
        newLine += str(rArcPointY) + " " 
        newLine += str(height - hArcPoint + offset) + ")" 
        newLine += "\n"

        file.write(newLine)

        file.write("\n")
       
def writePatch(patchName, patchFaces, file):

    file.write("    " + str(patchName) + "\n")
    file.write("    {" + "\n")
    file.write("        type " + patchFaces[0][0] + ";\n")
    
    if(patchFaces[0][0] == "regionCoupledOFFBEAT" ):
        file.write("        neighbourPatch " + patchFaces[0][1] + ";\n")
        file.write("        neighbourRegion region0;\n")
        file.write("        owner " + patchFaces[0][2] +";\n")
        if(patchName == "cladInner" or patchName == "fuelOuter"):
            file.write("        updateAMI true;\n")
        else:
            file.write("        updateAMI false;\n")

    file.write("        faces " + "\n")
    file.write("        (" + "\n")

    # Write face (there might be more than one per patch)
    for face in patchFaces:

        v1 = face[3]
        v2 = face[4]
        v3 = face[5]
        v4 = face[6]

        file.write("            (" + str(v1) + " "
                                    + str(v2) + " "
                                    + str(v3) + " "
                                    + str(v4)  + ")" + "\n")

    file.write("        );" + "\n")
    file.write("    }" + "\n\n")

################################
############# MAIN #############
################################

rodDict = sys.argv[1] if len(sys.argv) == 2 else 'rodDict'
  
# Reading the data from the rodDict file 
with open(rodDict) as f:
    print(f"Open {rodDict} dict")
    data = f.read() 
      
# Reconstructing the data as a dictionary 
rodDict = ast.literal_eval(data) 

# Read geometryType, wedgeAngle and convertToMeters
geometryType    = rodDict['geometryType']
wedgeAngle      = rodDict['wedgeAngle']
convertToMeters = rodDict['convertToMeters']

# Read block quantities and names
nBlocksFuel      = rodDict['nBlocksFuel']
nBlocksClad      = rodDict['nBlocksClad']

blockNameFuel = rodDict['blockNameFuel']
blockNameClad = rodDict['blockNameClad']

# Read geometrical quantities
rInnerFuel    = rodDict['rInnerFuel']
rInnerClad    = rodDict['rInnerClad']

rDishFuel     = [0.0 for rIn in rInnerClad]
chamferWidth  = [0.0 for rIn in rInnerClad]
rLandFuel     = [0.0 for rIn in rInnerClad]
rCurvatureDish    = [1e6 for rIn in rInnerClad]
heightChamferFuel = [0.0 for rIn in rInnerClad]
nCellsRDish = [0 for rIn in rInnerClad]
nCellsRLand = [0 for rIn in rInnerClad]
nCellsRChamfer = [0 for rIn in rInnerClad]

rOuterFuel    = rodDict['rOuterFuel']
rOuterClad    = rodDict['rOuterClad']

heightFuel        = rodDict['heightFuel']
heightClad        = rodDict['heightClad']

offsetFuel    = rodDict['offsetFuel']
offsetClad    = rodDict['offsetClad']

mergeCladPatchPairs = rodDict.get('mergeCladPatchPairs', True)
mergeFuelPatchPairs = rodDict.get('mergeFuelPatchPairs', True)

# Read number of sections per block depending on type of geometry
nSectionsFuel = []
if geometryType == '2Ddiscrete':
    nSectionsFuel = rodDict['nPelletsFuel']
elif geometryType == '2Dsmeared' or geometryType == '1D':
    nSectionsFuel = [1 for rIn in rInnerFuel]
else:
    print('ERROR: unknown geometryType. Use either 1D, 2Dsmeared or 2Ddiscrete')
    exit()

nSectionsClad = []
nSectionsClad = [1 for rIn in rInnerClad]

# Read number of radial and axial cells
nCellsRFuel = rodDict['nCellsRFuel']
nCellsRClad = rodDict['nCellsRClad']

nCellsZFuel = []
nCellsZFuel = rodDict['nCellsZFuel']

nCellsZClad = []
nCellsZClad = rodDict['nCellsZClad']

# Adjust dish, land and chamfer settings according to geometry type
if(geometryType == '2Ddiscrete'):
    rDishFuel     = rodDict['rDishFuel']
    chamferWidth = rodDict['chamferWidth']

# Automatically detect pellet type 'dishChamf', 'dish', 'chamf'
pelletType = []
if geometryType == '2Ddiscrete':
    for i in range(len(rDishFuel)):
        rDish = rDishFuel[i]
        wCham = chamferWidth[i]

        if rDish > 0.0 and wCham > 0.0:
            pelletType.append('dishedChamfered')
            rCurvatureDish[i] = rodDict['rCurvatureDish'][i]
            rLandFuel[i]     = rOuterFuel[i] - wCham
            heightChamferFuel[i] = rodDict['heightChamferFuel'][i]
            nCellsRDish[i] = rodDict['nCellsRDish'][i]
            nCellsRLand[i] = rodDict['nCellsRLand'][i]
            nCellsRChamfer[i] = rodDict['nCellsRChamfer'][i]
        elif rDish > 0.0:
            pelletType.append('dished')
            rCurvatureDish[i]    = rodDict['rCurvatureDish'][i]
            nCellsRDish[i] = rodDict['nCellsRDish'][i]
            nCellsRLand[i] = rodDict['nCellsRLand'][i]
        elif wCham > 0.0:
            pelletType.append('chamfered')
            rLandFuel[i]     = rOuterFuel[i] - wCham
            heightChamferFuel[i] = rodDict['heightChamferFuel'][i]
            nCellsRLand[i] = rodDict['nCellsRLand'][i]
            nCellsRChamfer[i] = rodDict['nCellsRChamfer'][i]
        else:
            pelletType.append('flat')
            nCellsRLand[i] = rodDict['nCellsRLand'][i]

# else:
#     rDishFuel     = [rIn + (rOut - rIn)/3.0 for rIn,rOut in zip(rInnerFuel, rOuterFuel)]
#     rLandFuel     = [rIn + (rOut - rIn)*2/3.0 for rIn,rOut in zip(rInnerFuel, rOuterFuel)]
#     rCurvatureDish    = [1e10 for hfuel in heightFuel]
#     heightChamferFuel = [hfuel*0.0 for hfuel in heightFuel]
#     nCellsRDish = [int(round(cellR/3.0)) for cellR in nCellsRFuel]
#     nCellsRLand = [int(round(cellR/3.0)) for cellR in nCellsRFuel]
#     nCellsRChamfer = [cellR - 2*int(round(cellR/3.0)) for cellR in nCellsRFuel]

# Adjust if bottom and top caps are present
bottomCap = False
bottomCapHeight = rodDict['bottomCapHeight']
if(bottomCapHeight > 0):
    bottomCap = True
    nCellsRBottomCap = rodDict['nCellsRBottomCap']
    nCellsZBottomCap = rodDict['nCellsZBottomCap']

topCap = False
topCapHeight = rodDict['topCapHeight']
if(topCapHeight > 0):
    topCap = True
    nCellsRTopCap = rodDict['nCellsRTopCap']
    nCellsZTopCap = rodDict['nCellsZTopCap']

if(bottomCap):
    nBlocksClad += 1
    blockNameClad.insert(0, 'cladding')
    rInnerClad.insert(0, 0.0)
    rOuterClad.insert(0, rOuterClad[0])
    heightClad.insert(0, bottomCapHeight)
    offsetClad -= float(heightClad[0])
    nSectionsClad.insert(0, 1)
    nCellsRClad.insert(0, nCellsRBottomCap)
    nCellsZClad.insert(0, nCellsZBottomCap)

if(topCap):
    nBlocksClad += 1
    blockNameClad.append('cladding')
    rInnerClad.append(0.0)
    rOuterClad.append(rOuterClad[nBlocksClad-2])
    heightClad.append(topCapHeight)
    nSectionsClad.append(1)
    nCellsRClad.append(nCellsRTopCap)
    nCellsZClad.append(nCellsZTopCap)

# Write file header and first lines 
file = open("blockMeshDict", "w+")
writeHeader(file)
file.write("\nconvertToMeters " + str(convertToMeters) + "; \n\n")
file.write("vertices \n(\n")

# Write fuel vertices
offset  = offsetFuel
for i in range(nBlocksFuel):

    n = nSectionsFuel[i]

    rIn  = rInnerFuel[i]
    rOut = rOuterFuel[i]
    h    = heightFuel[i]/(n)

    if(geometryType == '2Ddiscrete'):
        rDish = rDishFuel[i]
        # hDish = heightDishFuel[i]

        rLand = rLandFuel[i]
        hChamfer = heightChamferFuel[i]
        rCurvature = rCurvatureDish[i]

        # Calculate height of dish based on curvature radius
        angleDish = math.acos(rDish/rCurvature)
        hDish = rCurvature*(1 - math.sin(angleDish))

        pType = pelletType[i]

        for sectionI in range(n):
            writePelletVertices(pType, wedgeAngle, rIn, rDish, hDish, rLand, hChamfer, rOut, h, offset, file)
            offset = offset + h
    else:
        for sectionI in range(n):
            writeVertices(wedgeAngle, rIn, rOut, h, offset, file) 
            offset = offset + h     

# Write clad vertices
offset  = offsetClad
for i in range(nBlocksClad):

    if(i==0 and bottomCap):

        rIn  = rInnerClad[1]
        rOut = rOuterClad[i]
        h    = heightClad[i]

        writeCapVertices(wedgeAngle, rIn, rOut, h, offset, file)
        offset = offset + h

    elif(i==nBlocksClad-1 and topCap):

        rIn  = rInnerClad[nBlocksClad-2]
        rOut = rOuterClad[i]
        h    = heightClad[i]

        writeCapVertices(wedgeAngle, rIn, rOut, h, offset, file)
        offset = offset + h

    else:

        n = nSectionsClad[i]

        rIn  = rInnerClad[i]
        rOut = rOuterClad[i]
        h    = heightClad[i]/(n)

        for sectionI in range(n):
            writeVertices(wedgeAngle, rIn, rOut, h, offset, file)
            offset = offset + h

file.write(");\n\n")

# Write fuel blocks
startIndex = 0
file.write("blocks \n(\n")

for i in range(nBlocksFuel):

    meshR = nCellsRFuel[i]
    meshZ = nCellsZFuel[i] 

    name = blockNameFuel[i]
    n = nSectionsFuel[i]

    if(geometryType == '2Ddiscrete'):
        meshDishR = nCellsRDish[i] 
        meshLandR = nCellsRLand[i] 
        meshChamferR = nCellsRChamfer[i] 
        pType = pelletType[i]
        for sectionI in range(n):
            writePelletBlock(pType, startIndex, name, meshDishR, meshLandR, meshChamferR, meshZ, file)
            if pType == 'dishedChamfered':
                startIndex += 16
            elif pType == 'flat':
                startIndex += 8
            else:
                startIndex += 12

    else:
        for sectionI in range(n):
            writeBlock(startIndex, name, meshR, meshZ, file)
            startIndex += 8

# Write clad blocks
for i in range(nBlocksClad):

    if(i==0 and bottomCap):

        meshR = nCellsRClad[1] 
        meshRInnerPart = nCellsRClad[i] 
        meshZ = nCellsZClad[i] 

        name = blockNameClad[i]

        writeCapBlock(startIndex, name, meshR, meshRInnerPart, meshZ, file)
        startIndex += 12

    elif(i==nBlocksClad-1 and topCap):

        meshR = nCellsRClad[nBlocksClad-2] 
        meshRInnerPart = nCellsRClad[i] 
        meshZ = nCellsZClad[i] 

        name = blockNameClad[i]

        writeCapBlock(startIndex, name, meshR, meshRInnerPart, meshZ, file)
        startIndex += 12

    else:
        meshR = nCellsRClad[i] 
        meshZ = nCellsZClad[i] 

        name = blockNameClad[i]
        n = nSectionsClad[i]

        for sectionI in range(n):
            writeBlock(startIndex, name, meshR, meshZ, file)
            startIndex += 8

file.write(");\n\n")

# Write fuel edges
startIndex = 0
file.write("edges \n(\n")

offset = offsetFuel
if(geometryType == '2Ddiscrete'):

    for i in range(nBlocksFuel):
        pType = pelletType[i]

        if pType == 'dished' or pType == 'dishedChamfered':
            n = nSectionsFuel[i]

            rDish = rDishFuel[i]
            h     = heightFuel[i]/(n)

            rCurvature = rCurvatureDish[i]

            # Calculate height of dish based on curvature radius
            angleDish = math.acos(rDish/rCurvature)

            # Calculate height of arc point based on curvature radius
            angleArcPoint = math.acos(rDish/rCurvature/2.0)
            hArcPoint = rCurvature*(math.sin(angleArcPoint) - math.sin(angleDish))
            rArcPoint = rDish/2.0

            for sectionI in range(n):
                writeFuelEdge(startIndex, wedgeAngle, rArcPoint, hArcPoint, h, offset, file)
                offset = offset + h

                if pType == 'dishedChamfered':
                    startIndex += 16

                if pType == 'dished':
                    startIndex += 12

file.write(");\n\n")

# Collect boundaries into patchDict (defaultDict allows to append to entries)
startIndex = 0

patchDict = defaultdict(list)

# Prepare mergePatchDict (one slave patch for each master patch)
mergePatchDict = defaultdict(list)

globalFuelSectionIndex = 1
globalCladSectionIndex = 1

# Collect fuel patches
for iBlock in range(nBlocksFuel):

    n = nSectionsFuel[iBlock]

    for sectionI in range(n):

        i = startIndex

        # Collect bottom boundary
        bottomName  = 'fuelBottom'
        bottomType  = 'patch'
        bottomNeig  = 'none'
        bottomOwner = 'false'

        if globalFuelSectionIndex > 1:
            bottomName += '_' + str(globalFuelSectionIndex)

            if(not(mergeFuelPatchPairs)):
                bottomType  = 'regionCoupledOFFBEAT'
                bottomNeig  = 'fuelTop_' + str(globalFuelSectionIndex-1)

        elif(bottomCap):
            bottomType  = 'regionCoupledOFFBEAT'
            bottomNeig  = 'bottomCapInner'
            bottomOwner = 'true'

        if geometryType == '1D': bottomType = 'empty'

        # Collect bottom dish/chamfer boundaries
        dishChamfBottomName  = 'fuelBottom'
        dishChamfBottomType  = 'patch'
        dishChamfBottomNeig  = 'none'
        dishChamfBottomOnwer = 'false'

        if globalFuelSectionIndex > 1:
            dishChamfBottomName = 'dishChamferBottom_' + str(globalFuelSectionIndex)

            # if(not(mergeFuelPatchPairs)):
            #     dishChamfBottomType  = 'regionCoupledOFFBEAT'
            #     dishChamfBottomNeig  = 'dishChamferTop_' + str(globalFuelSectionIndex-1)

        elif(bottomCap):
            dishChamfBottomType  = 'regionCoupledOFFBEAT'
            dishChamfBottomNeig  = 'bottomCapInner'
            dishChamfBottomOnwer = 'true'

        # Collect top boundary
        topName  = 'fuelTop'
        topType  = 'patch'
        topNeig  = 'none'
        topOwner = 'true'

        if globalFuelSectionIndex < sum(nSectionsFuel):
            topName += '_' + str(globalFuelSectionIndex)

            if(not(mergeFuelPatchPairs)):
                topType  = 'regionCoupledOFFBEAT'
                topNeig  = 'fuelBottom_' + str(globalFuelSectionIndex+1)
            else:
                mergePatchDict[topName] = 'fuelBottom_' + str(globalFuelSectionIndex+1)

        if geometryType == '1D': topType = 'empty'

        # Collect top dish/chamfer boundaries
        dishChamfTopName  = 'fuelTop'
        dishChamfTopType  = 'patch'
        dishChamfTopNeig  = 'none'
        dishChamfTopOwner = 'false'

        if globalFuelSectionIndex < sum(nSectionsFuel):
            dishChamfTopName = 'dishChamferTop_' + str(globalFuelSectionIndex)

            # if(not(mergeFuelPatchPairs)):
            #     dishChamfTopType  = 'regionCoupledOFFBEAT'
            #     dishChamfTopNeig  = 'dishChamferBottom_' + str(globalFuelSectionIndex-1)

        globalFuelSectionIndex += 1

        # Collect inner boundary
        innerName  = 'fuelInner'
        innerType  = 'patch' 
        innerNeig  = 'none'
        innerOwner = 'false'

        # Collect outer boundary
        outerName  = 'fuelOuter'
        outerType  = 'regionCoupledOFFBEAT' 
        outerNeig  = 'cladInner' 
        outerOwner = 'true'

        # Collect wedge boundaries
        wedgeFrontName = 'fuelFront'
        wedgeBackName  = 'fuelBack'

        # Bottom
        if(geometryType == '2Ddiscrete'):

            pType = pelletType[iBlock]

            if pType == 'flat':
                patchDict[str(bottomName)].append(
                    [bottomType, bottomNeig, bottomOwner, i+0, i+3, i+2, i+1])

            elif pType == 'dished':   
                patchDict[str(dishChamfBottomName)].append(
                    [dishChamfBottomType, dishChamfBottomNeig, dishChamfBottomOnwer, i+0, i+3, i+2, i+1])             
                patchDict[str(bottomName)].append(
                    [bottomType, bottomNeig, bottomOwner, i+1, i+2, i+9, i+8])

            elif pType == 'chamfered':            
                patchDict[str(bottomName)].append(
                    [bottomType, bottomNeig, bottomOwner, i+0, i+3, i+2, i+1])
                patchDict[str(dishChamfBottomName)].append(
                    [dishChamfBottomType, dishChamfBottomNeig, dishChamfBottomOnwer, i+1, i+2, i+9, i+8]) 

            elif pType == 'dishedChamfered':
                patchDict[str(dishChamfBottomName)].append(
                    [dishChamfBottomType, dishChamfBottomNeig, dishChamfBottomOnwer, i+0, i+3, i+2, i+1]) 
                patchDict[str(bottomName)].append(
                    [bottomType, bottomNeig, bottomOwner, i+1, i+2, i+9, i+8])
                patchDict[str(dishChamfBottomName)].append(
                    [dishChamfBottomType, dishChamfBottomNeig, dishChamfBottomOnwer, i+8, i+9, i+13, i+12])
        else:
            patchDict[str(bottomName)].append(
                [bottomType, bottomNeig, bottomOwner, i+0, i+3, i+2, i+1 ])

        # Top
        if(geometryType == '2Ddiscrete'):

            pType = pelletType[iBlock]

            if pType == 'flat':
                patchDict[str(topName)].append(
                    [topType, topNeig, topOwner, i+4, i+5, i+6, i+7])

            elif pType == 'dished':   
                patchDict[str(dishChamfTopName)].append(
                    [dishChamfTopType, dishChamfTopNeig, dishChamfTopOwner, i+4, i+5, i+6, i+7])             
                patchDict[str(topName)].append(
                    [topType, topNeig, topOwner, i+5, i+10, i+11, i+6])

            elif pType == 'chamfered':            
                patchDict[str(topName)].append(
                    [topType, topNeig, topOwner, i+4, i+5, i+6, i+7])
                patchDict[str(dishChamfTopName)].append(
                    [dishChamfTopType, dishChamfTopNeig, dishChamfTopOwner, i+5, i+10, i+11, i+6]) 

            elif pType == 'dishedChamfered':
                patchDict[str(dishChamfTopName)].append(
                    [dishChamfTopType, dishChamfTopNeig, dishChamfTopOwner, i+4, i+5, i+6, i+7]) 
                patchDict[str(topName)].append(
                    [topType, topNeig, topOwner, i+5, i+10, i+11, i+6])
                patchDict[str(dishChamfTopName)].append(
                    [dishChamfTopType, dishChamfTopNeig, dishChamfTopOwner, i+10, i+14, i+15, i+11])
        else:
            patchDict[str(topName)].append(
                [topType, topNeig, topOwner, i+4, i+5, i+6, i+7]) 

        # Inner
        if(rInnerFuel[iBlock] > 0.0):
            patchDict[str(innerName)].append(
                [innerType, innerNeig, innerOwner, i+0, i+4, i+7, i+3])   

        # Outer
        if geometryType == '2Ddiscrete':

            pType = pelletType[iBlock]
            
            if pType == 'dishedChamfered':
                patchDict[str(outerName)].append(
                    [outerType, outerNeig, outerOwner, i+12, i+13, i+15, i+14])    
            
            elif pType == 'dished' or pType == 'chamfered':
                patchDict[str(outerName)].append(
                    [outerType, outerNeig, outerOwner, i+8, i+9, i+11, i+10])    
            
            elif pType == 'flat':
                patchDict[str(outerName)].append(
                    [outerType, outerNeig, outerOwner, i+1, i+2, i+6, i+5])    
        else:
            patchDict[str(outerName)].append(
                [outerType, outerNeig, outerOwner, i+1, i+2, i+6, i+5])             

        # Wedge
        if(geometryType == '2Ddiscrete'):

            pType = pelletType[iBlock]

            patchDict[str(wedgeFrontName)].append(
                ["wedge", "none", "false", i+0, i+1, i+5, i+4]) 

            if pType == 'dished' or pType == 'chamfered' or pType == 'dishedChamfered':
                patchDict[str(wedgeFrontName)].append(
                    ["wedge", "none", "false", i+1, i+8, i+10, i+5]) 

            if pType == 'dishedChamfered':
                patchDict[str(wedgeFrontName)].append(
                    ["wedge", "none", "false", i+8, i+12, i+14, i+10])  
        else:
            patchDict[str(wedgeFrontName)].append(
                ["wedge", "none", "false", i+0, i+1, i+5, i+4]) 

        if(geometryType == '2Ddiscrete'):

            pType = pelletType[iBlock]

            patchDict[str(wedgeBackName)].append(
                ["wedge", "none", "false", i+3, i+7, i+6, i+2])

            if pType == 'dished' or pType == 'chamfered' or pType == 'dishedChamfered':            
                patchDict[str(wedgeBackName)].append(
                    ["wedge", "none", "false", i+2, i+6, i+11, i+9])

            if pType == 'dishedChamfered':
                patchDict[str(wedgeBackName)].append(
                    ["wedge", "none", "false", i+9, i+11, i+15, i+13]) 
        else:
            patchDict[str(wedgeBackName)].append(
                ["wedge", "none", "false", i+3, i+7, i+6, i+2])

        if geometryType == '2Ddiscrete':

            pType = pelletType[iBlock]

            if pType == 'flat':
                startIndex += 8
            elif pType == 'dishedChamfered':
                startIndex += 16
            else:
                startIndex += 12
        else:
            startIndex += 8

# Collect clad patches
for iBlock in range(nBlocksClad):

    n = nSectionsClad[iBlock]

    for sectionI in range(n):

        i = startIndex

        # Collect bottom boundary
        bottomName  = 'cladBottom'
        bottomType  = 'patch'
        bottomNeig  = 'none'
        bottomOwner = 'false'

        if globalCladSectionIndex > 1:
            bottomName += '_' + str(globalCladSectionIndex)

            if(not(mergeCladPatchPairs)):
                bottomType  = 'regionCoupledOFFBEAT'
                bottomNeig  = 'cladTop_' + str(globalCladSectionIndex-1)

        if geometryType == '1D': bottomType = 'empty'

        # Collect top boundary
        topName  = 'cladTop'
        topType  = 'patch'
        topNeig  = 'none'
        topOwner = 'true'

        if globalCladSectionIndex < sum(nSectionsClad):
            topName += '_' + str(globalCladSectionIndex)

            if(not(mergeCladPatchPairs)):
                topType  = 'regionCoupledOFFBEAT'
                topNeig  = 'cladBottom_' + str(globalCladSectionIndex+1)
            else:
                mergePatchDict[topName] = 'cladBottom_' + str(globalCladSectionIndex+1)

        if geometryType == '1D': topType = 'empty'

        globalCladSectionIndex += 1

        # Collect bottom cap inner boundary
        bottomCapInnerName  = 'bottomCapInner'
        bottomCapInnerType  = 'regionCoupledOFFBEAT'
        bottomCapInnerNeig  = 'fuelBottom'
        bottomCapInnerOwner = 'false'

        # Collect top cap inner boundary
        topCapInnerName  = 'topCapInner'
        topCapInnerType  = 'patch'
        topCapInnerNeig  = 'none'
        topCapInnerOwner = 'true'

        # Collect inner boundary
        innerName  = 'cladInner'
        innerType  = 'regionCoupledOFFBEAT' 
        innerNeig  = 'fuelOuter'
        innerOwner = 'false'

        # Collect outer boundary
        outerName  = 'cladOuter'
        outerType  = 'patch' 
        outerNeig  = 'none' 
        outerOwner = 'false'

        # Collect wedge boundaries
        wedgeFrontName = 'cladFront'
        wedgeBackName  = 'cladBack'

        # Collect bottom
        if(bottomCap and iBlock==0):
            patchDict[str(bottomName)].append(
                [bottomType, bottomNeig, bottomOwner, i+0, i+3, i+2, i+1 ])
            patchDict[str(bottomName)].append(
                [bottomType, bottomNeig, bottomOwner, i+1, i+2, i+9, i+8 ])
        elif(topCap and iBlock==nBlocksClad-1):
            patchDict[str(bottomName)].append(
                [bottomType, bottomNeig, bottomOwner, i+1, i+2, i+9, i+8 ])
        else:
            patchDict[str(bottomName)].append(
                [bottomType, bottomNeig, bottomOwner, i+0, i+3, i+2, i+1 ])

        # Collect top cap inner
        if(topCap and iBlock==nBlocksClad-1):
            patchDict[str(topCapInnerName)].append(
                [topCapInnerType, topCapInnerNeig, topCapInnerOwner, i+0, i+3, i+2, i+1 ])

        # Collect top
        if(topCap and iBlock==nBlocksClad-1):
            patchDict[str(topName)].append(
                [topType, topNeig, topOwner, i+4, i+5, i+6, i+7]) 
            patchDict[str(topName)].append(
                [topType, topNeig, topOwner, i+5, i+10, i+11, i+6]) 
        elif(bottomCap and iBlock==0):
            patchDict[str(topName)].append(
                [topType, topNeig, topOwner, i+5, i+10, i+11, i+6]) 
        else:
            patchDict[str(topName)].append(
                [topType, topNeig, topOwner, i+4, i+5, i+6, i+7]) 

        # Collect top cap inner
        if(bottomCap and iBlock==0):
            patchDict[str(bottomCapInnerName)].append(
                [bottomCapInnerType, bottomCapInnerNeig, bottomCapInnerOwner, i+4, i+5, i+6, i+7 ])

        # Collect inner
        if(not((bottomCap and iBlock==0) or (topCap and iBlock==nBlocksClad-1))):
            patchDict[str(innerName)].append(
                [innerType, innerNeig, innerOwner, i+0, i+4, i+7, i+3])   

        # Collect outer
        if((bottomCap and iBlock==0) or (topCap and iBlock==nBlocksClad-1)):
            patchDict[str(outerName)].append(
                [outerType, outerNeig, outerOwner, i+8, i+9, i+11, i+10])   
        else:
            patchDict[str(outerName)].append(
                [outerType, outerNeig, outerOwner, i+1, i+2, i+6, i+5])    

        # Collect wedgeFront
        if((bottomCap and iBlock==0) or (topCap and iBlock==nBlocksClad-1)):
            patchDict[str(wedgeFrontName)].append(
                ["wedge", "none", "false", i+0, i+1, i+5, i+4]) 
            patchDict[str(wedgeFrontName)].append(
                ["wedge", "none", "false", i+1, i+8, i+10, i+5]) 
        else:
            patchDict[str(wedgeFrontName)].append(
                ["wedge", "none", "false", i+0, i+1, i+5, i+4]) 

        # Collect wedgeBackt
        if((bottomCap and iBlock==0) or (topCap and iBlock==nBlocksClad-1)):
            patchDict[str(wedgeBackName)].append(
                ["wedge", "none", "false", i+3, i+7, i+6, i+2])
            patchDict[str(wedgeBackName)].append(
                ["wedge", "none", "false", i+2, i+6, i+11, i+9])
        else:
            patchDict[str(wedgeBackName)].append(
                ["wedge", "none", "false", i+3, i+7, i+6, i+2])
        
        if((bottomCap and iBlock==0) or (topCap and iBlock==nBlocksClad-1)):
            startIndex += 12
        else:
            startIndex += 8

# Write boundaries
file.write("boundary \n(\n")
for patchName in patchDict:
    patchFaces = patchDict[patchName]
    writePatch(patchName, patchFaces, file)

file.write(");\n\n")

# Write mergePatchPairs
file.write("mergePatchPairs \n(\n")
for masterPatchName in mergePatchDict:
    slavePatchName = mergePatchDict[masterPatchName]
    file.write("\t(")
    file.write(masterPatchName + " " + slavePatchName)
    file.write(")\n")

file.write(");\n\n")
