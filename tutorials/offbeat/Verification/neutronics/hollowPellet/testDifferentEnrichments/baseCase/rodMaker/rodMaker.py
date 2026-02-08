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

# importing the module 
import ast 

# Function to write header of the dictionary:
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
    for ii in range(numberPoints):
        newLine = "    (" 
        newLine += str(pointList[ii][0]) + " "
        newLine += str(pointList[ii][1]) + " "

        newLine += str(pointList[ii][2]) + ") //" + str(i_point*numberPoints+ii)
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

       
def writePatch(patchName, patchFaces, file):

    file.write("    " + str(patchName) + "\n")
    file.write("    {" + "\n")
    file.write("        type " + patchFaces[0][0] + ";\n")
    
    if(patchFaces[0][0] == "regionCoupledOFFBEAT"):
        file.write("        neighbourPatch " + patchFaces[0][1] + ";\n")
        file.write("        neighbourRegion region0;\n")
        file.write("        owner " + patchFaces[0][2] +";\n")
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
  
# Reading the data from the rodDict file 
with open('rodDict') as f: 
    data = f.read() 
      
# Reconstructing the data as a dictionary 
rodDict = ast.literal_eval(data) 

# Read geometryType, wedgeAngle and convertToMeters
geometryType    = rodDict['geometryType']
wedgeAngle      = rodDict['wedgeAngle']
convertToMeters = rodDict['convertToMeters']

# Read block quantities and names
nBlocksFuel = rodDict['nBlocksFuel']
nBlocksClad = rodDict['nBlocksClad']

blockNameFuel = rodDict['blockNameFuel']
blockNameClad = rodDict['blockNameClad']

# Read geometrical quantities
rInnerFuel    = rodDict['rInnerFuel']
rInnerClad    = rodDict['rInnerClad']

rOuterFuel    = rodDict['rOuterFuel']
rOuterClad    = rodDict['rOuterClad']

heightFuel    = rodDict['heightFuel']
heightClad    = rodDict['heightClad']

offsetFuel    = rodDict['offsetFuel']
offsetClad    = rodDict['offsetClad']

# Read number of sections per block depending on type of geometry
nSectionsFuel = []
if geometryType == '1D':
    nSectionsFuel = rodDict['nSlicesFuel']
elif geometryType == '2Dsmeared':
    nSectionsFuel = [1 for rIn in rInnerFuel]
elif geometryType == '2Ddiscrete':
    nSectionsFuel = rodDict['nPelletsFuel']
else:
    print('ERROR: unknown geometryType. Use either 1D, 2Dsmeared or 2Ddiscrete')
    exit()

nSectionsClad = []
if geometryType == '1D':
    nSectionsClad = rodDict['nSlicesClad']
else:
    nSectionsClad = [1 for rIn in rInnerClad]

# Read number of radial and axial cells
nCellsRFuel = rodDict['nCellsRFuel']
nCellsRClad = rodDict['nCellsRClad']

nCellsZFuel = []
if geometryType == '1D':
    nCellsZFuel = [1 for rIn in rInnerFuel]
else:
    nCellsZFuel = rodDict['nCellsZFuel']

nCellsZClad = []
if geometryType == '1D':
    nCellsZClad = [1 for rIn in rInnerClad]
else:
    nCellsZClad = rodDict['nCellsZClad']

# Write file header and first lines 
file = open("blockMeshDict", "w+")
writeHeader(file)
file.write("\nconvertToMeters " + str(convertToMeters) + "; \n\n")
file.write("vertices \n(\n")

i_point = 0
# Write fuel vertices
offset  = offsetFuel
for i in range(len(rInnerFuel)):

    n = nSectionsFuel[i]

    rIn  = rInnerFuel[i]
    rOut = rOuterFuel[i]
    h    = heightFuel[i]/(n)

    for sectionI in range(n):
        i_point = i
        writeVertices(wedgeAngle, rIn, rOut, h, offset, file)
        offset = offset + h

# Write clad vertices
offset  = offsetClad
for i in range(len(rInnerClad)):

    n = nSectionsClad[i]

    rIn  = rInnerClad[i]
    rOut = rOuterClad[i]
    h    = heightClad[i]/(n)

    for sectionI in range(n):
        i_point = i + len(rInnerFuel)
        writeVertices(wedgeAngle, rIn, rOut, h, offset, file)
        offset = offset + h

file.write(");\n\n")

# Write fuel blocks
startIndex = 0
file.write("blocks \n(\n")

for i in range(len(rInnerFuel)):
    meshR = nCellsRFuel[i] 
    meshZ = nCellsZFuel[i] 

    name = blockNameFuel[i]
    n = nSectionsFuel[i]

    for sectionI in range(n):
        writeBlock(startIndex, name, meshR, meshZ, file)
        startIndex += 8

# Write clad blocks
for i in range(len(rInnerClad)):
    meshR = nCellsRClad[i] 
    meshZ = nCellsZClad[i] 

    name = blockNameClad[i]
    n = nSectionsClad[i]

    for sectionI in range(n):
        writeBlock(startIndex, name, meshR, meshZ, file)
        startIndex += 8

file.write(");\n\n")

# Collect boundaries into patchDict (defaultDict allows to append to entries)
startIndex = 0

patchDict = defaultdict(list)

globalFuelSectionIndex = 1
globalCladSectionIndex = 1

# Collect fuel patches
for iBlock in range(nBlocksFuel):

    n = nSectionsFuel[iBlock]

    for sectionI in range(n):

        i = startIndex

        # Collect bottom boundary
        bottomName  = 'bottom'
        bottomType  = 'patch'
        bottomNeig  = 'none'
        bottomOnwer = 'false'

        if globalFuelSectionIndex > 1:
            bottomName += '_' + str(globalFuelSectionIndex)
            bottomType  = 'regionCoupledOFFBEAT'
            bottomNeig  = 'top_' + str(globalFuelSectionIndex-1)

        if geometryType == '1D': bottomType = 'empty'

        # Collect top boundary
        topName  = 'top'
        topType  = 'patch'
        topNeig  = 'none'
        topOwner = 'true'

        if globalFuelSectionIndex < sum(nSectionsFuel):
            topName += '_' + str(globalFuelSectionIndex)
            topType  = 'regionCoupledOFFBEAT'
            topNeig  = 'bottom_' + str(globalFuelSectionIndex+1)

        if geometryType == '1D': topType = 'empty'

        globalFuelSectionIndex += 1

        # Collect inner boundary
        innerName  = 'fuelInnerSurface'
        innerType  = 'patch' 
        innerNeig  = 'none'
        innerOwner = 'false'

        # Collect outer boundary
        outerName  = 'fuelOuterSurface'
        outerType  = 'regionCoupledOFFBEAT' 
        outerNeig  = 'claddingInnerSurface' 
        outerOwner = 'false'

        # Collect wedge boundaries
        wedgeFrontName = 'wedgeFuelFront'
        wedgeBackName  = 'wedgeFuelBack'

        patchDict[str(bottomName)].append(
            [bottomType, bottomNeig, bottomOnwer, i+0, i+3, i+2, i+1])

        patchDict[str(topName)].append(
            [topType, topNeig, topOwner, i+4, i+5, i+6, i+7]) 

        patchDict[str(innerName)].append(
            [innerType, innerNeig, innerOwner, i+0, i+4, i+7, i+3])   

        patchDict[str(outerName)].append(
            [outerType, outerNeig, outerOwner, i+1, i+2, i+6, i+5])    

        patchDict[str(wedgeFrontName)].append(
            ["wedge", "none", "false", i+0, i+1, i+5, i+4]) 
        patchDict[str(wedgeBackName)].append(
            ["wedge", "none", "false", i+3, i+7, i+6, i+2])
        
        startIndex += 8

# Collect clad patches
for iBlock in range(nBlocksClad):

    n = nSectionsClad[iBlock]

    for sectionI in range(n):

        i = startIndex

        # Collect bottom boundary
        bottomName  = 'bottomCladding'
        bottomType  = 'patch'
        bottomNeig  = 'none'
        bottomOnwer = 'false'

        if globalCladSectionIndex > 1:
            bottomName += '_' + str(globalCladSectionIndex)
            bottomType  = 'regionCoupledOFFBEAT'
            bottomNeig  = 'topCladding' + str(globalCladSectionIndex-1)

        if geometryType == '1D': bottomType = 'empty'

        # Collect top boundary
        topName  = 'topCladding'
        topType  = 'patch'
        topNeig  = 'none'
        topOwner = 'true'

        if globalCladSectionIndex < sum(nSectionsClad):
            topName += '_' + str(globalCladSectionIndex)
            topType  = 'regionCoupledOFFBEAT'
            topNeig  = 'bottomCladding' + str(globalCladSectionIndex+1)

        if geometryType == '1D': topType = 'empty'

        globalCladSectionIndex += 1

        # Collect inner boundary
        innerName  = 'claddingInnerSurface'
        innerType  = 'regionCoupledOFFBEAT' 
        innerNeig  = 'fuelOuterSurface'
        innerOwner = 'true'

        # Collect outer boundary
        outerName  = 'claddingOuterSurface'
        outerType  = 'patch' 
        outerNeig  = 'none' 
        outerOwner = 'false'

        # Collect wedge boundaries
        wedgeFrontName = 'wedgeCladdingFront'
        wedgeBackName  = 'wedgeCladdingBack'

        patchDict[str(bottomName)].append(
            [bottomType, bottomNeig, bottomOnwer, i+0, i+3, i+2, i+1])

        patchDict[str(topName)].append(
            [topType, topNeig, topOwner, i+4, i+5, i+6, i+7]) 

        patchDict[str(innerName)].append(
            [innerType, innerNeig, innerOwner, i+0, i+4, i+7, i+3])   

        patchDict[str(outerName)].append(
            [outerType, outerNeig, outerOwner, i+1, i+2, i+6, i+5])    

        patchDict[str(wedgeFrontName)].append(
            ["wedge", "none", "false", i+0, i+1, i+5, i+4]) 
        patchDict[str(wedgeBackName)].append(
            ["wedge", "none", "false", i+3, i+7, i+6, i+2])
        
        startIndex += 8

# Write boundaries
file.write("boundary \n(\n")
for patchName in patchDict:
    patchFaces = patchDict[patchName]
    writePatch(patchName, patchFaces, file)

file.write(");\n\n")
