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

# For Kernel and IPyC
def writeVertices(wedgeAngle, rInner, thick, file):

    #### Prepare input data ####
    wedgeAngleDegree = float(wedgeAngle)
    wedgeAngleRadiant = wedgeAngleDegree/180*math.pi
    rInner = float(rInner)
    thick = float(thick)
    rOuter = rInner + thick

    # Prepare list of points
    pointList = []

    sinAngle = math.sin(wedgeAngleRadiant/2)
    cosAngle = math.cos(wedgeAngleRadiant/2)
    #correction = math.sqrt(sinAngle/(wedgeAngleRadiant/2))

    # Coordinates for vertices
    rInnerX = rInner * cosAngle / math.sqrt(1 + sinAngle*sinAngle)
    rInnerY = rInner * sinAngle / math.sqrt(1 + sinAngle*sinAngle)
    rInnerZ = rInner * sinAngle / math.sqrt(1 + sinAngle*sinAngle)

    rOuterX = rOuter * cosAngle / math.sqrt(1 + sinAngle*sinAngle)
    rOuterY = rOuter * sinAngle / math.sqrt(1 + sinAngle*sinAngle)
    rOuterZ = rOuter * sinAngle / math.sqrt(1 + sinAngle*sinAngle)

    # Define block vertices for each block
    pointList.append([rInnerX,  rInnerY, -rInnerZ])
    pointList.append([rInnerX,  rInnerY,  rInnerZ])
    pointList.append([rInnerX, -rInnerY,  rInnerZ])
    pointList.append([rInnerX, -rInnerY, -rInnerZ])
    pointList.append([rOuterX,  rOuterY, -rOuterZ])
    pointList.append([rOuterX,  rOuterY,  rOuterZ])
    pointList.append([rOuterX, -rOuterY,  rOuterZ])
    pointList.append([rOuterX, -rOuterY, -rOuterZ])


    # Write vertices
    numberPoints = len(pointList)
    for i in range(numberPoints):
        newLine = "    ("
        newLine += str(pointList[i][0]) + " "
        newLine += str(pointList[i][1]) + " "
        newLine += str(pointList[i][2]) + ")"
        newLine += "\n"

        file.write(newLine)

    file.write("\n")

# For Buffer, SiC and OPyC
def writeHalfVertices(wedgeAngle, rInner, thick, file):

    #### Prepare input data ####
    wedgeAngleDegree = float(wedgeAngle)
    wedgeAngleRadiant = wedgeAngleDegree/180*math.pi
    rInner = float(rInner)
    thick = float(thick)

    # Prepare list of points
    pointList = []

    sinAngle = math.sin(wedgeAngleRadiant/2)
    cosAngle = math.cos(wedgeAngleRadiant/2)
    #correction = math.sqrt(sinAngle/(wedgeAngleRadiant/2))

    # Coordinates for vertices
    rX = (rInner + thick) * cosAngle / math.sqrt(1 + sinAngle*sinAngle)
    rY = (rInner + thick) * sinAngle / math.sqrt(1 + sinAngle*sinAngle)
    rZ = (rInner + thick) * sinAngle / math.sqrt(1 + sinAngle*sinAngle)

    # Define block vertices for each block
    pointList.append([rX,  rY, -rZ])
    pointList.append([rX,  rY,  rZ])
    pointList.append([rX, -rY,  rZ])
    pointList.append([rX, -rY, -rZ])


    # Write vertices
    numberPoints = len(pointList)
    for i in range(numberPoints):
        newLine = "    ("
        newLine += str(pointList[i][0]) + " "
        newLine += str(pointList[i][1]) + " "
        newLine += str(pointList[i][2]) + ")"
        newLine += "\n"

        file.write(newLine)

    file.write("\n")



def writeBlock(startIndex, zoneName, meshX, file):

    index = startIndex
    newLine = "    hex ( "

    for i in range(8):
        newLine += str(index) + " "
        index += 1

    newLine += ") " + str(zoneName) + " ("
    newLine +=  "1 1 " + str(meshX)
    newLine += ") simpleGrading (1 1 1)"

    file.write(newLine + "\n")


def writePatch(patchName, patchFaces, file):

    file.write("    " + str(patchName) + "\n")
    file.write("    {" + "\n")
    file.write("        type " + patchFaces[0][0] + ";\n")

    if(patchFaces[0][0] == "regionCoupledOFFBEAT" ):
        file.write("        neighbourPatch " + patchFaces[0][1] + ";\n")
        file.write("        neighbourRegion region0;\n")
        file.write("        owner " + patchFaces[0][2] +";\n")
        file.write("        updateAMI true;\n")

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

#********************************************************************************
#***************************************** MAIN *********************************
#********************************************************************************

# Reading the data from the rodDict file
with open('sphereDict') as f:
    data = f.read()

# Reconstructing the data as a dictionary
sphereDict = ast.literal_eval(data)

# Read wedgeAngle and convertToMeters
wedgeAngle      = sphereDict['wedgeAngle']
convertToMeters = sphereDict['convertToMeters']
rKernel         = sphereDict['rKernel']
tBuffer         = sphereDict['tBuffer']
tIPyC           = sphereDict['tIPyC']
tSiC            = sphereDict['tSiC']
tOPyC           = sphereDict['tOPyC']
meshKernel      = sphereDict['meshKernel']
meshBuffer      = sphereDict['meshBuffer']
meshIPyC        = sphereDict['meshIPyC']
meshSiC         = sphereDict['meshSiC']
meshOPyC        = sphereDict['meshOPyC']
Kernel          = sphereDict['Kernel']
Buffer          = sphereDict['Buffer']
IPyC            = sphereDict['IPyC']
SiC             = sphereDict['SiC']
OPyC            = sphereDict['OPyC']
gap             = sphereDict['gap']
detached        = sphereDict['detached']

partList = []
partList.append([Kernel, Buffer, IPyC, SiC, OPyC])
rList = []
rList.append([rKernel, tBuffer, tIPyC, tSiC, tOPyC])
#*******************************Write******************************************

#----------------------------------vertice----------------------------------------

# Write file header and first lines
file = open("blockMeshDict", "w+")
writeHeader(file)
file.write("\nconvertToMeters " + str(convertToMeters) + "; \n\n")
file.write("vertices \n(\n")

first_part = True
rIn = 0

for i in range(len(partList[0])):
    if partList[0][i] == True:
        if first_part == True:
            writeVertices(wedgeAngle, rIn, rList[0][i], file)
            first_part = False
        else:
            writeHalfVertices(wedgeAngle, rIn, rList[0][i], file)
        if i == 1 and detached:
            first_part = True
            rIn += gap
    rIn += rList[0][i]

file.write(");\n\n")


#----------------------------------block----------------------------------------

file.write("blocks \n(\n")

i_block = 0
# Kernel block
if Kernel == True:
    writeBlock(i_block, 'Kernel', meshKernel, file)
    i_block += 4

# Buffer block
if Buffer == True:
    writeBlock(i_block, 'Buffer', meshBuffer, file)
    if detached:
        i_block += 8
    else:
        i_block += 4


# IPyC block
if IPyC == True:
    writeBlock(i_block, 'IPyC', meshIPyC, file)
    i_block += 4

# SiC block
if SiC == True:
    writeBlock(i_block, 'SiC', meshSiC, file)
    i_block += 4

# OpyC block
if OPyC == True:
    writeBlock(i_block, 'OPyC', meshOPyC, file)

file.write(");\n\n")

#----------------------------------edge----------------------------------------
file.write("edges \n(\n")
file.write(");\n\n")

#--------------------------------boundary--------------------------------------

# Collect boundaries into patchDict (defaultdict allows to append to entries)

patchDict = defaultdict(list)

# For the wedge boundaries
i_patch = 0
for j in range(5):
    if partList[0][j] == True:
        patchDict[str('front')].append(
        ["wedge", "none", "false", i_patch+1, i_patch+2, i_patch+6, i_patch+5])
        patchDict[str('back')].append(
        ["wedge", "none", "false", i_patch+0, i_patch+4, i_patch+7, i_patch+3])
        patchDict[str('top')].append(
        ["wedge", "none", "false", i_patch+0, i_patch+1, i_patch+5, i_patch+4])
        patchDict[str('bottom')].append(
        ["wedge", "none", "false", i_patch+2, i_patch+3, i_patch+7, i_patch+6])
        i_patch += 4
        if j == 1 and detached:
            i_patch += 4


# For the outer boundaries
patchDict[str('outer')].append(
["patch", "none", "false", i_patch, i_patch+1, i_patch+2, i_patch+3])

# Buffer outer - IPyC inner boundary
if Buffer == True and detached:
    if Kernel == True:
        patchDict[str('bufferOuter')].append(
        ["regionCoupledOFFBEAT", "ipycInner", "true", 8, 9, 10, 11])
        patchDict[str('ipycInner')].append([
        "regionCoupledOFFBEAT", "bufferOuter", "false", 12, 13, 14, 15])
    else:
        patchDict[str('bufferOuter')].append(
        ["regionCoupledOFFBEAT", "ipycInner", "true", 4, 5, 6, 7])
        patchDict[str('ipycInner')].append([
        "regionCoupledOFFBEAT", "bufferOuter", "false", 8, 9, 10, 11])

# In case the Kernel is not modeled, the inner patch will be created
if Kernel == False:
    patchDict[str('inner')].append(["patch", "none", "false", 0, 1, 2, 3])

# Write boundaries
file.write("boundary \n(\n")
for patchName in patchDict:
    patchFaces = patchDict[patchName]
    writePatch(patchName, patchFaces, file)

file.write(");\n\n")
