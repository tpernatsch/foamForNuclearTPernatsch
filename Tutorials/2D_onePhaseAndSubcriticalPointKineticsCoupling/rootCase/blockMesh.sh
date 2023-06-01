#!/bin/sh

cp system/blockMeshDict system/fluidRegion
cp system/blockMeshDict system/neutroRegion
cp system/blockMeshDict system/thermoMechanicalRegion

blockMesh -region fluidRegion
blockMesh -region neutroRegion
blockMesh -region thermoMechanicalRegion
