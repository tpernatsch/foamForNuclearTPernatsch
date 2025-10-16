#!/bin/bash
# Command to build an application for SerpentToFoamXS
# Must be executed in the current folder.
# Author: Thomas Guilbaud, EPFL/Transmutex SA, 18/03/2022

# Create build directory
buildFolder=build
mkdir $buildFolder

# Copy the source file to the build directory
cp SerpentToFoamXS.py $buildFolder
cp -r src $buildFolder
cd $buildFolder

# Install PyInstaller to generate the SerpentToFoamXS executable
python3 -m pip install pyinstaller

# Build the Application SerpentToFoamXS
pyinstaller SerpentToFoamXS.py src/*.py src/*/*.py --onefile --distpath bin/ --clean --hidden-import='PIL._tkinter_finder'

# End
echo "End of building"
