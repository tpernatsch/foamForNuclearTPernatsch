# 3D Nuclear Fuel Rod Mesh Generator

This Python script, `rod3Dmaker.py`, is designed to create a mesh for a 3D nuclear fuel rod using blockMesh. The mesh is generated based on the user's inputs provided in a JSON file named `rod3Ddict.json`. The resulting mesh information is stored in a `blockMeshDict` file.

## Prerequisites

Before running the script, you need to have the following software and libraries installed:

- JSON library - Most Python installations come with a built-in JSON library, so no additional installation should be necessary.

## Usage

1. Modify the JSON input file named `rod3Ddict.json` in the same directory as the script. The JSON file should contain the necessary parameters for mesh generation. Below is an example JSON structure:

```json
{
    "unit" : "mm",
    "fuelInnerRadius" : [0.825000],
    "fuelOuterRadius" : [3.780000],
    "cladInnerRadius" : [3.877500],
    "cladOuterRadius" : [4.560000 ],
    "rodHeight" : [2],
    "nRadialCellsFuel":[20],
    "nRadialCellsClad": [3],
    "nAzimuthalCellsPerQuarter":[20],
    "nZetaCells": [1],
    "center": [0.0,0.0,1209]
}
```

3. Run the script by executing the following command in your terminal:

```bash
python rod3Dmaker.py
```

4. The script will read the input parameters from `rod3Ddict.json` and generate a `blockMeshDict` file in the same directory.
