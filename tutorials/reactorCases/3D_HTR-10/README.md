# Pebble Bed Gas-Cooled Reactor Model (HTR-10 Design)

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]()


This tutorial provides an example of a pebble bed gas-cooled reactor, based on the HTR-10 design.


## Key Features of the Model

### Pebble-Wise Power Tallying

The power is tallied using a Voronoi mesh, ensuring accurate pebble-wise power capture. The Voronoï, written with the polyMesh format is shared with Serpent through the [`ifc` interface](http://serpent.vtt.fi/mediawiki/index.php?title=Input_syntax_manual#ifc). Most of the information about the pebbles is included in the data_peb.csv file while `data_mesh.csv` includes TH mesh information.


### Power Mapping

The power data is mapped to a finer thermal-hydraulic mesh ([`long_geo_HTR_ori.msh`](./rootCase/long_geo_HTR_ori.msh)) by overlapping the Voronoi and thermal-hydraulic meshes, ensuring power conservation.


### Packing Fraction Mapping

The equivalent mapping is done for the packing fraction, which represents the ratio of the pebble volume to the Voronoi cell volume. This packing fraction information is stored in the file [`./0/fluidRegion/alpha.structure`](./0/fluidRegion/alpha.structure).


### Defueling Chute Flow

There is an upward flow in the defueling chute, accounting for 3% of the total flow rate.


## Model Options

The model can be run with different packing fraction configurations:
- Pebble-wise packing fraction
- Homogeneous packing fraction
- Zoned packing fraction

The corresponding mappings can be generated using the [`PF_mapping.ipynb`](./PF_mapping.ipynb) notebook.


## Temperature Mapping

If you want to map temperatures back to Serpent, it is recommended to average the temperatures over several thermal-hydraulic cells corresponding to each pebble.


## Performance

This model was run on 6 nodes, each with 28 processors, achieving a total computation time of approximately 200 seconds.
