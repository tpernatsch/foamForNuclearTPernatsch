# Post-processing {#postprocessing}

![Image title](https://upload.wikimedia.org/wikipedia/commons/2/20/UnderCon_icon.svg)
**This page is a work in progress!!!**

## Graphical Visualisation using ParaView
Several tutorials focused on ParaView usage for visualising CFD simulation results. This guide assumes that the reader has a working understanding of how to use ParaView to postprocess general OpenFOAM simulation results. A fairly comprehensive overview is provided in the [OpenFOAM User Guide](https://cfd.direct/openfoam/user-guide/v9-paraview/).
The recommended approach for visualising the solution fields from OFFBEAT simulations is to use `paraFoam`, the ParaView plugin that is shipped with OpenFOAM. Even though ParaView has its own builtin OpenFOAM reader, we recommend `paraFoam` primarily because it has better support for cases containing `cellZones` (OFFBEAT uses `cellZones` to define different material regions in the mesh). `paraFoam` allows us to focus only on the fuel or cladding materials, or to view the gap width or the fuel-cladding contact pressure.

The user launches ParaView from the commandline as follows:
```sh
paraFoam -case <case_directory>
```

In some cases `paraFoam` is not available. This may be the case if, for example, the user is running simulations on a cluster. In this case, the ParaView builtin OpenFOAM reader can be used:
* In the case folder, create the file `.foam`
* Run ParaView and open the `.foam` file.

### List of Fields

The list of fields available for visualization using ParaView is quite extensive. The following is a brief overview of these:

- **Basic solution fields:** displacement vector `D` [m], temperature `T` [K].
- **Strain:** total strain tensor `epsilon` [-] and several strain components, e.g., creep strain `epsilonCreep` [-].
- **Stress:** stress tensor `sigma` [Pa], von Mises stress `sigmaEq` [Pa], hydrostatic stress `sigmaHyd` [Pa], and yield stress `sigmaY` [Pa].
- **Burnup** `BU` in MWd/MTU and **fast fluence** `fastFluence`.
- **Gap width** `gapWidth` [m] and **gap contact pressure** `interfaceP` [Pa], which are zero everywhere except on the boundaries of the gap.

#### TUBRNP-related fields
- Isotope concentrations, e.g., `N_Am241` and `N_Pu240`
- Burnup in MWd/kgU, `BUkgU`
- Radial form factor, `formFactorTUBRNP`

#### Sciantix-related fields
- `Grain_radius`
- `Gas_produced`
- `Gas_grain`
- `Gas_grain_solution`
- `Gas_grain_bubbles`
- `Gas_boundary`
- `Gas_released`
- `Effective_burn_up`
- `Oxygen_to_metal_ratio`
- `Helium_produced`
- `Helium_grain`
- `Helium_grain_solution`
- `Helium_grain_bubbles`
- `Helium_boundary`
- `Helium_released`
- `Intragranular_bubble_concentration`
- `Intragranular_bubble_radius`
- `Intragranular_gas_swelling` (tensor)
- `Intergranular_bubble_concentration`
- `Intergranular_atoms_per_bubble`
- `Intergranular_vacancies_per_bubble`
- `Intergranular_bubble_radius`
- `Intergranular_bubble_area`
- `Intergranular_bubble_volume`
- `Intergranular_fractional_coverage`
- `Intergranular_saturation_fractional_coverage`
- `Intergranular_gas_swelling` (tensor)
- `Intergranular_fractional_intactness`

### Tips and Tricks for Graphical Post-processing
#### Dealing with the high aspect ratio of fuel rod models.
In general, fuel rod models have very high aspect ratio, i.e. height >> radius. As a consequence, it can be very difficult to view the results in ParaView since they appear as extremely thin vertical lines and radial behaviour cannot be easily distinguished. In these cases, the user should use the **Transform** filter and scale the vertical-/z-axis by an appropriate factor. For typical LWR fuel rods with a \f$ radius/height \approx 1e-3 \f$, a z-axis scaling of \f$ \approx 5e-3 \f$ generally works well. The option *Tranform all Input Vectors* should be disabled so that vector and tensor fields are not distorted by the filter.
If the user wishes to visualise the displacement of the model and the closing of the gap, the **Warp By Vector** filter can be used with the displacement field `D` or `DD`. A warp scaling factor of 1 will give a realistic representation of the actual configuration, while larger values can be used to visualise regions of high deformation. Note that this filter should be applied before the **Transform** filter or there will be an unrealistic distortion of the model in the z-direction.

An example of visual post-processing using ParaView with the **Warp By Scalar** and **Transform** filter applied is shown in the figure below.

#### Saving ParaView States
You can save yourself some time repeatedly opening the case and applying these filters each time you run ParaView by saving a ParaView [state file](https://www.paraview.org/Wiki/Advanced_State_Management#Save_State_/_Load_State).


\image html paraview_transform.png "Screenshot of ParaView showing the Transform filter applied to the 'generic_pwr2D' tutorial case" width=50%


## Post-processing using the OpenFOAM Utilities
Graphical post-processing using ParaView (see previous section) is useful for qualitative analysis, i.e. to identify trends and areas of interest for further analysis, but is not ideal for quantitative analysis of results. OpenFOAM provides several utilities for sampling and extracting data. These are preferred over ParaView's sampling filters since they are better adapted for OpenFOAM models. Comprehensive descriptions of the available utilities are available in the [6.2](https://cfd.direct/openfoam/user-guide/v9-post-processing-cli) and [6.3](https://cfd.direct/openfoam/user-guide/v9-graphs-monitoring) of the OpenFOAM User Guide.

Of particular relevance in fuel behaviour applications are the following:

### Probing Data
Probes can be defined and enabled in the model `controlDict` to extract time-dependent point values at selected locations in the mesh, see [6.3.1](https://cfd.direct/openfoam/user-guide/v9-graphs-monitoring) of the OpenFOAM User Guide. These are useful for validation exercises where, for example, thermocouple measurements are available. Examples of probes can be found in the *generic_pwr2D* test case.

As an example, to probe the fuel centreline temperature and burnup at three axial locations:

In `system/probes`:
```c

points
(
    (0.0 0.0 1.0)   // Location of probe (x,y,z)
    (0.0 0.0 1.5)
    (0.0 0.0 2.0)
);

fields  (T BU);

#includeEtc "caseDicts/postProcessing/probes/probes.cfg"

```

In `system/controlDict`:
```c

functions 
{ 
    #includeFunc  probes
}

```

While executing, OFFBEAT will create separate text files for each field temperature `T` and burnup `BU` in `postProcessing/probes/0`. Each text file will contain the probed values sampled at each time point in a simple ASCII table, which can be readily plotted using your favourite graphing tool, e.g.

```
# Probe 0 (0 0 1)
# Probe 1 (0 0 1.5)
# Probe 2 (0 0 2)
# Time        0             1             2            
0             300           300           300          
1             332.354       332.351       332.345      
302.139       661.691       657.762       646.909      
3600          1571.06       1503.88       1322.83      
6241.26       1572.51       1505.23       1323.91      
```

### Line Graphs
Radial and axial profiles of, for example, temperature `T` can be extracted using the sampling tools, see [6.3.2](https://cfd.direct/openfoam/user-guide/v9-graphs-monitoring) of the OpenFOAM User Guide.

As an example, radial and axial profiles of temperature can be extracted at each output time as follows:

In `system/sampleLines`:
```c

type                sets;
interpolationScheme cellPointFace;
setFormat           csv;
libs                ("libsampling.so");

sets
(
    radialProfile
    {
        type    lineCellFace;
        axis    distance;
        start   (0 0 1.6);
        end     (5.6e-3 0 1.6);
    }
    axialProfile
    {
        type    lineCellFace;
        axis    distance;
        start   (0 0 0);
        end     (0 0 3.2);
    }
);

fields (T);


```

After the simulation is completed, the following can be run from the command line:
```sh

postProcess -func sampleLines

```

This will create separate comma-separated value (CSV) text files `postProcessing/sampleLines/<time>/axialProfile_T.csv` and `radialProfile_T.csv` for each output time, which can be readily plotted using your favourite graphing tool.

***

Return to [User Manual](@ref userManual)
