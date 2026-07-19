# Power-temperature-momentum controller

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Multicode-FMI-red.svg)]()

## Description

This case demonstrates the use of an **FMU** as an external model for GeN-Foam.  
The FMU is a Modelica model that can receive and send information for:

- Pump momentum
- Heat-exchanger temperature (input and output)
- Heat-exchanger power (input and output)
- Heat-exchanger heat conductivity

## Model

The case consists of two channels:

- One with a fixed-temperature heat-exchanger `hx1` (green in Figure 1)
- One with a fixed-power heat-exchanger `hx2` (blue in Figure 1)

Two pumps are placed in the white and red zones (see Figure 1).

![](./images/powerTemperatureMomentumControl_mesh.png)

*Figure 1: GeN-Foam mesh with colorized cellZones.*

![](./images/powerTemperatureMomentumControl_momentumSourceTest.png)

*Figure 2: Modelica control model.*

The `Allrun` script is expected to output, at the end of the simulation, the same outlet temperature for both heat exchangers, using conservation of energy:

$$
q''' = \frac{\dot{m}}{\alpha V} c_p (T_{out} - T_{in})
$$

where $q'''$ is the power density imposed to `hx2`, $\alpha$ is the volume fraction of structure, $V$ is the volume of the heat-exchanger mesh, $\dot{m}$ is the mass flow rate, $c_p$ is the specific heat capacity of the fluid, $T_{out}$ is the measured `hx1` outlet temperature, and $T_{in}$ is the inlet temperature.

## FMU generation

Run the following command to build the FMU from the `momentumSourceTest.mo` file:
```bash
omc FMUGen.mos
```

## Run the case

```bash
./Allclean

./Allrun

# or
./Alltest
```

![](./images/powerTemperatureMomentumControl_temperatureDist.png)

*Figure 3: Temperature distribution using fix-temperature and fix-power FMI inputs.*

## Results visualization

```bash
python3 Allplot.py log.GeN-Foam
```