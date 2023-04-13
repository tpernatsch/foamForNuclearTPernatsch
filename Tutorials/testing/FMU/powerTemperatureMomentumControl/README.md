# Power-temperature-momentum controller

This case is a demonstration of use of the FMU as an external model for GeN-Foam.
The FMU is a Modelica model that is able to receive and send information from:

- Core power
- Pump momentum
- Heat-exchanger temperature
- Heat-exchanger heat conductivity


## FMU generation

Run the following command to build the FMU from the `momentumSourceTest.mo` file:
```bash
omc FMUGen.mos
```


## Run the case

```bash
./Allclean

runApplication GeN-Foam
```


## Results vizualisation

```bash
python3 plot.py log.GeN-Foam
# or 
python3 plotMassFlowRate.py log.GeN-Foam
```