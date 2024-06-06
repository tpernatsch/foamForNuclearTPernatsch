# Residual plot

This is a gnuplot utility (which should come as a default application on Ubuntu systems) for plotting the residuals of various GeN-Foam log-printed variables in real time. By default, it plots all available residuals from the log file.


## Usage

To use it, place it in the same folder as the log.GeN-Foam log file and run:

```bash
gnuplot residuals -
```

from the command line. 


## Additional options

To specify a different log name and/or path, run it with:

```bash
gnuplot -e "log='pathToYourLog'" residuals -
```

-------------------------------------------------------------------------------

The two-phase enthalpy residuals are called `h.fluid1Name`, `h.fluid2Name` in the log file. As the residuals script requires knowledge of the phase names, these can be passed via the fluid1, fluid2 arguments as follows:

```bash
gnuplot -e "fluid1='fluid1Name';fluid2='fluid2Name'" residuals -
```

If not provided, these default to 'liquid' and 'vapour respectively'.

