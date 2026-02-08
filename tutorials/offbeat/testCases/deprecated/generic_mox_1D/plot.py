import numpy as np
import matplotlib.pyplot as plt

timeDependentOutputs = {
    "rodInternalPressure":{
        "path":"postProcessing/rodPressure/0/writePressure.dat",
        "name":"Rod Internal Pressure",
        "unit":"bar",
        "scale":1e-5
    },
    "centralVoidRadius":{
        "path":"postProcessing/centralVoidRadius/0/centralVoid_PPN.dat",
        "name":"Central Void Radius",
        "unit":"mm",
        "scale":1e3
    },
    "gapSize":{
        "path":"postProcessing/gapSize_PPN/0/gapWidth",
        "name":"Radial Gap Size",
        "unit":r"$\mu m$",
        "scale":1e6
    }
}

for output, dict in timeDependentOutputs.items():
    
    # Read data
    data = np.genfromtxt(dict["path"], comments="#")

    # Scale data if needed
    if "scale" in dict.keys():
        data[:,1] *= dict["scale"]

    # Scale time
    data[:,0] /= 3600*24

    # Create plot
    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(10*cm,8*cm))

    ax.plot(data[:,0], data[:,1], "-k")

    ax.set_xlabel("Time (days)")
    ax.set_ylabel(f"{dict['name']} ({dict['unit']})")

    ax.grid(alpha=0.2)

    ax.tick_params(direction='in')
    fig.tight_layout()
    fig.savefig(f"plt_{output}")
    plt.close()
