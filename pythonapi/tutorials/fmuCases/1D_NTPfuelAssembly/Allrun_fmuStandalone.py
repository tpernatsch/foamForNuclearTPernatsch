import fmpy
import numpy as np
import pandas as pd
import tqdm


# extract the FMU to a temporary directory
unzipdir = fmpy.extract("modelica/NTPTurbomachine.fmu")

# read the model description
model_description = fmpy.read_model_description(unzipdir)

# instantiate the FMU
fmu_instance = fmpy.instantiate_fmu(unzipdir, model_description)

end_time = 10
step_size = 0.001

gfT_in: float = 300
gfMassFlowInlet_in: float = 30
gfpNozzleChamber_in: float = 30e5
gfExtReact_in: float = 0
gfTnozzleChamber_out: float = 300
gfpInlet_out: float = 30e5
gfMassFlowOutlet_out: float = 30
gfPower_out: float = 10e3


fmu_instance.setupExperiment(stopTime=end_time)
fmu_instance.enterInitializationMode()
fmu_instance.exitInitializationMode()

valueRef = {
    variable.name: [variable.valueReference]
    for variable in model_description.modelVariables if variable.type == 'Real'
}

getVar = lambda key: fmu_instance.getReal(valueRef[key])[0]
setVar = lambda key, value: fmu_instance.setReal(valueRef[key], [value])

setVar("Tnozzle_in", gfTnozzleChamber_out)
setVar("pTurbine_in", gfpInlet_out)
setVar("mFlowCoreOutlet_in", gfMassFlowOutlet_out)
setVar("corePower_in", gfPower_out)


import psutil
process = psutil.Process()

outputFilename = "Turbomachine_standalone.csv"

# Extract all the parameters and variables composing the FMU
parameterList = [key for key in valueRef.keys() if key[:3] != "_D_"]

# Create the output buffer (write into a csv)
results = pd.DataFrame(columns=['time'] + parameterList)

writeInterval = 100*step_size
previous_time = 0
for t in tqdm.tqdm(np.arange(0, end_time-step_size, step_size)):
    # Make a step with the FMU
    fmu_instance.doStep(
        currentCommunicationPoint=t,
        communicationStepSize=step_size,
        noSetFMUStatePriorToCurrentPoint=True
    )
    if psutil.Process().memory_info().rss > 1 * 1024**3:
        # Save
        state = fmu_instance.getFMUstate()

        # Kill
        fmu_instance.terminate()
        fmu_instance.freeInstance()

        # Restart
        fmu_instance = fmpy.instantiate_fmu(unzipdir, model_description)
        fmu_instance.setFMUstate(state)
        fmu_instance.setupExperiment(startTime=t, stopTime=end_time)
        fmu_instance.enterInitializationMode()
        fmu_instance.exitInitializationMode()

        print("Restart")

    # Write the results
    if (previous_time + writeInterval < t):
        isFirstWrite = pd.isnull(results.index.max())
        results.loc[0] = [t] + [getVar(key) for key in parameterList]
        if (isFirstWrite):
            results.to_csv(outputFilename, index=False)
        else:
            results.to_csv(outputFilename, index=False, mode='a', header=False)

        previous_time = t


    # omega = getVar('Tturbine_out')
    # data.append(omega)

# free the FMU instance and unload the shared library
# fmu_instance.terminate()
# fmu_instance.freeInstance()

# state = fmu_instance.getFMUstate()
# fmu_instance.terminate()
# fmu_instance.freeInstance()


# # Restart
# fmu_instance = fmpy.instantiate_fmu(unzipdir, model_description)
# fmu_instance.setFMUstate(state)

# fmu_instance.setupExperiment(startTime=end_time, stopTime=2*end_time)
# fmu_instance.enterInitializationMode()
# fmu_instance.exitInitializationMode()

# for t in tqdm.tqdm(np.arange(end_time, 2*end_time-step_size, step_size)):
#     # Make a step with the FMU
#     fmu_instance.doStep(
#         currentCommunicationPoint=t,
#         communicationStepSize=step_size,
#         noSetFMUStatePriorToCurrentPoint=True
#     )

#     omega = getVar('Tturbine_out')
#     data.append(omega)


# import matplotlib.pyplot as plt

# plt.plot(data)
# plt.show()