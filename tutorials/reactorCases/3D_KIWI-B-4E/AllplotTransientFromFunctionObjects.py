"""
Simple script to extract point-kinetics parameters from the log.GeN-Foam file.

Author: Thomas Guilbaud, EPFL, 2023/07/17

From LA-3185-MS (page 66, 327):
- 20800 s :
    LH2 flow rate at 68.4 lb/s
    Thermal power = 905 MW
    Fluid at core inlet = 169 R, 573 psia
    Fluid at core outlet = 3520 R, 489 psia
- 20828 s : End of full power hold -> start transient
- 20850 s : Cooldown start
- 20858 s : LH2 flow rate at 28 lb/s
- 20910 s : (see page 314)
    Core power: Total = 52 MW, decay = 15 MW
    Core inlet at 57 R and 93 psia
    Chamber temperature at 484 R and 78.8 psia
    Core material at outlet at 523 R
    LH2 flow rate at 29.7 lb/s
    Reactivity: Total = -0.7$, hydrogen = 1.4$
- 20926 s : LH2 flow cut to 28 lb/s
- 20937 s : LH2 flow cut to 14 lb/s
- 20962.4 s : LH2 flow cut to 14 lb/s
- 20964 s : LH2 flow cut to 0 lb/s
- 20965 s : GH2 flow at
"""

#==============================================================================*
# Imports
#==============================================================================*

import sys
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd


#==============================================================================*
# User help
#==============================================================================*

if (len(sys.argv) < 3):
    print(f"\n    Usage: python3 {sys.argv[0]} path/to/transientFolder <timeFolder>\n")
    print("Don't forget to execute:")
    #print(">>> ./Allrunpostprocess_parallel path/to/transientFolder <timeFolder>")
    #print("or")
    #print(">>> ./Allrunpostprocess path/to/transientFolder <timeFolder>")
    #print("or")
    print('>>> reconstructPar -regions "(uniform neutroRegion fluidRegion)" -case transientShutdown -newTimes')
    print(">>> GeN-Foam -postProcess -case transientShutdown/ -region fluidRegion -dict system/functions")
    print(">>> postProcess -case transientShutdown/ -region fluidRegion -func sampleFluidDict\n")
    sys.exit(0)


#==============================================================================*
# Functions
#==============================================================================*

def applyTicker(ax):
    ax.xaxis.set_ticks(np.arange(0, 21000, 20))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))


#==============================================================================*
# Data extraction
#==============================================================================*

folderName: str = sys.argv[1] # "transientShutdown"
tFunction: float = sys.argv[2] # 20838


# Post process from postProcessing folder if functionObjects not used in controlDict
# >>> GeN-Foam -postProcess -case transientShutdown/ -region fluidRegion -dict system/functions
# >>> postProcess -case transientShutdown/ -region fluidRegion -func sampleFluidDict
# Remove TBulk and massFlow

fluidField = ["Time", "areaAverage(magU)", "areaAverage(p)", "areaAverage(T)", "areaAverage(alphaRhoPhi)", "areaAverage(alphaPhi)"]

dataInlet = pd.read_csv(
    f"{folderName}/postProcessing/fluidRegion/vInlet/{tFunction}/surfaceFieldValue.dat",
    sep='\s+', # Remove spaces between columns
    skiprows=5, # Not include column header because of "#"
    names=fluidField
)
dataOutlet = pd.read_csv(
    f"{folderName}/postProcessing/fluidRegion/vOutlet/{tFunction}/surfaceFieldValue.dat",
    sep='\s+', # Remove spaces between columns
    skiprows=5, # Not include column header because of "#"
    names=fluidField
)
dataPower = pd.read_csv(
    f"{folderName}/postProcessing/neutroRegion/totalPowerNeutro/{tFunction}/volFieldValue.dat",
    sep='\s+', # Remove spaces between columns
    skiprows=4, # Not include column header because of "#"
    names=["Time", "volIntegrate(powerDensity)"]
)

if (
    len(dataPower[dataPower['Time'] > 20000]) <= 1
    or np.isnan(dataInlet['areaAverage(alphaRhoPhi)'][0])
):
    dataInlet = pd.read_csv(
        f"{folderName}/postProcessing/fluidRegion/vInlet/{tFunction}/surfaceFieldValue_{tFunction}.dat",
        sep='\s+', # Remove spaces between columns
        skiprows=5, # Not include column header because of "#"
        names=fluidField
    )
    dataOutlet = pd.read_csv(
        f"{folderName}/postProcessing/fluidRegion/vOutlet/{tFunction}/surfaceFieldValue_{tFunction}.dat",
        sep='\s+', # Remove spaces between columns
        skiprows=5, # Not include column header because of "#"
        names=fluidField
    )
    dataPower = pd.read_csv(
        f"{folderName}/postProcessing/fluidRegion/totalPowerFluid/{tFunction}/volFieldValue_{tFunction}.dat",
        sep='\s+', # Remove spaces between columns
        skiprows=4, # Not include column header because of "#"
        names=["Time", "volIntegrate(powerDensity)"]
    )


# Extract axial temperature for station 26
timeFolders = [
    f"{t:.6g}" for t in sorted(
        [float(t) for t in os.listdir(f"{folderName}/postProcessing/sampleFluidDict/fluidRegion/")]
    )
]

radiusSamples = [0.0438, 0.0938, 0.1447, 0.1951, 0.2464, 0.2973, 0.3480, 0.3995]
timeStation26 = []
Tstation26 = []
for timeFolder in timeFolders:
    if (float(timeFolder) > 20000):
        tempTstation26 = []
        try:
            for radius in range(1, 9):
                axialFileName: str = f"FluidAxialR{radius}_T_Tmatrix.lumpedNuclearStructure_Tmax.lumpedNuclearStructure_Tsurface.lumpedNuclearStructure_magU_p_powerDensityNeutronics.xy"
                fieldNames = axialFileName.split('.xy')[0].split('_')[1:]
                dataAxialSample = pd.read_csv(
                    f"{folderName}/postProcessing/sampleFluidDict/fluidRegion/{timeFolder}/{axialFileName}",
                    sep='\s+', # Remove spaces between columns
                    names=['pos']+fieldNames
                )
                dataAxialSample = dataAxialSample[0 <= dataAxialSample['pos']]
                dataAxialSample = dataAxialSample[dataAxialSample['pos'] <= 0.001]

                tempTstation26.append(list(dataAxialSample['Tmatrix.lumpedNuclearStructure'])[0])

            # Compute radial average
            Tstation26.append(sum([T*r for T, r in zip(tempTstation26, radiusSamples)])/sum(radiusSamples))
            timeStation26.append(float(timeFolder))
        except:
            pass
            #print(f"Cannot extract samples at time {timeFolder}")

surface = 0.56602329 # from file

# Compute mass flow rate
dataInlet['rho'] = dataInlet['areaAverage(alphaRhoPhi)'] / dataInlet['areaAverage(alphaPhi)']
dataOutlet['rho'] = dataOutlet['areaAverage(alphaRhoPhi)'] / dataOutlet['areaAverage(alphaPhi)']
dataInlet['massFlow'] = surface * dataInlet['areaAverage(magU)'] * dataInlet['rho']

# Filter to shutdown transient
dataInlet = dataInlet[dataInlet['Time'] > 20000]
dataOutlet = dataOutlet[dataOutlet['Time'] > 20000]
dataPower = dataPower[dataPower['Time'] > 20000]


#==============================================================================*
# Extract reactivity
#==============================================================================*

time, totalPower, fissionPow, decayPow = [], [], [], []
totalRho, extRho, doppler, axial, radial, reactTcool, reactRhoCool = [], [], [], [], [], [], []
Tfuel, Tclad, Tcool, rhoCool = [], [], [], []

try:
    with open(folderName+"/log.GeN-Foam", 'r') as file:
        for line in file.readlines():
            if ("totalPower = " in line):
                totalPower.append(float(line.split()[2]))
            elif ("-> fission = " in line):
                fissionPow.append(float(line.split()[3]))
            elif ("-> decay   = " in line):
                decayPow.append(float(line.split()[3]))
            elif ("totalReactivity" in line):
                totalRho.append(float(line.split()[2]))
            elif ("-> extReactivity  =" in line):
                extRho.append(float(line.split()[3]))
            elif ("Doppler" in line):
                doppler.append(float(line.split()[4]))
            elif ("-> TFuel " in line):
                axial.append(float(line.split()[3]))
            elif ("-> TStruct " in line):
                radial.append(float(line.split()[3]))
            elif ("-> TCool " in line):
                reactTcool.append(float(line.split()[3]))
            elif ("-> rhoCool " in line):
                reactRhoCool.append(float(line.split()[3]))
            elif ("TFuel = " in line):
                Tfuel.append(float(line.split()[2]))
            elif ("TClad = " in line):
                Tclad.append(float(line.split()[2]))
            elif ("TCool = " in line):
                Tcool.append(float(line.split()[2]))
            elif ("rhoCool = " in line):
                rhoCool.append(float(line.split()[2]))
            elif ("Time = " == line[:7]):
                time.append(float(line.split()[2]))
                if (time[-1] > 21100):
                    break

except:
    print(f"{folderName}/log.GeN-Foam does not exist")


#==============================================================================*
# Plot the results
#==============================================================================*

figPower, axPower = plt.subplots(figsize=(5, 4), dpi=200) # Start at 20828
figTemperature, axTemperature = plt.subplots(figsize=(5, 4), dpi=200)
figMassflow, axMassflow = plt.subplots(figsize=(5, 4), dpi=200)
figReactivity, axReactivity = plt.subplots(figsize=(5, 4), dpi=200)

# Summary figure setup
figSummary, (axPowerSummary, axTemperatureSummary) = plt.subplots(nrows=2, figsize=(5.5, 8), dpi=400)
axReactivitySummary = axPowerSummary.twinx()
axMassflowSummary = axTemperatureSummary.twinx()

axPowerSummary.set_xlabel('Time [s]')
axPowerSummary.set_ylabel("Thermal Power [MW]")
axPowerSummary.set_ylim((0, 1000))

axReactivitySummary.set_ylabel(r'Reactivity $\rho$ [pcm]')
axReactivitySummary.set_ylim((-2800, 1500))

axTemperatureSummary.set_xlabel('Time [s]')
axTemperatureSummary.set_ylabel('Temperature $T$ [K]')
axTemperatureSummary.yaxis.set_ticks(np.arange(0, 2500, 500))
axTemperatureSummary.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
axTemperatureSummary.set_ylim((-750, 2000))

axMassflowSummary.set_ylabel("Mass flow rate $\dot{m}$ [kg/s]")
axMassflowSummary.yaxis.set_ticks(np.arange(10, 40, 10))
axMassflowSummary.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
axMassflowSummary.set_ylim((10, 40))


# Custom x-ticks
applyTicker(axPower)
applyTicker(axTemperature)
applyTicker(axReactivity)
applyTicker(axPowerSummary)
applyTicker(axTemperatureSummary)


#==============================================================================*
# Power
#==============================================================================*

axPower.plot(
    dataPower['Time'],
    [power/1e6 for power in dataPower['volIntegrate(powerDensity)']],
    label="GeN-Foam"
)
# For summary
pPowerGenfoam, = axPowerSummary.plot(
    dataPower['Time'],
    [power/1e6 for power in dataPower['volIntegrate(powerDensity)']],
    label="Power GeN-Foam",
    color='tab:red',
    marker='o', markevery=5, markersize=4
)

timeExpTotalPower = [
    20837.8797806378, 20841.3680179241, 20844.7231268237, 20847.4080345977,
    20849.2882451021, 20850.7615477803, 20852.8993525948, 20854.2372472905,
    20855.0419445602, 20856.5236817424, 20858.4729639954, 20860.3536304189,
    20862.7036655898, 20865.7911499583, 20867.8703691642, 20871.2204629533,
    20874.8324822843, 20879.0431234368, 20883.7206257816, 20888.1302756357,
    20890.2012882972
]
expTotalPower = [
    531.70034799815, 432.622193867089, 349.783645553741, 291.767904079125,
    254.250964692217, 241.148344673356, 239.821134488679, 233.243475421761,
    223.26591482428, 188.509600167652, 159.190616854695, 137.851740905784,
    116.42835850457, 94.7192169838883, 84.1173086510379, 72.8801773153253,
    67.8849739510681, 65.4709064009655, 63.7516537064256, 62.6597024844917,
    62.3133465069776
]
timeExpFissionPower = [
    20822.1979250018, 20824.255423806, 20825.6494399953, 20826.9116210329,
    20828.241306383, 20829.4391565814, 20830.8383863341, 20832.7726637002,
    20834.6401167838, 20836.84169128, 20839.5764999762, 20842.2449377433,
    20845.3129611454, 20848.1159541841, 20849.1837826443, 20849.717923551,
    20850.2511577511, 20851.0494223145, 20852.1115838579, 20853.0396434829,
    20853.8363213096, 20854.501617338, 20856.1731311075, 20857.0423001335,
    20858.4465167729, 20860.6512647425, 20862.6546331578, 20864.8571143606,
    20868.3247233978, 20871.7887056082, 20875.1129189835, 20880.491911043,
    20885.403813155, 20889.6517792985
]
expFissionPower = [
    930.508041386843, 930.809206461036, 928.089308411643, 913.787858080139,
    885.676265226114, 847.678452365787, 786.212300643977, 689.120084798204,
    607.823738473734, 519.543633373145, 430.369631203254, 356.497653155966,
    289.802618515234, 234.09793978393, 215.750298401745, 206.472495586059,
    200.095548912406, 195.147500216208, 194.567108255079, 197.680311801365,
    197.084157459385, 192.812114436133, 162.733887726621, 149.034221472047,
    128.985350941052, 105.501660035986, 89.0483021761169, 75.1632576138329,
    59.9627852106056, 50.3055238408666, 46.5253591191142, 44.6989775036693,
    44.452978009723, 44.3429826812366
]

expFissionPowerInterp = np.interp(timeExpTotalPower, timeExpFissionPower, expFissionPower)

expDecayPower = [tot-fiss for tot, fiss in zip(expTotalPower, expFissionPowerInterp)]

# plt.plot( # Experimental data from LANL Report LA-3185-MS, Chapter XVII, Figure 2, page 334
#     timeExpTotalPower,
#     expTotalPower,
#     ls='--',
#     label="Experimental (Total Shutdown)"
# )
axPower.plot( # Experimental data from LANL Report LA-3185-MS, Chapter XVII, Figure 2, page 334
    timeExpFissionPower,
    expFissionPower,
    ls='--',
    label="Experimental" # (Fission Shutdown)"
)
# For summary
pPowerExp, = axPowerSummary.plot(
    timeExpFissionPower,
    expFissionPower,
    ls='--',
    marker='x',
    markevery=3,
    color='black',
    label="Power Exp." # (Fission Shutdown)"
)
# plt.plot(
#     timeExpTotalPower,
#     expDecayPower,
#     ls='--',
#     label="Experimental (Decay Shutdown)"
# )
axPower.set_xlabel("Time [s]")
axPower.set_ylabel("Thermal Power [MW]")
# axPower.set_ylim((1e-2, 1e3))
axPower.set_ylim((0, 1000))
# plt.legend(bbox_to_anchor=(0, 1.01, 1, 0.2), loc="lower left", mode='expand', ncol=2)
axPower.legend()


#==============================================================================*
# Temperature
#==============================================================================*

# axTemperature.plot(
#     dataInlet['Time'],
#     dataInlet['areaAverage(T)'],
#     label="Inlet"
# )

# Experimental data from LANL Report LA-3185-MS, Chapiter III, Figure 10, page 63
timeExpStation26 = [
    20822.6459790389, 20825.9590417786, 20828.697520356, 20830.1223414074,
    20830.7821441024, 20831.471993765, 20832.9894169042, 20834.8642246242,
    20836.9568472219, 20840.0407121029, 20842.575117419, 20844.889451771,
    20847.0378612302, 20849.3519904834, 20851.608897184, 20853.8071456406,
    20856.7150357897, 20858.8521648174, 20860.3879442952, 20861.7069343889,
    20863.4145866235, 20865.8337264393, 20868.8002748323, 20871.9284410444,
    20875.3837397707, 20879.0545972886, 20882.4506224747, 20885.5170539615,
    20890.0056402158
]
temperatureExpStation26 = [
    1785.2989239822, 1786.74658051972, 1784.82206899522, 1777.64189766803,
    1757.64304389524, 1717.37016196649, 1630.64237485489, 1533.146894489,
    1445.10989261029, 1315.37115299957, 1199.43727980393, 1091.20999849457,
    995.394620107881, 888.727799741575, 800.762582484566, 742.422194938162,
    701.557944459257, 691.567917939987, 673.512376425948, 638.196051709497,
    562.475239924523, 490.18641671965, 419.697336530138, 369.566033221895,
    330.501525758351, 318.060567155484, 313.302272231678, 316.202712779702,
    332.208974260104
]

# Experimental data from LANL Report LA-3185-MS, Chapiter III, Figure 10, page 63
timeExpOutlet = [
    20824.8428855157, 20827.4301540999, 20829.8714356646, 20831.7326704142,
    20833.6024391758, 20834.7599615194, 20835.7675955809, 20836.7805938786,
    20838.7054901404, 20839.6499503476, 20844.8224038871, 20849.7859403816,
    20859.9014487871, 20864.8834054476, 20869.9596352319, 20874.7647494325,
    20879.8839152717, 20884.8536804866, 20889.7928340899
]
temperatureExpOutlet = [
    1908.26865670308, 1904.45507095333, 1893.86943962264, 1829.38329071499,
    1706.24515844091, 1629.69759637537, 1573.19032398213, 1479.81609061569,
    1381.84210734856, 1355.47112740814, 1160.14381678252, 986.504982075472,
    758.303596822244, 660.088247368422, 520.020001142007, 425.073439473685,
    393.957889920556, 379.531088133069, 373.46934875869
]

# Plot
axTemperature.plot(timeStation26, Tstation26, ls='-.', label="Station 26 GeN-Foam")
axTemperature.plot(
    timeExpStation26, temperatureExpStation26,
    ls='--',
    label="Station 26 Exp."
)
axTemperature.plot(dataOutlet['Time'], dataOutlet['areaAverage(T)'], label="Nozzle Chamber GeN-Foam")
axTemperature.scatter(
    timeExpOutlet, temperatureExpOutlet,
    marker='x',
    label="Nozzle Chamber Exp.",
    color='black'
)

# For summary
pTStation26Genfoam, = axTemperatureSummary.plot(
    timeStation26, Tstation26, ls='-', label=r"$T_{Station 26}$ GeN-Foam"
)
pTStation26Exp, = axTemperatureSummary.plot(
    timeExpStation26, temperatureExpStation26,
    ls=':',
    label=r"$T_{Station 26}$ Exp."
)
pTOutletGenfoam, = axTemperatureSummary.plot(
    dataOutlet['Time'], dataOutlet['areaAverage(T)'],
    label=r"$T_{Chamber}$ GeN-Foam",
    ls='-.'
)
sTOutletExp = axTemperatureSummary.scatter(
    timeExpOutlet, temperatureExpOutlet,
    label=r"$T_{Chamber}$ Exp.",
    color='black',
    marker='x'
)

# axTemperature.hlines([83.3, 1972], xmin=min(time), xmax=max(time), ls='--', color='grey')
axTemperature.set_xlabel("Time [s]")
axTemperature.set_ylabel("Temperature [K]")
axTemperature.set_ylim(0)
axTemperature.legend(loc='upper right')


#==============================================================================*
# Mass flow rate
#==============================================================================*

# Experimental data from LANL Report LA-3185-MS, Chapiter III, Figure 9, page 61
timeExpMassFlowRate = [
    20822.9602774092, 20830.4146364074, 20832.1971992561, 20834.1850959171,
    20835.7713343801, 20838.8368569558, 20843.3830541666, 20847.558850965,
    20849.7798629177, 20852.8192262901, 20854.8809654861, 20856.8160802331,
    20858.5171526736, 20860.6469695302, 20865.3325666147, 20868.1548057035,
    20870.7108176851, 20873.7990521272, 20879.230548619, 20882.7979918539,
    20887.430343517, 20889.9328783235, 20926, 20937, 20962.4, 20964, 20991,
    21013
]
massFlowRateExp = [
    31.082840078101, 31.0423682460298, 30.9092973877877, 29.9607191133102,
    28.7675359187008, 26.332393087713, 22.8033544527425, 19.5231124633684,
    17.8822687573944, 15.4966250794625, 13.9060032378315, 13.0070686586673,
    12.7016902894025, 12.6901269088107, 12.6646874715088, 12.6740432976239,
    12.6848445463131, 12.668077644455, 12.6879456347445, 12.6685769722533,
    12.6434266194661, 12.6298396472708, 12.6298396472708, 6.35029, 6.35029,
    0.907185, 3.85554, 2.58548
]

axMassflow.plot(
    dataInlet["Time"],
    dataInlet["massFlow"],
    label="GeN-Foam"
)
axMassflow.plot(
    timeExpMassFlowRate, massFlowRateExp,
    ls='--',
    label="Experimental"
)
axMassflow.set_xlabel("Time [s]")
axMassflow.set_ylabel("Mass flow rate [kg/s]")
axMassflow.set_ylim(0)
axMassflow.legend(loc='upper right')

pMassFlowGenfoam, = axMassflowSummary.plot(
    dataInlet["Time"],
    dataInlet["massFlow"],
    label=r"$\dot{m}$ GeN-Foam",
    marker='o', markersize=4, markevery=5,
    color='tab:red'
)
pMassFlowExp, = axMassflowSummary.plot(
    [t for t in timeExpMassFlowRate if t < float(tFunction)],
    [massFlow for t, massFlow in zip(timeExpMassFlowRate, massFlowRateExp) if t < float(tFunction)],
    ls='--',
    label=r"$\dot{m}$ Exp.",
    color='black'
)


#==============================================================================*
# Reactivity
#==============================================================================*

time = time[:len(totalRho)]
axReactivity.plot(time, totalRho, label="Total")
axReactivity.plot(time, doppler, label="Doppler")
axReactivity.plot(time, axial, label="Axial")
# axReactivity.plot(time, radial, label="Radial")
# axReactivity.plot(time, reactTcool, label="T Coolant")
axReactivity.plot(time, reactRhoCool, label=r"$\rho$ Coolant")
axReactivity.plot(time, extRho, label="Control Drum from Exp.")
axReactivity.set_xlabel("Time [s]")
axReactivity.set_ylabel("Reactivity [pcm]")
axReactivity.legend(loc='upper right')

# For summary
pRhoTot, = axReactivitySummary.plot(time, totalRho, label=r"$\rho_{Total}$", ls='-')
pRhoDoppler, = axReactivitySummary.plot(time, doppler, label=r"$\rho_{Doppler}$", ls='--')
pRhoAx, = axReactivitySummary.plot(time, axial, label=r"$\rho_{Axial}$", ls='-.')
pRhoCool, = axReactivitySummary.plot(time, reactRhoCool, label=r"$\rho_{Coolant}$", ls=':')
pRhoExt, = axReactivitySummary.plot(time, extRho, label=r"$\rho_{Control Drum}$ from Exp.", ls=(0, (3, 1, 1, 1, 1, 1)))


#==============================================================================*

# Legend and grid
axPowerSummary.legend(
    handles=[pPowerGenfoam, pRhoTot, pRhoDoppler, pRhoAx, pPowerExp, pRhoCool, pRhoExt],
    bbox_to_anchor=(0, 1.01, 1, 0.2),
    loc="lower left",
    mode='expand',
    ncol=2
)
axTemperatureSummary.legend(
    handles=[pTStation26Genfoam, pTStation26Exp, pMassFlowGenfoam, pTOutletGenfoam, sTOutletExp, pMassFlowExp],
    bbox_to_anchor=(0, 1.01, 1, 0.2),
    loc="lower left",
    mode='expand',
    ncol=2
)

# Tight layout
for fig in [figPower, figTemperature, figMassflow, figReactivity, figSummary]:
    fig.tight_layout()

# Save figures
figPower.savefig("fig_results_transientPower.png")
figTemperature.savefig("fig_results_transientTemperature.png")
figMassflow.savefig("fig_results_transientMassflowrate.png")
figReactivity.savefig("fig_results_transientReactivity.png")
figSummary.savefig("fig_results_transientSummary.png")


#==============================================================================*
