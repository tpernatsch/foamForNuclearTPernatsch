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

import sys
import matplotlib.pyplot as plt


#==============================================================================*
# User help

if (len(sys.argv) < 2):
    print(f"\n    Usage: python3 {sys.argv[0]} path1/to/log.GeN-Foam path2/to/log.GeN-Foam ...\n")
    sys.exit(1)


#==============================================================================*

time, totalPower, fissionPow, decayPow = [], [], [], []
totalRho, extRho, doppler, axial, radial, reactTcool, reactRhoCool = [], [], [], [], [], [], []
Tmin, Tavg, Tmax = [], [], []
Tfuel, Tclad, Tcool, rhoCool = [], [], [], []
Tinlet, Toutlet = [], []
massFlow = []

# Loop over the log files
for filename in sys.argv[1:]:
    print(f"Analyze {filename} ...")

    try:
        with open(filename, 'r') as file:
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
                elif ("patch inlet TBulk" in line):
                    Tinlet.append(float(line.split()[4]))
                elif ("patch outlet TBulk" in line):
                    Toutlet.append(float(line.split()[4]))
                elif ("T (avg min max)" in line):
                    line = line.split()
                    Tavg.append(float(line[5]))
                    Tmin.append(float(line[6]))
                    Tmax.append(float(line[7]))
                elif ("patch inlet massFlow =" in line):
                    massFlow.append(float(line.split()[4]))
    except:
        print(f"{filename} does not exist")




if (len(time) == 0):
    print("No data extracted, end of script.")
    sys.exit(0)


# Remove last time if necessary
while (len(totalPower) < len(time)):
    time = time[:-1]

tStartTransient = 20820

# Plot the results

plt.figure("Power", figsize=(5, 4), dpi=200) # Start at 20828
plt.plot(
    time, # [t-tStartTransient for t in time],
    [power / 1e6 for power in totalPower],
    label="Total"
)
# plt.plot(
#     time, # [t-tStartTransient for t in time],
#     [power / 1e6 for power in fissionPow],
#     label="Fission"
# )
# plt.plot(
#     time, # [t-tStartTransient for t in time],
#     [power / 1e6 for power in decayPow],
#     label="Decay"
# )

# timeExpTotalPower = [
#     20837.8797806378, 20841.3680179241, 20844.7231268237, 20847.4080345977,
#     20849.2882451021, 20850.7615477803, 20852.8993525948, 20854.2372472905,
#     20855.0419445602, 20856.5236817424, 20858.4729639954, 20860.3536304189,
#     20862.7036655898, 20865.7911499583, 20867.8703691642, 20871.2204629533,
#     20874.8324822843, 20879.0431234368, 20883.7206257816, 20888.1302756357,
#     20890.2012882972
# ]
# expTotalPower = [
#     531.70034799815, 432.622193867089, 349.783645553741, 291.767904079125,
#     254.250964692217, 241.148344673356, 239.821134488679, 233.243475421761,
#     223.26591482428, 188.509600167652, 159.190616854695, 137.851740905784,
#     116.42835850457, 94.7192169838883, 84.1173086510379, 72.8801773153253,
#     67.8849739510681, 65.4709064009655, 63.7516537064256, 62.6597024844917,
#     62.3133465069776
# ]
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

# expFissionPowerInterp = np.interp(timeExpTotalPower, timeExpFissionPower, expFissionPower)

# expDecayPower = [tot-fiss for tot, fiss in zip(expTotalPower, expFissionPowerInterp)]

# # def decayFitFunc(x, a, b, c):
# #     return a * np.exp(-b * (x-c))

# # popt, pcov = curve_fit(decayFitFunc, [t-tStartTransient for t in timeExpTotalPower], expDecayPower, bounds=(0, [60, 1/20, 20]))

# plt.plot( # Experimental data from LANL Report LA-3185-MS, Chapter XVII, Figure 2, page 334
#     timeExpTotalPower, # [t-tStartTransient for t in timeExpTotalPower],
#     expTotalPower,
#     ls='--',
#     label="Experimental (Total Shutdown)"
# )
plt.plot( # Experimental data from LANL Report LA-3185-MS, Chapter XVII, Figure 2, page 334
    timeExpFissionPower, # [t-tStartTransient for t in timeExpFissionPower],
    expFissionPower,
    ls='--',
    label="Experimental (Fission Shutdown)"
)
# plt.plot(
#     timeExpTotalPower, # [t-tStartTransient for t in timeExpTotalPower],
#     expDecayPower,
#     ls='--',
#     label="Experimental (Decay Shutdown)"
# )

# # plt.plot(
# #     [t-tStartTransient for t in time],
# #     [decayFitFunc(t-tStartTransient, *popt) for t in time],
# #     'r-',
# #     label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
# # )

# if (len(time) > 0 and max(time) > 209000):
#     plt.plot(
#         [t-tStartTransient for t in [
#          20851.0001311515, 20851.5216810494, 20852.2397646473, 20853.4841776556,
#          20854.907762341, 20856.1229418327, 20857.5571598078, 20860.4201493022,
#          20863.2971490185, 20866.5985755744, 20871.7733729473, 20878.2356670769,
#          20890.6568713372, 20903.004886131, 20907.9033036382, 20910.5347953853,
#          20914.7342973716, 20917.3188029441, 20921.5504621451, 20928.1662936304,
#          20933.5398258335, 20936.3967735885, 20944.4835734444, 20952.2110099508,
#          20961.7729789723, 20974.2260905334, 20999.1628362864, 21024.1573288362,
#          21084.8266075061, 21166.6592158714, 21238.6803786148, 21340.5255089433,
#          21445.5323893577, 21610.0446031019, 21830.8622633804, 21851.3199591816,
#          22345.3640372776, 23291.854864465, 24186.3547847777, 25459.5819916949,
#          27290.4665521563, 29699.7799453018, 32876.9351989888, 36745.3597181598,
#          41858.7646050091, 53311.3007372218, 65960.4394553233, 81478.5782278182,
#          101427.027561752]],
#         [243.080401775879, 241.338861298086, 239.53179466761, 232.545065368636,
#          219.251384369495, 200.712558411584, 171.663017798266, 134.257882431574,
#          108.546127590241, 86.752259361002, 70.9556503355749, 64.6140931719678,
#          62.3281681798919, 59.7255768330461, 59.1049181268265, 55.2496190966926,
#          45.0986439446288, 39.8360734694229, 37.6674922041408, 36.8566757882959,
#          34.8521770978535, 30.7835984964583, 23.2194093290841, 16.7371137172998,
#          12.7680966198811, 10.664805798128, 8.9147606375939, 8.24846036583876,
#          6.9027881285238, 5.647400283847, 4.82977834738185, 3.9488298290438,
#          3.3012034386284, 2.59455152977054, 1.9937551552353, 2.0168515576565,
#          1.32467864439411, 0.764570247867596, 0.533951052546154,
#          0.368738965026171, 0.248975398218595, 0.169995209393258,
#          0.115401864413866, 0.0810236854620977, 0.0565658483429295,
#          0.0322634802638148, 0.0212959767125009, 0.0145372853699244,
#          0.0101496075645513],
#         ls='--',
#         label="Experimental (Cooldown)" # LA-3185-MS, Chapter XVII, Figure 3, page 335
#     )
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.xlim(1)#, max(time)-tStartTransient)
# plt.xlabel("Time [s]")
# plt.ylabel("Thermal Power [MW]")
# # plt.ylim((1e-2, 1e3))
# plt.ylim(0)
# # plt.grid()
# plt.legend(bbox_to_anchor=(0, 1.01, 1, 0.2),
#             loc="lower left",
#             mode='expand',
#             ncol=2)
# # plt.xscale('log')
# # plt.yscale('log')
plt.tight_layout()
plt.savefig("fig_results_transientPower.png")


plt.figure("Temperature", figsize=(5, 4), dpi=200)
# plt.plot(time, Tmin, label="Coolant min")
# plt.plot(time, Tavg, label="Coolant avg")
# plt.plot(time, Tmax, label="Coolant max")
# plt.plot(time, Tfuel, label="Fuel avg")

# if (Tinlet != []):
#     plt.plot(
#         np.arange(tStartTransient, len(Tinlet)+tStartTransient, 1),
#         Tinlet,
#         label="Inlet"
#     )
# else:
#     plt.plot(
#         dataInlet['Time'],
#         dataInlet['areaAverage(T)'],
#         label="Inlet"
#     )

# if (Toutlet != []):
#     plt.plot(
#         np.arange(tStartTransient, len(Toutlet)+tStartTransient, 1),
#         Toutlet,
#         label="Outlet"
#     )
# else:
#     plt.plot(
#         dataOutlet['Time'],
#         dataOutlet['areaAverage(T)'],
#         label="Outlet"
#     )

# plt.plot(time, Tclad, label="Tclad")
# plt.plot(time, Tcool, label="Tcool")
plt.plot( # Experimental data from LANL Report LA-3185-MS, Chapiter III, Figure 10, page 63
    [20822.6459790389, 20825.9590417786, 20828.697520356, 20830.1223414074,
     20830.7821441024, 20831.471993765, 20832.9894169042, 20834.8642246242,
     20836.9568472219, 20840.0407121029, 20842.575117419, 20844.889451771,
     20847.0378612302, 20849.3519904834, 20851.608897184, 20853.8071456406,
     20856.7150357897, 20858.8521648174, 20860.3879442952, 20861.7069343889,
     20863.4145866235, 20865.8337264393, 20868.8002748323, 20871.9284410444,
     20875.3837397707, 20879.0545972886, 20882.4506224747, 20885.5170539615,
     20890.0056402158],
    [1785.2989239822, 1786.74658051972, 1784.82206899522, 1777.64189766803,
     1757.64304389524, 1717.37016196649, 1630.64237485489, 1533.146894489,
     1445.10989261029, 1315.37115299957, 1199.43727980393, 1091.20999849457,
     995.394620107881, 888.727799741575, 800.762582484566, 742.422194938162,
     701.557944459257, 691.567917939987, 673.512376425948, 638.196051709497,
     562.475239924523, 490.18641671965, 419.697336530138, 369.566033221895,
     330.501525758351, 318.060567155484, 313.302272231678, 316.202712779702,
     332.208974260104],
    ls='--',
    label="Station 26 Exp."
)
plt.scatter( # Experimental data from LANL Report LA-3185-MS, Chapiter III, Figure 10, page 63
    [20824.8428855157, 20827.4301540999, 20829.8714356646, 20831.7326704142,
     20833.6024391758, 20834.7599615194, 20835.7675955809, 20836.7805938786,
     20838.7054901404, 20839.6499503476, 20844.8224038871, 20849.7859403816,
     20859.9014487871, 20864.8834054476, 20869.9596352319, 20874.7647494325,
     20879.8839152717, 20884.8536804866, 20889.7928340899],
    [1908.26865670308, 1904.45507095333, 1893.86943962264, 1829.38329071499,
     1706.24515844091, 1629.69759637537, 1573.19032398213, 1479.81609061569,
     1381.84210734856, 1355.47112740814, 1160.14381678252, 986.504982075472,
     758.303596822244, 660.088247368422, 520.020001142007, 425.073439473685,
     393.957889920556, 379.531088133069, 373.46934875869],
    marker='x',
    label="Nozzle Chamber Exp."
)
# plt.hlines([83.3, 1972], xmin=min(time), xmax=max(time), ls='--', color='grey')
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")
plt.ylim(0)
# plt.grid()
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("fig_results_transientTemperature.png")


plt.figure("Reactivity", figsize=(5, 4), dpi=200)
plt.plot(time, totalRho, label="Total")
plt.plot(time, doppler, label="Doppler")
plt.plot(time, axial, label="Axial")
# plt.plot(time, radial, label="Radial")
# plt.plot(time, reactTcool, label="T Coolant")
plt.plot(time, reactRhoCool, label=r"$\rho$ Coolant")
plt.plot(time, extRho, label="Control Drum from Exp.")
plt.xlabel("Time [s]")
plt.ylabel("Reactivity [pcm]")
# plt.grid()
plt.legend(loc='right')
plt.tight_layout()
plt.savefig("fig_results_transientReactivity.png")


plt.figure("Massflow", figsize=(5, 4), dpi=200)
# if (massFlow != []):
#     plt.plot(
#         np.arange(tStartTransient, len(massFlow)+tStartTransient, 1),
#         massFlow,
#         label="GeN-Foam"
#     )
# else:
#     plt.plot(
#         dataInlet["Time"],
#         dataInlet["massFlow"],
#         label="GeN-Foam"
#     )
plt.plot( # Experimental data from LANL Report LA-3185-MS, Chapiter III, Figure 9, page 61
    [20822.9602774092, 20830.4146364074, 20832.1971992561, 20834.1850959171,
    20835.7713343801, 20838.8368569558, 20843.3830541666, 20847.558850965,
    20849.7798629177, 20852.8192262901, 20854.8809654861, 20856.8160802331,
    20858.5171526736, 20860.6469695302, 20865.3325666147, 20868.1548057035,
    20870.7108176851, 20873.7990521272, 20879.230548619, 20882.7979918539,
    20887.430343517, 20889.9328783235, 20926, 20937, 20962.4, 20964, 20991,
    21013],
    [31.082840078101, 31.0423682460298, 30.9092973877877, 29.9607191133102,
    28.7675359187008, 26.332393087713, 22.8033544527425, 19.5231124633684,
    17.8822687573944, 15.4966250794625, 13.9060032378315, 13.0070686586673,
    12.7016902894025, 12.6901269088107, 12.6646874715088, 12.6740432976239,
    12.6848445463131, 12.668077644455, 12.6879456347445, 12.6685769722533,
    12.6434266194661, 12.6298396472708, 12.6298396472708, 6.35029, 6.35029,
    0.907185, 3.85554, 2.58548],
    ls='--',
    label="Experimental"
)
plt.xlabel("Time [s]")
plt.ylabel("Mass flow rate [kg/s]")
plt.ylim(0)
# plt.grid()
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("fig_results_transientMassflowrate.png")


plt.figure("CoolantDensity", figsize=(5, 4), dpi=200)
plt.plot(time, rhoCool, label="Average")
# if (massFlow == []):
#     plt.plot(dataInlet['Time'], dataInlet['rho'], label="Inlet")
#     plt.plot(dataOutlet['Time'], dataOutlet['rho'], label="Outlet")
plt.xlabel(r"Time [$s$]")
plt.ylabel(r"Averaged hydrogen density [$kg/m^3$]")
# plt.ylim(0)
# plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("fig_results_transientHydrogenDensity.png")


#==============================================================================*
