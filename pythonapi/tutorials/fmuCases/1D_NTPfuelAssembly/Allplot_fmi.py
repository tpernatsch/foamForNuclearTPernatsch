#==============================================================================*
# Imports

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate
import foamForNuclear as ffn
import sys


#==============================================================================*
# General inputs

dpi: int = 200
figsize = (5, 4)

prefix = "model.root.system1." if False else ""



#==============================================================================*
# Functions

def doublingTime(time: list[float], power: list[float]):
    doubleTime = []
    for ti, tf, power_i, power_f in zip(time[:-1], time[1:], power[:-1], power[1:]):
        dt = tf-ti
        dlogPower = np.log(power_f) - np.log(power_i)
        dlogPdt = dlogPower/dt
        if (dlogPdt > 0):
            doubleTime.append(np.log(2)/dlogPdt)
        else:
            doubleTime.append(float('nan'))

    return(time[1:], doubleTime)

def doublingTime2(time: list[float], power: list[float]):
    doubleTime = []
    for ti, tf, power_i, power_f in zip(time[:-1], time[1:], power[:-1], power[1:]):
        dt = tf-ti
        dPower = power_f - power_i
        powerMean = (power_f+power_i)/2
        dPdt = dPower/dt
        if (dPdt > 0):
            doubleTime.append(powerMean*np.log(2)/dPdt)
        else:
            doubleTime.append(float('nan'))

    return(time[1:], doubleTime)


def intersection(x, xList, yList):
    for xi, xf, yi, yf in zip(xList[:-1], xList[1:], yList[:-1], yList[1:]):
        if (xi <= x and x <= xf):
            m = (yf-yi) / (xf-xi)
            return(m * (x-xi) + yi)

    return(None)


def plotPressure(time, data, caseFolderName: str, restartPoints: list[float]=None):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    for variable in ['pTurbine_in', 'SinkP1.in_p0', 'pNozzleChamber_out']:
        ax.plot(time, [val*1e-5 for val in data[prefix+variable]], label=variable)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_xlim(0)
    ax.set_ylim(29)
    if (restartPoints is not None):
        ax.set_xticks(restartPoints, minor=True)
        ax.xaxis.grid(True, which='minor')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_pressure.png")
    plt.close()


def plotTemperature(time, data, caseFolderName: str, restartPoints: list[float]=None):
    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    ax.plot(time, data[prefix+'sensTnozzleWall_in.T'], label='Nozzle Wall inlet', color='tab:orange', ls='-')
    ax.plot(time, data[prefix+'sensTnozzleWall_out.T'], label='Nozzle Wall outlet', color='tab:orange', ls='--')
    ax.plot(time, data[prefix+'Tturbine_out'], label='Tturbine_out', color='tab:green')
    ax.plot(time, data[prefix+'Tnozzle_in'], label='Tnozzle_in', color='tab:purple')
    # ax.plot(time, data[prefix+'sourceMassFlow.in_T'], label='Tnozzle srcMFlow', color='tab:purple', ls='--')
    ax.plot(time, data[prefix+'NozzleChamber.T[1]'], label='Nozzle Chamber inlet', color='tab:red', ls='-')
    ax.plot(time, data[prefix+'NozzleChamber.T[11]'], label='Nozzle Chamber outlet', color='tab:red', ls='--')
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [K]")
    ax.set_xlim(0)
    ax.set_ylim(100)
    if (restartPoints is not None):
        ax.set_xticks(restartPoints, minor=True)
        ax.xaxis.grid(True, which='minor')
    ax.legend(
        bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
        ncols=2, mode="expand", borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_temperature.png")
    plt.close()


def plotMassFlowRate(time, data, caseFolderName: str, restartPoints: list[float]=None):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.plot(time, data[prefix+'sensW.w'], label='Turbine outlet')
    ax.plot(time, [e*1500 for e in data[prefix+'mFlowTurbine_out']], label='Turbine outlet damped')
    ax.plot(time, data[prefix+'mFlowCoreOutlet_in'], label='Core outlet')
    ax.plot(time, data[prefix+'controlSystem.PI_massFlow.u_s'], label='Cmd', alpha=0.8, color="tab:grey")
    ax.plot(time, data[prefix+'controlSystem.PI_massFlow.u_m'], label='Measure', alpha=0.8, color="tab:grey", ls='--')
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mass flow rate [kg/s]")
    ax.set_xlim(0)
    # ax.set_ylim((90, 110))
    # ax.set_ylim((25, 35))
    ax.set_ylim((0, 35))
    if (restartPoints is not None):
        ax.set_xticks(restartPoints, minor=True)
        ax.xaxis.grid(True, which='minor')
    ax.legend(
        bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
        ncols=2, mode="expand", borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_massflow.png")
    plt.close()


def plotReactivity(time, data, caseFolderName: str, restartPoints: list[float]=None):
    model = ffn.case.Case(caseFolder=caseFolder)
    model.settings.application = "GeN-Foam"
    dataPK = model.get_parameters_from_point_kinetics(writeInterval=1)
    dataPK['time'] = [t-tmin for t in dataPK['time']]
    print(f"End read, N points = {len(dataPK['time'])}")

    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    ax.plot(time, [e*1e5 for e in data[prefix+'controlSystem.extReactivity']], label='Control drum')
    ax.plot(dataPK['time'], dataPK['totalReactivity'], label='Total')
    ax.plot(dataPK['time'], dataPK['dopplerReactivity'], label='Doppler')
    ax.plot(dataPK['time'], dataPK['TFuelReactivity'], label='Fuel temperature')
    ax.plot(dataPK['time'], dataPK['rhoCoolReactivity'], label='Coolant density')
    # ax.hlines(y=[xmin], xmin=0, xmax=max(time), color='tab:grey', alpha=0.5)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Reactivity [pcm]")
    ax.set_xlim(0)
    ax.set_yticks([0.1], minor=True)
    ax.yaxis.grid(True, which='minor')
    if (restartPoints is not None):
        ax.set_xticks(restartPoints, minor=True)
        ax.xaxis.grid(True, which='minor')
    ax.legend(
        bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
        ncols=2, mode="expand", borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_reactivity.png")
    plt.close()


    # Time step
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    deltaT = [tf-ti for ti, tf in zip(dataPK['time'][:-1], dataPK['time'][1:])]
    ax.plot(dataPK['time'][:-1], deltaT)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Time step [s]")
    ax.set_xlim(0)
    ax.set_yscale('log')
    # ax.legend()
    if (restartPoints is not None):
        ax.set_xticks(restartPoints, minor=True)
        ax.xaxis.grid(True, which='minor')
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_timeStep.png")
    plt.close()


def plotPower(time, data, caseFolderName: str, restartPoints: list[float]=None):
    power = [e*1e-6 for e in data[prefix+'controlSystem.PI_Power.u_m']]

    fig, axPower = plt.subplots(figsize=figsize, dpi=dpi)
    axControlDrum = axPower.twinx()
    pPowerCmd, = axPower.plot(
        time,
        [e*1e-6 for e in data[prefix+'controlSystem.PI_Power.u_s']],
        label='Cmd',
        ls='--',
        color='tab:blue'
    )
    pPowerMeasure, = axPower.plot(
        time,
        power,
        label='Measure',
        color='tab:blue'
    )
    pControlDrum, = axControlDrum.plot(
        time,
        data[prefix+'controlSystem.controlDrumReactivity.angle'],
        color='tab:grey',
        label="Control drum",
        ls='-.'
    )
    pControlDrum2, = axControlDrum.plot(
        time,
        data[prefix+'controlSystem.PI_Power.y'],
        color='tab:orange',
        label="PI Power CD Cmd",
        ls='-.'
    )
    pControlDrum3, = axControlDrum.plot(
        time,
        data[prefix+'controlSystem.PI_DoublingTime.y'],
        color='tab:red',
        label="PI Doubling Time CD Cmd",
        ls='-.'
    )
    axPower.set_xlabel("Time [s]")
    axPower.set_ylabel("Power [MW]")
    axPower.set_xlim(0)
    axPower.set_ylim(0)
    axControlDrum.set_ylabel("Control drum orientation [deg]")
    lines = [pPowerCmd, pPowerMeasure, pControlDrum, pControlDrum2, pControlDrum3]
    axPower.legend(
        handles=lines,
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncols=2,
        mode="expand",
        borderaxespad=0.
    )
    if (restartPoints is not None):
        axPower.set_xticks(restartPoints, minor=True)
        axPower.xaxis.grid(True, which='minor')
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_power.png")
    plt.close()


def plotDoublingTime(time, data, caseFolderName: str, restartPoints: list[float]=None):
    promptGenerationTime = 3.571485e-05

    Beta = [0.000228513, 0.00118003, 0.00112685, 0.00252764, 0.0010375, 0.000434565]

    doublingTimeCompute = data[prefix+'controlSystem.doublingTimeCalculator.doublingTime']
    power = [e*1e-6 for e in data[prefix+'controlSystem.PI_Power.u_m']]

    fig, axDoublingTime = plt.subplots(figsize=figsize, dpi=dpi)
    axPower = axDoublingTime.twinx()
    pDoublingTime, = axDoublingTime.plot(
        time,
        [val if val > 0 else float('nan') for val in doublingTimeCompute],
        # doublingTimeCompute,
        label='Positive $t_d$',
        # color='tab:grey',
        color='black',
        # ls='-.',
        # linewidth=1

    )
    # pDoublingTimeLimited, = axDoublingTime.plot(
    #     time,
    #     data[prefix+'controlSystem.PI_DoublingTime.u_m'],
    #     label='Doubling time limited',
    #     color='black',
    # )
    pPowerMeasure, = axPower.plot(
        time,
        power,
        label='Power Measure',
        # color='tab:blue',
        color='tab:grey',
        ls='--'
    )
    axDoublingTime.set_xlabel("Time [s]")
    axDoublingTime.set_ylabel("Doubling time [s]")
    axPower.set_ylabel("Power [MW]")
    axDoublingTime.set_xlim(0)
    axDoublingTime.set_ylim(1)
    axDoublingTime.set_yscale('log')
    axDoublingTime.set_yticks([20], minor=True)
    axDoublingTime.yaxis.grid(True, which='minor')
    axPower.set_ylim(0)
    if (restartPoints is not None):
        axDoublingTime.set_xticks(restartPoints, minor=True)
        axDoublingTime.xaxis.grid(True, which='minor')
    # lines = [pDoublingTime, pDoublingTimeLimited, pPowerMeasure]
    lines = [pDoublingTime, pPowerMeasure]
    axDoublingTime.legend(
        handles=lines,
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncols=2,
        mode="expand",
        borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_doublingTime.png")
    plt.close()


def plotRotationSpeedAndValveOpening(time, data, caseFolderName: str, restartPoints: list[float]=None):
    fig, axRotation = plt.subplots(figsize=figsize, dpi=dpi)
    axValve = axRotation.twinx()
    pRotation, = axRotation.plot(
        time,
        [e * 60/(2*np.pi) for e in data['Turbine1.omega']],
        label="Shaft speed",
        color="black"
    )
    # pRotation, = axRotation.plot(time, data[prefix+'Turbine1.omega'])
    pValve, = axValve.plot(
        time,
        [e*100 for e in data[prefix+'valveLin.cmd']],
        label="Tank valve",
        color="tab:grey",
        ls='--'
    )
    axRotation.set_xlabel("Time [s]")
    axRotation.set_ylabel("Arbor rotation speed [rpm]")
    # axRotation.set_ylabel("Arbor rotation speed [rad/s]")
    axValve.set_ylabel("Valve opening [%]")
    axRotation.set_xlim(0)
    if (restartPoints is not None):
        axRotation.set_xticks(restartPoints, minor=True)
        axRotation.xaxis.grid(True, which='minor')
    lines = [pRotation, pValve]
    axRotation.legend(
        handles=lines,
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncols=2,
        mode="expand",
        borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_rotationAndValve.png")
    plt.close()


def plotThrust(time, data, caseFolderName: str, restartPoints: list[float]=None):
    fig, axThrust = plt.subplots(figsize=figsize, dpi=dpi)
    axV = axThrust.twinx()
    pF, = axThrust.plot(time, [e*1e-3 for e in data[prefix+'controlSystem.thrustCalculator.thrust']], label="Thrust", color='tab:blue')
    pFcmd, = axThrust.plot(time, [e*1e-3 for e in data[prefix+'controlSystem.PI_Thrust.u_s']], label="Thrust Cmd", color='tab:blue', ls='--')
    pV, = axV.plot(
        time, data[prefix+'controlSystem.thrustCalculator.exhaustVelocity'],
        color='tab:orange',
        label='$v_{out}$'
    )
    axThrust.set_xlabel("Time [s]")
    axThrust.set_ylabel("Thrust [kN]")
    axV.set_ylabel("Exhaust velocity [m/s]")
    axThrust.set_xlim(0)
    if (restartPoints is not None):
        axThrust.set_xticks(restartPoints, minor=True)
        axThrust.xaxis.grid(True, which='minor')
    axThrust.legend(
        handles=[pF, pFcmd, pV],
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncols=2,
        mode="expand",
        borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_rocketPerformance.png")
    plt.close()


def plotSummaries(time, data, caseFolderName: str):
    xmin = 0 # min(time)
    xmax = int(max(time)/10+1)*10

    moonTransfer = 3119e-3 # km/s

    deltaV = [e/1e3 for e in data[prefix+'controlSystem.deltaVCalculator.deltaV']]
    remainMass = [e/1e3 for e in data[prefix+'controlSystem.deltaVCalculator.totalMassPropellant.y']]
    initialMass = [e/1e3 for e in data[prefix+'controlSystem.deltaVCalculator.totalMassPropellant.u1']]

    power = [e*1e-6 for e in data[prefix+'controlSystem.PI_Power.u_m']]
    doublingTimeCompute = list(data[prefix+'controlSystem.doublingTimeCalculator.doublingTime'])
    exhaustTemperature = list(data[prefix+'Tnozzle_in'])
    thrust = [e*1e-3 for e in data[prefix+'controlSystem.thrustCalculator.thrust']]

    lastDeltaV = list(deltaV)[-1]

    vStartDeltaV = None # intersection(tStartBurn, time, deltaV)
    if (vStartDeltaV is None):
        vStartDeltaV = 0
    vEndDeltaV = vStartDeltaV + moonTransfer

    print(tabulate(
        [
            ["Delta-v progress", f"{100*lastDeltaV/moonTransfer:.1f}", "%"],
            ["Remaining propellant", f"{remainMass[-1]:.3f} / {initialMass[-1]:.3f} ({remainMass[-1]/initialMass[-1]*100:.1f}%)", "t"],
            ["Power", f"{power[-1]:.4g}", "MW"],
            ["Exhaust temperature", f"{exhaustTemperature[-1]:.5g}", "K"],
            ["Thrust", f"{thrust[-1]:.4g}", "kN"],
            ["Doubling time", f"{doublingTimeCompute[-1]:.4g}", "s"],
        ],
        ["Parameter", "Value", "Unit"],
        tablefmt="simple_grid"
    ))


    # Delta v
    fig, axDeltaV = plt.subplots(figsize=(6, 4), dpi=dpi)
    axPropellantMass = axDeltaV.twinx()
    axThrust = axDeltaV.twinx()
    axThrust.spines.right.set_position(("axes", 1.2))
    deltaVLine, = axDeltaV.plot(
        time,
        deltaV,
        color="black",
        label="$\Delta v$"
    )
    remainMassLine, = axPropellantMass.plot(
        time,
        remainMass,
        label="Remaining H$_2$",
        ls='-.'
    )
    initMassLine, = axPropellantMass.plot(
        time,
        initialMass,
        label="Initial H$_2$",
        ls='--',
        color="tab:grey"
    )
    thrCmd, = axThrust.plot(time, [e*1e-3 for e in data[prefix+'controlSystem.PI_Thrust.u_s']], label="Thrust Cmd", color='tab:green', ls='--')
    thrEst, = axThrust.plot(time, thrust, label="Thrust Est.", color='tab:green')
    axDeltaV.set_xlabel("Time [s]")
    axDeltaV.set_ylabel("$\Delta v$ [km/s]")
    axDeltaV.set_xlim((xmin, xmax))
    axDeltaV.set_ylim((0, 4.0))
    axPropellantMass.set_ylabel("Propellant mass [t]")
    axPropellantMass.set_ylim((0, 35))
    axThrust.set_ylim((0, 270))
    axThrust.set_ylabel("Thrust [kN]")
    axDeltaV.set_yticks([vStartDeltaV, vEndDeltaV], minor=True)
    axDeltaV.yaxis.grid(True, which='minor')
    axDeltaV.text(
        x=(0+0.05*xmax)/2, y=(vStartDeltaV+vEndDeltaV)/2, s=f"Moon transfer $\Delta v$ = {moonTransfer:.2f} km/s",
        ha="center", va="center", rotation='vertical', fontsize=8
    )
    axDeltaV.fill_between(
        x=[0, 0.05*xmax],
        y1=[vStartDeltaV, vStartDeltaV],
        y2=[vEndDeltaV, vEndDeltaV],
        color="tab:grey",
        alpha=0.3
    )
    lines = [deltaVLine, remainMassLine, initMassLine, thrCmd, thrEst]
    axDeltaV.legend(
        handles=lines,
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncols=2,
        mode="expand",
        borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_summaryDeltaV.png")
    plt.close()


    # Summary
    fig, axPower = plt.subplots(figsize=(6, 4), dpi=dpi)
    axTemperature = axPower.twinx()
    axThrust = axPower.twinx()
    axThrust.spines.right.set_position(("axes", 1.2))
    powCmd, = axPower.plot(time, [e*1e-6 for e in data[prefix+'controlSystem.PI_Power.u_s']], label='Power Cmd', color='tab:blue', ls='--')
    powMea, = axPower.plot(time, power, label='Power Measure', color='tab:blue')
    Tnozzle, = axTemperature.plot(time, exhaustTemperature, label='Exhaust temperature', color='tab:red')
    thrCmd, = axThrust.plot(time, [e*1e-3 for e in data[prefix+'controlSystem.PI_Thrust.u_s']], label="Thrust Cmd", color='tab:green', ls='--')
    thrEst, = axThrust.plot(time, thrust, label="Thrust Est.", color='tab:green')
    axPower.set_xlim((xmin, xmax))
    axPower.set_ylim((0, 1000))
    axPower.set_xlabel("Time [s]")
    axPower.set_ylabel("Power [MW]")
    axTemperature.set_ylim((0, 3000))
    axTemperature.set_ylabel("Temperature [K]")
    axThrust.set_ylim((0, 270))
    axThrust.set_ylabel("Thrust [kN]")
    lines = [powCmd, powMea, Tnozzle, thrCmd, thrEst]
    axPower.legend(
        handles=lines,
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncols=2,
        mode="expand",
        borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_summaryRocketPerformance.png")
    plt.close()


    # Summary
    fig, axPower = plt.subplots(figsize=(6, 4), dpi=dpi)
    axDoublingTime = axPower.twinx()
    axControlDrum = axPower.twinx()
    axControlDrum.spines.right.set_position(("axes", 1.2))
    powCmd, = axPower.plot(time, [e*1e-6 for e in data[prefix+'controlSystem.PI_Power.u_s']], label='Power Cmd', color='tab:blue', ls='--')
    powMea, = axPower.plot(time, power, label='Power Measure', color='tab:blue')

    pDoublingTime, = axDoublingTime.plot(
        time,
        [val if val > 0 else float('nan') for val in doublingTimeCompute],
        label='Positive $t_d$',
        color='black',
        ls=":"
    )
    pControlDrum, = axControlDrum.plot(
        time,
        data[prefix+'controlSystem.controlDrumReactivity.angle'],
        color='tab:grey',
        label="Control drum",
        ls='-.',
        # marker=".",
        # markevery=200
    )
    pControlDrum2, = axControlDrum.plot(
        time,
        [d+99.9 for d in data[prefix+'controlSystem.PI_Power.y']],
        color='tab:green',
        label="PI Power CD Cmd",
        ls='',
        marker="o",
        markersize=3,
        markevery=100
    )
    pControlDrum3, = axControlDrum.plot(
        time,
        [d+99.9 for d in data[prefix+'controlSystem.PI_DoublingTime.y']],
        color='tab:green',
        label="PI $t_d$ CD Cmd",
        ls='',
        marker="v",
        markersize=3,
        markevery=100
    )
    axPower.set_xlabel("Time [s]")
    axPower.set_ylabel("Power [MW]")
    axPower.set_xlim((xmin, xmax))
    axPower.set_ylim((0, 1000))
    axDoublingTime.set_ylabel("Doubling time [s]")
    axDoublingTime.set_ylim((1, 1e6))
    axDoublingTime.set_yscale('log')
    axDoublingTime.set_yticks([20], minor=True)
    axDoublingTime.yaxis.grid(True, which='minor')
    axControlDrum.set_ylabel("Control drum orientation [deg]")
    axControlDrum.set_ylim((0, 180))
    lines = [powCmd, powMea, pDoublingTime, pControlDrum, pControlDrum2, pControlDrum3]
    axPower.legend(
        handles=lines,
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncols=2,
        mode="expand",
        borderaxespad=0.
    )
    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_{caseFolderName}_summaryPower.png")
    plt.close()


def extractRestartPoints(caseFolder):
    restartPoints = []
    with open(f"{caseFolder}/log.GeN-Foam", 'r') as f:
        for line in f:
            if ("Restart" in line):
                restartPoints.append(float(line.split()[-1]))
                continue

    print(f"Number of FMU restart = {len(restartPoints)}")

    return(restartPoints)


def plotCompare(
        caseFolders,
        parameterInFile: str,
        parameterName: str,
        unit: str='-',
        scale: float=1,
        offset: float=0,
        ylim: tuple[float]=None
    ):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    linestyles = ['-', '--', '-.', ':']

    for caseFolder, ls in zip(caseFolders, linestyles):
        data = pd.read_csv(f'{caseFolder}/Turbomachine.csv', sep=',')

        time = data['time']

        tmin = min(time)
        time = [t-tmin for t in data['time']]

        xmin = 0 # min(time)
        xmax = int(max(time)/10+1)*10

        ax.plot(
            time,
            [e*scale+offset for e in data[parameterInFile]],
            label=caseFolder.replace('.', '').replace('/', ''),
            ls=ls,
            color="black"
        )

    ax.set_xlabel('Time [s]')
    ax.set_ylabel(f'{parameterName.capitalize()} [{unit}]')
    ax.set_xlim((xmin, xmax))
    if (ylim is not None):
        ax.set_ylim(ylim)
    ax.legend()

    fig.tight_layout()
    fig.savefig(f"fig_results_NTP_compare_{parameterName}.png")
    plt.close()


#==============================================================================*
# Main

if __name__ == "__main__":

    if (len(sys.argv[1:]) < 1):
        print(f"Usage:\n    python3 {sys.argv[0]} <path/to/folder1> <path/to/folder2>\n")
        sys.exit(0)

    caseFolders = sys.argv[1:]

    isAddRestartPoints: bool = False
    isPlotReactivity: bool = False


    for caseFolder in caseFolders:
        # Data extraction
        data = pd.read_csv(f'{caseFolder}/Turbomachine.csv', sep=',')

        caseFolderName = caseFolder.replace('/', '').replace('.', '')

        time = data['time']
        print(f"Last time in .CSV: {list(time)[-1]:.6f} s")

        tmin = min(time)
        time = [t-tmin for t in data['time']]

        print(f"N points = {len(time)}")

        # Extract restart of the FMU
        restartPoints = extractRestartPoints(caseFolder) if isAddRestartPoints else None

        # Plot
        plotPressure(time, data, caseFolderName=caseFolderName, restartPoints=restartPoints)
        plotTemperature(time, data, caseFolderName=caseFolderName, restartPoints=restartPoints)
        plotMassFlowRate(time, data, caseFolderName=caseFolderName, restartPoints=restartPoints)
        plotPower(time, data, caseFolderName=caseFolderName, restartPoints=restartPoints)
        plotDoublingTime(time, data, caseFolderName=caseFolderName, restartPoints=restartPoints)
        plotRotationSpeedAndValveOpening(time, data, caseFolderName=caseFolderName, restartPoints=restartPoints)
        plotThrust(time, data, caseFolderName=caseFolderName, restartPoints=restartPoints)

        plotSummaries(time, data, caseFolderName=caseFolderName)

        if (isPlotReactivity):
            plotReactivity(time, data, caseFolderName=caseFolderName, restartPoints=restartPoints)


    if (len(caseFolders) >= 2):
        plotCompare(
            caseFolders=caseFolders,
            parameterInFile='controlSystem.PI_Power.u_m',
            parameterName='power',
            unit='MW',
            scale=1e-6,
            ylim=(0, 600)
        )
        plotCompare(
            caseFolders=caseFolders,
            parameterInFile='controlSystem.thrustCalculator.thrust',
            parameterName='thrust',
            unit='kN',
            scale=1e-3,
            ylim=(0, 250)
        )
        plotCompare(
            caseFolders=caseFolders,
            parameterInFile='controlSystem.deltaVCalculator.deltaV',
            parameterName='deltaV',
            unit='km/s',
            scale=1e-3,
            ylim=(0, 4)
        )
        plotCompare(
            caseFolders=caseFolders,
            parameterInFile='Tnozzle_in',
            parameterName='Tnozzle',
            unit='K',
            ylim=(0, 2500)
        )

#==============================================================================*
