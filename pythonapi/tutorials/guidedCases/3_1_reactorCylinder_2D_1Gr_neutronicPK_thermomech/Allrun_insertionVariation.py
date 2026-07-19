"""
Fuchs Reactivity Insertion Experiment
Varies reactivity insertion
"""
#==============================================================================*
# Imports

from Allrun import *

import os
import shutil
from tabulate import tabulate


#==============================================================================*
# Parameters

reactivityInsertions = [1.5, 2, 2.5, 3] # in $

# Storage for results
results = {
    'reactivityInsertion': [],
    'times': [],
    'powers': [],
    'Pmax': [],
    'Energy': [],
    'FWHM': [],
    'Pmax_error': [],
    'Energy_error': [],
    'FWHM_error': []
}


#==============================================================================*
# Dataset

steadyStateCaseFolderName = model.caseFolder
transientCaseFolderName = "transient"

# Clean and export steady-state
ffn.allclean()
model.export_to_openfoam()

# Preprocessing and run
ffn.run(model, is_preprocessing=True)


for iteration, reactivityInsertion in enumerate(reactivityInsertions):
    print(f"\nRunning with reactivity insertion = {reactivityInsertion:g} $ ({iteration+1}/{len(reactivityInsertions)})")

    reactivityInsertion_pcm = reactivityInsertion * beta_total

    # --- Restart for transient
    if os.path.exists(transientCaseFolderName):
        shutil.rmtree(transientCaseFolderName)

    ffn.duplicateFolder(steadyStateCaseFolderName, transientCaseFolderName)
    model.caseFolder = transientCaseFolderName

    # Reset all settings to Allrun.py defaults for transient, then set sweep
    # param
    settings.endTime = tStartReactitivity + 3
    settings.writeControl = 'adjustableRunTime'
    settings.writeInterval = settings.endTime
    settings.maxPowerVariation = 0.01
    settings.maxDeltaT = 3e-3

    (lastTime, _) = pointKineticsData.externalReactivityTimeProfile.table[-1]
    pointKineticsData.externalReactivityTimeProfile.table[-1] = (lastTime, reactivityInsertion_pcm)


    neutronicsSolver.nuclearData = pointKineticsData
    neutronicsSolver.solver = "pointKinetics"
    neutronicsSolver.eigenvalueNeutronics = False
    neutronicsSolver.fastNeutrons = True


    model.export_to_openfoam()

    ffn.run(model)

    # Post-processing
    res = model.get_parameters_from_point_kinetics()

    # Extract transient data
    time = res['time']
    power = res['totalPower']
    TFuel = res["TFuel"]

    powerScaled = [p * 1e-6 * scalingPower for p in power]

    # Calculate metrics
    metrics = calculate_metrics(
        time, powerScaled,
        promptGenerationTime, beta_total, reactivityInsertion_pcm, gamma, initialPowerScaled,
        wedgeAngle
    )

    # Store results
    results['reactivityInsertion'].append(reactivityInsertion)
    results['times'].append(time)
    results['powers'].append(powerScaled)
    results['Pmax'].append(metrics['Pmax_sim'])
    results['Energy'].append(metrics['E_sim'])
    results['FWHM'].append(metrics['FWHM_sim'])
    results['Pmax_error'].append(metrics['Pmax_rel_error'])
    results['Energy_error'].append(metrics['E_rel_error'])
    results['FWHM_error'].append(metrics['FWHM_rel_error'])

    print(f"FWHM   error: {metrics['FWHM_rel_error']:.2f}%")
    print(f"Energy error: {metrics['E_rel_error']:.2f}%")
    print(f"Pmax   error: {metrics['Pmax_rel_error']:.2f}%")

    # Convergence
    #------------

    for parameterName, ylabel, ext in [
        ("FWHM", "FWHM [s]", "FWHM"),
        ("FWHM_error", "FWHM Relative Error [%]", "relErr_FWHM"),
        ("Energy", "Energy [MJ]", "Energy"),
        ("Energy_error", "Energy Relative Error [%]", "relErr_Energy"),
        ("Pmax", "Pmax [MW]", "Pmax"),
        ("Pmax_error", "Pmax Relative Error [%]", "relErr_Pmax"),
    ]:
        fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
        ax.plot(
            results['reactivityInsertion'],
            results[parameterName],
            marker="o",
            color='black'
        )
        ax.set_xlabel("Reactivity insertion [$]")
        ax.set_ylabel(ylabel)
        fig.tight_layout()
        fig.savefig(f"fig_results_reactivityInsertion_{ext}.png")
        plt.close()


    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    for reactivityInsertion_i, times_i, powers_i in zip(
        results['reactivityInsertion'],
        results['times'],
        results['powers']
    ):
        ax.plot(
            [t - tStartReactitivity for t in times_i],
            powers_i,
            label=f"{reactivityInsertion_i:.1f} $"
        )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Power [MW]")
    ax.set_xlim((0, tEndPlot - tStartReactitivity))
    ax.set_ylim(0)
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"fig_results_reactivityInsertion_power.png")
    plt.close()


#==============================================================================*
# Print summary table of all runs

print("\nStudy Summary Table:\n")

print(tabulate(
    [
        [
            f"{results['reactivityInsertion'][i]:.3e}",
            f"{results['Pmax_error'][i]:.2f}",
            f"{results['Energy_error'][i]:.2f}",
            f"{results['FWHM_error'][i]:.2f}"
        ]
        for i, _ in enumerate(results['reactivityInsertion'])
    ],
    headers=["Reactivity insertion [$]", "Pmax err. [%]", "Energy err. [%]", "FWHM err. [%]"],
    tablefmt="pipe"
))

print()

print(tabulate(
    [
        [
            f"{results['reactivityInsertion'][i]:.3e}",
            f"{results['Pmax'][i]:.2f}",
            f"{results['Energy'][i]:.2f}",
            f"{results['FWHM'][i]:.2f}"
        ]
        for i, _ in enumerate(results['reactivityInsertion'])
    ],
    headers=["Reactivity insertion [$]", "Pmax [MW]", "Energy [MJ]", "FWHM [s]"],
    tablefmt="pipe"
))


#==============================================================================*
