"""
Fuchs Reactivity Insertion Experiment - Convergence Study
Varies solver parameters (maxPowerVariation and maxDeltaT) to assess convergence
"""
#==============================================================================*
# Imports

from Allrun import *

# In this convergence study, we use a shorter t_end (default: 2.5 s after pulse start) to focus on the pulse and immediate aftermath.
# This avoids long runtimes due to small maxDeltaT required for accuracy during the pulse. You can change t_end if needed.
import os
import shutil
from tabulate import tabulate


#==============================================================================*
# Parameters

# Dataset 1: Vary maxPowerVariation (keep maxDeltaT fixed)
maxPowerVariation_values = [0.1, 0.05, 0.01, 0.005, 0.001]  # Power variation threshold

# Storage for results
results_powerVar = {
    'maxPowerVariation': [],
    'Pmax': [],
    'Energy': [],
    'FWHM': [],
    'Pmax_error': [],
    'Energy_error': [],
    'FWHM_error': []
}


#==============================================================================*
# Dataset 1: Vary maxPowerVariation

print("\n" + "="*80)
print("DATASET 1: Varying maxPowerVariation")
print("="*80 + "\n")

steadyStateCaseFolderName = model.caseFolder
transientCaseFolderName = "transient"

# Clean and export
ffn.allclean()
model.export_to_openfoam()

# Preprocessing and run
ffn.run(model, is_preprocessing=True)


for iteration, maxPowerVar in enumerate(maxPowerVariation_values):
    print(f"\nRunning with maxPowerVariation={maxPowerVar:.1e} ({iteration+1}/{len(maxPowerVariation_values)})")

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
    settings.maxPowerVariation = maxPowerVar
    settings.maxDeltaT = 3e-3

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
        promptGenerationTime, beta_total, reactivityInsertion, gamma, initialPowerScaled,
        wedgeAngle
    )

    # Store results
    results_powerVar['maxPowerVariation'].append(maxPowerVar)
    results_powerVar['Pmax'].append(metrics['Pmax_sim'])
    results_powerVar['Energy'].append(metrics['E_sim'])
    results_powerVar['FWHM'].append(metrics['FWHM_sim'])
    results_powerVar['Pmax_error'].append(metrics['Pmax_rel_error'])
    results_powerVar['Energy_error'].append(metrics['E_rel_error'])
    results_powerVar['FWHM_error'].append(metrics['FWHM_rel_error'])

    print(f"FWHM   error: {metrics['FWHM_rel_error']:.2f}%")
    print(f"Energy error: {metrics['E_rel_error']:.2f}%")
    print(f"Pmax   error: {metrics['Pmax_rel_error']:.2f}%")



    # FWHM convergence
    #-----------------

    for parameterName, ylabel, ext in [
        ("FWHM_error", "FWHM Relative Error [%]", "FWHM"),
        ("Energy_error", "Energy Relative Error [%]", "Energy"),
        ("Pmax_error", "Pmax Relative Error [%]", "Pmax"),
    ]:
        fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
        ax.semilogx(
            results_powerVar['maxPowerVariation'],
            results_powerVar[parameterName],
            marker="o",
            color='black'
        )
        ax.axhline(0, color='tab:gray')
        ax.set_xlabel("maxPowerVariation")
        ax.set_ylabel(ylabel)
        fig.tight_layout()
        fig.savefig(f"fig_convergence_maxPowerVariation_{ext}.png")
        plt.close()


#==============================================================================*
# Print summary table of all runs

print("\nConvergence Study Summary Table:\n")

# Dataset 1 table
table1 = []
for i, _ in enumerate(results_powerVar['maxPowerVariation']):
    table1.append([
        f"{results_powerVar['maxPowerVariation'][i]:.3e}",
        f"{results_powerVar['Pmax_error'][i]:.2f}",
        f"{results_powerVar['Energy_error'][i]:.2f}",
        f"{results_powerVar['FWHM_error'][i]:.2f}"
    ])

print("Dataset 1: Varying maxPowerVariation (maxDeltaT fixed)")
print(tabulate(
    table1,
    headers=["maxPowerVariation", "Pmax err. [%]", "Energy err. [%]", "FWHM err. [%]"],
    tablefmt="pipe"
))


#==============================================================================*
