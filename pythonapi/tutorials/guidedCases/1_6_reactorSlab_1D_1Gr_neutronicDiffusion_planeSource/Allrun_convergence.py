"""
Mesh refinement study for 1D slab reactor with plane source boundary condition.
Runs multiple simulations with different mesh resolutions.
"""

#==============================================================================*
# Imports

from Allrun import *

from tabulate import tabulate


#==============================================================================*
# Parameters

nz_values = [21, 51, 101, 201, 501, 1001]  # Different mesh resolutions to test

# Storage for results
results = {
    'nz_values': nz_values,
    'flux0': [],
    'z_mesh': []
}


#==============================================================================*
# Loop over different mesh resolutions

for nz in nz_values:
    print(f"\n{'='*60}")
    print(f"Running simulation with nz={nz}")
    print(f"{'='*60}\n")

    # Recreate mesh with new nz value
    nMesh = createMesh(nz=nz)

    # Update solver with new mesh
    neutronicsSolver.mesh = nMesh

    # Clean and export
    ffn.allclean()
    model.export_to_openfoam()

    # Run
    ffn.run(model, is_preprocessing=True)

    # Extract flux data
    points, flux0 = model.sample_over_line(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux0",
        point1=(0, 0, fuelLength/nz * 0.5),
        point2=(0, 0, fuelLength/nz * (nz-0.5)),
    )

    z = points[:, 2]

    # Store results
    results['flux0'].append(flux0)
    results['z_mesh'].append(z)


#==============================================================================*
# Calculate analytical solution and errors for each mesh

L2_errors = []
Linf_errors = []

for i, nz in enumerate(nz_values):
    z = results['z_mesh'][i]
    flux_sim = results['flux0'][i]

    # Analytical solution at this mesh
    flux_ana = flux_analytical(z, nuSigmaf, removalXS, D, fuelLength, externalSourceCurrent)

    # Normalize both solutions
    flux_sim_norm = flux_sim / np.max(flux_sim)
    flux_ana_norm = flux_ana / np.max(flux_ana)

    # Calculate errors on normalized flux
    diff = flux_sim_norm - flux_ana_norm

    # L2 error (relative)
    L2_error = np.sqrt(np.trapezoid(diff**2, z) / np.trapezoid(flux_ana_norm**2, z))
    L2_errors.append(L2_error)

    # L∞ error (maximum absolute difference)
    Linf_error = np.max(np.abs(diff))
    Linf_errors.append(Linf_error)

# Store errors in results
results['L2_errors'] = L2_errors
results['Linf_errors'] = Linf_errors

# Use most refined mesh for analytical plotting
z_analytical = results['z_mesh'][-1]
flux_analytical_refined = flux_analytical(z_analytical, nuSigmaf, removalXS, D, fuelLength, externalSourceCurrent)
flux_analytical_refined_norm = flux_analytical_refined / np.max(flux_analytical_refined)

print(f"\n{'='*60}")

# Summary plots

# Plot 1: L2 and L∞ Error convergence
fig1, (ax_L2, ax_Linf) = plt.subplots(nrows=2, figsize=(5, 6), dpi=200)

# L2 Error subplot
ref_nz = results['nz_values'][0]
ref_L2 = results['L2_errors'][0]
ax_L2.loglog(
    results['nz_values'],
    results['L2_errors'],
    'bo-',
    label="L2 error",
    markersize=8,
    linewidth=2
)
second_order_L2 = ref_L2 * (ref_nz / np.array(results['nz_values']))**2
first_order_L2 = ref_L2 * (ref_nz / np.array(results['nz_values']))**1
ax_L2.loglog(results['nz_values'], second_order_L2, 'r--', label="2nd order", alpha=0.7)
ax_L2.loglog(results['nz_values'], first_order_L2, 'g--', label="1st order", alpha=0.7)
ax_L2.set_xlabel("Number of axial cells (nz) [-]")
ax_L2.set_ylabel("L2 Error [-]")
ax_L2.set_title("L2 Error Convergence - Norm. Flux ($k_{inf}$ = "+f"{k_inf:.4f})")
ax_L2.legend()
ax_L2.grid(True, alpha=0.3)

# L∞ Error subplot
ref_Linf = results['Linf_errors'][0]
ax_Linf.loglog(
    results['nz_values'],
    results['Linf_errors'],
    'mo-',
    label="L∞ error",
    markersize=8,
    linewidth=2
)
second_order_Linf = ref_Linf * (ref_nz / np.array(results['nz_values']))**2
first_order_Linf = ref_Linf * (ref_nz / np.array(results['nz_values']))**1
ax_Linf.loglog(results['nz_values'], second_order_Linf, 'r--', label="2nd order", alpha=0.7)
ax_Linf.loglog(results['nz_values'], first_order_Linf, 'g--', label="1st order", alpha=0.7)
ax_Linf.set_xlabel("Number of axial cells (nz) [-]")
ax_Linf.set_ylabel("L∞ Error [-]")
ax_Linf.set_title("L∞ Error Convergence - Norm. Flux ($k_{inf}$ = "+f"{k_inf:.4f})")
ax_Linf.legend()
ax_Linf.grid(True, alpha=0.3)

fig1.tight_layout()
fig1.savefig("fig_results_error_convergence.png")
plt.close()


# Plot 2: Relative errors vs axial position
fig2, ax_rel_err = plt.subplots(figsize=(5, 4), dpi=200)

for i, nz in enumerate(nz_values):
    z = results['z_mesh'][i]
    flux_sim = results['flux0'][i]

    # Analytical solution at this mesh
    flux_ana = flux_analytical(z, nuSigmaf, removalXS, D, fuelLength, externalSourceCurrent)

    # Normalize both
    flux_sim_norm = flux_sim / np.max(flux_sim)
    flux_ana_norm = flux_ana / np.max(flux_ana)

    # Relative error
    rel_error = (flux_sim_norm - flux_ana_norm) / flux_ana_norm * 100
    ax_rel_err.plot(z, rel_error, label=f'nz={nz}', alpha=0.7)

ax_rel_err.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
ax_rel_err.set_xlabel("Axial position [m]")
ax_rel_err.set_ylabel("Relative error [%]")
ax_rel_err.set_title("Rel. Error vs Position - Norm. Flux ($k_{inf}$ = "+f"{k_inf:.4f})")
ax_rel_err.legend()
ax_rel_err.grid(True, alpha=0.3)
ax_rel_err.set_xlim(0, fuelLength)

fig2.tight_layout()
fig2.savefig("fig_results_relative_errors.png")
plt.close()

# Plot 3: Flux distribution comparison
fig3, ax_flux = plt.subplots(figsize=(5, 4), dpi=200)

for i, nz in enumerate(nz_values):
    z = results['z_mesh'][i]
    flux0 = results['flux0'][i]
    flux0_norm = flux0 / np.max(flux0)
    ax_flux.plot(z, flux0_norm, 'o-', label=f'nz={nz}', alpha=0.7, markersize=3)

ax_flux.plot(
    z_analytical,
    flux_analytical_refined_norm,
    'k--',
    linewidth=2,
    label='Analytical ($k_{inf}$'+f'={k_inf:.4f})'
)
ax_flux.set_xlabel("Axial position [m]")
ax_flux.set_ylabel("Normalized neutron flux [-]")
ax_flux.set_title("Norm. Flux Dist. Comparison ($k_{inf}$ = "+f"{k_inf:.4f})")
ax_flux.legend()
ax_flux.grid(True, alpha=0.3)
ax_flux.set_xlim(0, fuelLength)

fig3.tight_layout()
fig3.savefig("fig_results_flux_comparison.png")
plt.close()

# Print summary table
print("\nMesh Study Summary:")

table_data = []
for i, nz in enumerate(nz_values):
    table_data.append([
        nz,
        f"{results['L2_errors'][i]:.6e}",
        f"{results['Linf_errors'][i]:.6e}"
    ])

print(tabulate(
    table_data,
    headers=['nz', 'L2 error', 'L∞ error'],
    tablefmt='pipe'
))


#==============================================================================*
