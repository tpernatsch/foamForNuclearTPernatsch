"""
Mesh refinement study for 2-group 1D slab reactor.
Runs multiple simulations with different mesh resolutions.
"""

#==============================================================================*
# Imports

from Allrun import *

from tabulate import tabulate


#==============================================================================*
# Parameters

nz_values = [11, 51, 101, 201, 501, 1001, 10001]  # Different mesh resolutions to test

# Storage for results
results = {
    'nz_values': nz_values,
    'keff': [],
    'flux0': [],
    'flux1': [],
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
    ffn.run(model=model, is_preprocessing=True)

    # Extract results
    keff = model.keff()

    # Extract flux data
    points, flux0 = model.sample_over_line(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux0",
        point1=(0, 0, fuelLength/nz * 0.5),
        point2=(0, 0, fuelLength/nz * (nz-0.5)),
    )

    points, flux1 = model.sample_over_line(
        region=nMesh,
        time=settings.endTime,
        fieldName="flux1",
        point1=(0, 0, fuelLength/nz * 0.5),
        point2=(0, 0, fuelLength/nz * (nz-0.5)),
    )

    z = points[:, 2]

    # Store results
    results['keff'].append(keff)
    results['flux0'].append(flux0)
    results['flux1'].append(flux1)
    results['z_mesh'].append(z)


#==============================================================================*
# Calculate analytical solution and errors for each mesh

# Geometric buckling
B_squared = (np.pi / fuelLength)**2

# Calculate analytical k-effective
denom_fast = removal_xs[0] + diffusion_coefficient[0] * B_squared
denom_thermal = removal_xs[1] + diffusion_coefficient[1] * B_squared
k_analytical = (nu_sigma_f[0] / denom_fast +
                (scattering_matrix[0][1] / denom_fast) * (nu_sigma_f[1] / denom_thermal))

# Flux ratio φ₀/φ₁ (fast/thermal) from thermal group balance
flux_ratio_fast_to_thermal = (diffusion_coefficient[1] * B_squared + removal_xs[1]) / scattering_matrix[0][1]

# Amplitude calculation
integral_sin = 2 * fuelLength / np.pi
power_factor = kappa_J * (sigma_fission[0] * flux_ratio_fast_to_thermal + sigma_fission[1])
amplitude_1 = target_power / (cross_section_area * power_factor * integral_sin)
amplitude_0 = amplitude_1 * flux_ratio_fast_to_thermal

# Calculate errors for each mesh
L2_errors_0 = []
Linf_errors_0 = []
L2_errors_1 = []
Linf_errors_1 = []

for i, nz in enumerate(nz_values):
    z = results['z_mesh'][i]
    flux_sim_0 = results['flux0'][i]
    flux_sim_1 = results['flux1'][i]

    # Analytical solution at mesh points
    flux_analytical_shape = np.sin(z * np.pi / fuelLength)
    flux_ana_0 = amplitude_0 * flux_analytical_shape
    flux_ana_1 = amplitude_1 * flux_analytical_shape

    # Calculate errors over full domain (using absolute flux values)
    diff_0 = flux_sim_0 - flux_ana_0
    diff_1 = flux_sim_1 - flux_ana_1

    # L2 error (absolute)
    L2_error_0 = np.sqrt(np.trapezoid(diff_0**2, z) / np.trapezoid(flux_ana_0**2, z))
    L2_errors_0.append(L2_error_0)

    L2_error_1 = np.sqrt(np.trapezoid(diff_1**2, z) / np.trapezoid(flux_ana_1**2, z))
    L2_errors_1.append(L2_error_1)

    # L∞ error (maximum absolute difference)
    Linf_error_0 = np.max(np.abs(diff_0))
    Linf_errors_0.append(Linf_error_0)

    Linf_error_1 = np.max(np.abs(diff_1))
    Linf_errors_1.append(Linf_error_1)

# Store errors in results
results['L2_errors_0'] = L2_errors_0
results['Linf_errors_0'] = Linf_errors_0
results['L2_errors_1'] = L2_errors_1
results['Linf_errors_1'] = Linf_errors_1

# Use most refined mesh for flux plotting
z_analytical = results['z_mesh'][-1]
flux_analytical_shape = np.sin(z_analytical * np.pi / fuelLength)
flux_analytical_0 = amplitude_0 * flux_analytical_shape
flux_analytical_1 = amplitude_1 * flux_analytical_shape

# Calculate k_eff relative errors
k_rel_errors = [(k_sim - k_analytical) / k_analytical * 100000 for k_sim in results['keff']]
results['k_rel_errors'] = k_rel_errors

print(f"\n{'='*60}")


#==============================================================================*
# Summary plots

# Plot 1: L2 Error convergence for both groups (stacked)
#-------------------------------------------------------

fig1, axes = plt.subplots(2, 1, figsize=(5, 6), dpi=200)

# L2 error convergence - Fast (top)
ax_L2_fast = axes[0]
ref_nz = results['nz_values'][0]
ref_L2 = results['L2_errors_0'][0]
ax_L2_fast.loglog(
    results['nz_values'],
    results['L2_errors_0'],
    'bo-',
    label="L2 error (fast)",
    markersize=8,
    linewidth=2
)
second_order_L2 = ref_L2 * (ref_nz / np.array(results['nz_values']))**2
first_order_L2 = ref_L2 * (ref_nz / np.array(results['nz_values']))**1
ax_L2_fast.loglog(results['nz_values'], second_order_L2, 'r--', label="2nd order", alpha=0.7)
ax_L2_fast.loglog(results['nz_values'], first_order_L2, 'g--', label="1st order", alpha=0.7)
ax_L2_fast.set_xlabel("Number of axial cells (nz) [-]")
ax_L2_fast.set_ylabel("L2 Error [-]")
ax_L2_fast.set_title("L2 Error Convergence (Group 0 - Fast)")
ax_L2_fast.legend()

# L2 error convergence - Thermal (bottom)
ax_L2_thermal = axes[1]
ref_L2_1 = results['L2_errors_1'][0]
ax_L2_thermal.loglog(
    results['nz_values'],
    results['L2_errors_1'],
    'bo-',
    label="L2 error (thermal)",
    markersize=8,
    linewidth=2
)
second_order_L2 = ref_L2_1 * (ref_nz / np.array(results['nz_values']))**2
first_order_L2 = ref_L2_1 * (ref_nz / np.array(results['nz_values']))**1
ax_L2_thermal.loglog(results['nz_values'], second_order_L2, 'r--', label="2nd order", alpha=0.7)
ax_L2_thermal.loglog(results['nz_values'], first_order_L2, 'g--', label="1st order", alpha=0.7)
ax_L2_thermal.set_xlabel("Number of axial cells (nz) [-]")
ax_L2_thermal.set_ylabel("L2 Error [-]")
ax_L2_thermal.set_title("L2 Error Convergence (Group 1 - Thermal)")
ax_L2_thermal.legend()

fig1.tight_layout()
fig1.savefig("fig_results_convergence_L2.png")
plt.close()


# Plot 2: k_eff convergence and relative error (stacked)
#-------------------------------------------------------

fig2, axes = plt.subplots(2, 1, figsize=(5, 6), dpi=200)

# k_eff convergence
ax_keff = axes[0]
ax_keff.plot(results['nz_values'], results['keff'], 'o-', linewidth=2, markersize=8, label='GeN-Foam')
ax_keff.axhline(y=k_analytical, color='r', linestyle='--', linewidth=2, label=f'Analytical ({k_analytical:.6f})')
ax_keff.set_xlabel("Number of axial cells (nz) [-]")
ax_keff.set_ylabel(r"k$_\text{eff}$")
ax_keff.set_title("k-effective Convergence")
ax_keff.legend()
ax_keff.set_xscale('log')
ax_keff.ticklabel_format(axis='y', style='plain', useOffset=False)

# k_eff relative error
ax_kerr = axes[1]
ax_kerr.plot(results['nz_values'], k_rel_errors, 'o-', linewidth=2, markersize=8, color='#2ca02c')
ax_kerr.axhline(y=0, color='k', linestyle='--', linewidth=0.8)
ax_kerr.set_xlabel("Number of axial cells (nz) [-]")
ax_kerr.set_ylabel("Relative error [pcm]")
ax_kerr.set_title("k-effective Relative Error")
ax_kerr.set_xscale('log')

fig2.tight_layout()
fig2.savefig("fig_results_convergence_keff.png")
plt.close()


# Plot 3: Relative errors for both groups (stacked)
#--------------------------------------------------

fig3, axes = plt.subplots(2, 1, figsize=(5, 8), dpi=200)

ax_flux, ax_rel_err = axes.flatten()

colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink"]

for i, nz in enumerate(nz_values):
    z = results['z_mesh'][i]
    flux0 = results['flux0'][i]
    flux1 = results['flux1'][i]
    ax_flux.plot(z, flux0, label=f'nz={nz}', alpha=0.7, color=colors[i])
    ax_flux.plot(z, flux1, ls='--', alpha=0.7, color=colors[i])

ax_flux.plot(z_analytical, flux_analytical_0, 'k-', linewidth=2, label='Analytical Fast')
ax_flux.plot(z_analytical, flux_analytical_1, 'r--', linewidth=2, label='Analytical Thermal')
ax_flux.set_ylabel("Neutron flux [n/(m$^2$ s)]")
ax_flux.set_title("Flux Distribution Comparison")
ax_flux.set_ylim(0)


# Relative errors vs axial position (group 0 = fast, top)
for i, nz in enumerate(nz_values):
    z = results['z_mesh'][i]
    flux_sim_0 = results['flux0'][i]

    # Analytical solution at this mesh
    flux_analytical_shape = np.sin(z * np.pi / fuelLength)
    flux_ana_0 = amplitude_0 * flux_analytical_shape

    # Relative error (no masking)
    rel_error_0 = (flux_sim_0 - flux_ana_0) / flux_ana_0 * 100
    ax_rel_err.plot(z, rel_error_0, label=f'nz={nz}', alpha=0.7, color=colors[i])

ax_rel_err.set_ylabel("Relative error [%]")


# Relative errors vs axial position (group 1 = thermal, bottom)
for i, nz in enumerate(nz_values):
    z = results['z_mesh'][i]
    flux_sim_1 = results['flux1'][i]

    # Analytical solution at this mesh
    flux_analytical_shape = np.sin(z * np.pi / fuelLength)
    flux_ana_1 = amplitude_1 * flux_analytical_shape

    # Relative error (no masking)
    rel_error_1 = (flux_sim_1 - flux_ana_1) / flux_ana_1 * 100
    ax_rel_err.plot(z, rel_error_1, ls='--', alpha=0.7, color=colors[i])


for ax in axes.flatten():
    ax.set_xlabel("Axial position [m]")
    ax.legend()
    ax.set_xlim(0, fuelLength * 1.05)

fig3.tight_layout()
fig3.savefig("fig_results_convergence_flux.png")
plt.close()


#==============================================================================*
# Print summary table

print("\nMesh Study Summary:")
print(f"Analytical k_eff = {k_analytical:.6f}\n")

table_data = []
for i, nz in enumerate(nz_values):
    table_data.append([
        nz,
        f"{results['keff'][i]:.6f}",
        f"{k_rel_errors[i]:.1f}",
        f"{results['L2_errors_0'][i]:.6e}",
        f"{results['L2_errors_1'][i]:.6e}"
    ])

print(tabulate(
    table_data,
    headers=['nz', 'k_eff', 'k_err [pcm]', 'L2 err g=0', 'L2 err g=1'],
    tablefmt='pipe'
))


#==============================================================================*
