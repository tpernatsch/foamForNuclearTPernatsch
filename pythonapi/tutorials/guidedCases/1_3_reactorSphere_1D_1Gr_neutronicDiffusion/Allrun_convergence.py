"""
Mesh refinement study for 1D sphere reactor.
Runs multiple simulations with different mesh resolutions.
"""

#==============================================================================*
# Imports

from Allrun import *

from tabulate import tabulate


#==============================================================================*
# Parameters

# Different mesh resolutions to test
nr_values = [11, 21, 51, 101, 201, 501, 1001, 5001, 10001]

# Storage for results
results = {
    'nr_values': nr_values,
    'keff': [],
    'flux0': [],
    'r_mesh': []
}


#==============================================================================*
# Loop over different mesh resolutions

for nr in nr_values:
    print(f"\n{'='*60}")
    print(f"Running simulation with nr={nr}")
    print(f"{'='*60}\n")

    # Recreate mesh with new nr value
    nMesh, _ = createMesh(dr=dr, nr=nr)

    # Update solver with new mesh
    neutronicsSolver.mesh = nMesh

    # Clean and export
    ffn.allclean()
    model.export_to_openfoam()

    # Run
    ffn.run(model, is_preprocessing=True)

    # Extract results
    keff = model.keff()
    print(f"keff = {keff}")

    # Extract mesh and flux data
    X, Y, Z = fluidfoam.readmesh('.', region=nMesh.region, structured=False)
    flux0 = fluidfoam.readscalar(
        '.',
        time_name=str(settings.endTime),
        name='flux0',
        region=nMesh.region,
        structured=False
    )

    # Calculate radial coordinate from Cartesian coordinates
    r_mesh = np.sqrt(X**2 + Y**2 + Z**2)

    # Store results
    results['keff'].append(keff)
    results['flux0'].append(flux0)
    results['r_mesh'].append(r_mesh)


#==============================================================================*
# Calculate analytical solution and errors for each mesh

# Analytical solution function
fluxTh = lambda r_: np.sin(r_ * np.pi/dr)/r_/np.pi*dr

# Threshold for relative error plot only (5% of max analytical flux)
threshold = 0.05

# Calculate errors for each mesh
L2_errors = []
Linf_errors = []

for i, nr in enumerate(nr_values):
    r = results['r_mesh'][i]
    flux_sim = results['flux0'][i]

    # Calculate analytical solution at this mesh's points
    flux_ana = fluxTh(r)

    # Normalize both solutions
    flux_sim_norm = flux_sim / np.max(flux_sim)
    flux_ana_norm = flux_ana / np.max(flux_ana)

    # Calculate errors over full domain
    diff = flux_sim_norm - flux_ana_norm

    # L2 error (normalized)
    L2_error = np.sqrt(np.trapezoid(diff**2, r) / np.trapezoid(flux_ana_norm**2, r))
    L2_errors.append(L2_error)

    # L∞ error (maximum absolute difference)
    Linf_error = np.max(np.abs(diff))
    Linf_errors.append(Linf_error)

# Store errors in results
results['L2_errors'] = L2_errors
results['Linf_errors'] = Linf_errors

# Use most refined mesh for flux plotting
r_analytical = results['r_mesh'][-1]
flux_analytical_norm = fluxTh(r_analytical) / np.max(fluxTh(r_analytical))

print(f"\n{'='*60}")


#==============================================================================*
# Summary plots

# Plot 1: k_eff convergence and relative errors
#----------------------------------------------

fig1, ax_keff = plt.subplots(figsize=(5, 4), dpi=200)

# k_eff convergence
ax_keff.plot(results['nr_values'], results['keff'], 'o')
ax_keff.set_xlabel("Number of radial cells (nr) [-]")
ax_keff.set_ylabel(r"k$_\text{eff}$")
ax_keff.set_title("k-effective Convergence")
ax_keff.set_xscale('log')

fig1.tight_layout()
fig1.savefig("fig_results_keff.png")
plt.close()


# Plot 2: Flux distribution (standalone)
#---------------------------------------

fig_flux, (ax_flux, ax_rel_err) = plt.subplots(2, 1, figsize=(5, 6), dpi=200)

# Flux vs radial position
for i, nr in enumerate(nr_values):
    r = results['r_mesh'][i]
    flux = results['flux0'][i]
    normFlux = flux / max(flux)
    ax_flux.plot(r, normFlux, 'o-', label=f'nr={nr}', alpha=0.7, markersize=4)

ax_flux.plot(r_analytical, flux_analytical_norm, 'k--', linewidth=2, label='Analytical')
ax_flux.set_xlabel("Radial position [m]")
ax_flux.set_ylabel("Normalized neutron flux [a.u]")
ax_flux.set_title("Flux Distribution Comparison")
ax_flux.legend(ncols=2)
ax_flux.set_xlim((0, dr))

# Relative errors vs radial position
for i, nr in enumerate(nr_values):
    r = results['r_mesh'][i]
    flux_sim_norm = results['flux0'][i] / np.max(results['flux0'][i])
    flux_ana_norm = fluxTh(r) / np.max(fluxTh(r))

    # Mask regions where analytical flux is too small (< 5% of max)
    mask = flux_ana_norm > threshold
    rel_error = (flux_sim_norm - flux_ana_norm) / flux_ana_norm * 100
    # Only plot where analytical flux is significant
    ax_rel_err.plot(r[mask], rel_error[mask], label=f'nr={nr}', alpha=0.7)

ax_rel_err.set_xlabel("Radial position [m]")
ax_rel_err.set_ylabel("Relative error [%]")
ax_rel_err.set_title("Relative Error vs Position")
ax_rel_err.legend(ncols=2)
ax_rel_err.set_xlim((0, dr))

fig_flux.tight_layout()
fig_flux.savefig("fig_results_flux_comparison.png")
plt.close()


# Plot 3: Error convergence (L2 and L∞)
#--------------------------------------

fig2, (ax_L2, ax_Linf) = plt.subplots(2, 1, figsize=(5, 6), dpi=200)

# L2 error convergence
ax_L2.loglog(
    results['nr_values'], results['L2_errors'],
    'bo-',
    label="L2 error",
    markersize=8,
    linewidth=2
)

# Add reference lines for convergence order
ref_nr = results['nr_values'][0]
ref_L2 = results['L2_errors'][0]
second_order_L2 = ref_L2 * (ref_nr / np.array(results['nr_values']))**2
first_order_L2 = ref_L2 * (ref_nr / np.array(results['nr_values']))**1
ax_L2.loglog(results['nr_values'], second_order_L2, 'r--', label="2nd order", alpha=0.7)
ax_L2.loglog(results['nr_values'], first_order_L2, 'g--', label="1st order", alpha=0.7)

ax_L2.set_xlabel("Number of radial cells (nr) [-]")
ax_L2.set_ylabel("L2 Error [-]")
ax_L2.set_title("L2 Error Convergence")
ax_L2.legend()

# L∞ error convergence
ax_Linf.loglog(
    results['nr_values'], results['Linf_errors'],
    'rs-',
    label="L∞ error",
    markersize=8,
    linewidth=2
)

# Add reference lines for convergence order
ref_Linf = results['Linf_errors'][0]
second_order_Linf = ref_Linf * (ref_nr / np.array(results['nr_values']))**2
first_order_Linf = ref_Linf * (ref_nr / np.array(results['nr_values']))**1
ax_Linf.loglog(results['nr_values'], second_order_Linf, 'r--', label="2nd order", alpha=0.7)
ax_Linf.loglog(results['nr_values'], first_order_Linf, 'g--', label="1st order", alpha=0.7)

ax_Linf.set_xlabel("Number of radial cells (nr) [-]")
ax_Linf.set_ylabel("L∞ Error [-]")
ax_Linf.set_title("L∞ Error Convergence")
ax_Linf.legend()

fig2.tight_layout()
fig2.savefig("fig_results_error_convergence.png")
plt.close()


#==============================================================================*
# Print summary table

print("\nMesh Study Summary:")
table_data = []
for i, nr in enumerate(nr_values):
    table_data.append([
        nr,
        f"{results['keff'][i]:.6f}",
        f"{results['L2_errors'][i]:.6e}",
        f"{results['Linf_errors'][i]:.6e}",
        len(results['r_mesh'][i])
    ])

print(tabulate(
    table_data,
    headers=['nr', 'k_eff', 'L2 error', 'L∞ error', '# cells'],
    tablefmt='pipe'
))


#==============================================================================*
