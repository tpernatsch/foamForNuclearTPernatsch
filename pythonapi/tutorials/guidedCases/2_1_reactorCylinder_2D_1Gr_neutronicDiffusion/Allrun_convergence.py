"""
2D Cylindrical Reactor (R-Z geometry) - Mesh Convergence Study
"""
#==============================================================================*
# Imports

from Allrun import *

from tabulate import tabulate


#==============================================================================*
# Parameters

# Mesh resolution configurations
nr_values = [5, 11, 21, 51, 201]  # Radial cells
nz_values = [5, 11, 21, 51, 201]  # Axial cells
wedgeAngle = 1  # Wedge angle [degrees]

# Storage for results
results = {
    'nr': [],
    'nz': [],
    'keff': [],
    'L2_error': [],
    'Linf_error': [],
    'radial_profiles': [],
    'axial_profiles': [],
    'radial_errors': [],
    'axial_errors': [],
    'r_coords': [],
    'z_coords': []
}


#==============================================================================*
# Run simulations for different mesh resolutions

for nr in nr_values:
    for nz in nz_values:
        print(f"\n{'='*80}")
        print(f"Running simulation with nr={nr}, nz={nz}")
        print(f"{'='*80}\n")

        # Recreate mesh with new nr and nz values
        nMesh, wedgeFaces = createMesh(nr=nr, nz=nz)

        # Update solver with new mesh
        neutronicsSolver.mesh = nMesh

        # Clean and export
        ffn.allclean()
        model.export_to_openfoam()

        # Run
        ffn.run(model, is_preprocessing=True)

        # Post-processing
        X, Y, Z = fluidfoam.readmesh(".", region=nMesh.region, structured=False)
        flux0 = fluidfoam.readscalar(".", time_name="1", name="flux0", region=nMesh.region, structured=False)
        R = np.sqrt(X**2 + Y**2)

        # Normalize and compute errors
        flux_anal = analytical_flux(R, Z)
        flux0_norm, flux_anal_norm = flux0 / np.max(flux0), flux_anal / np.max(flux_anal)
        mask = flux_anal_norm > threshold
        rel_error = (flux0_norm[mask] - flux_anal_norm[mask]) / flux_anal_norm[mask] * 100
        L2_error = np.sqrt(np.mean(rel_error**2))
        Linf_error = np.max(np.abs(rel_error))

        # Reshape to 2D grid and extract profiles
        sort_idx = np.lexsort((Z, R))
        flux0_grid = flux0_norm[sort_idx].reshape(nr, nz).T
        flux_anal_grid = flux_anal_norm[sort_idx].reshape(nr, nz).T
        R_grid = R[sort_idx].reshape(nr, nz).T
        Z_grid = Z[sort_idx].reshape(nr, nz).T

        mid_z_idx, mid_r_idx = nz // 2, nr // 2

        # Radial profile at mid-height
        r_coords = R_grid[mid_z_idx, :]
        flux_r = flux0_grid[mid_z_idx, :]
        flux_r_anal = flux_anal_grid[mid_z_idx, :]
        mask_r = flux_r_anal > threshold
        rel_error_r = np.full_like(flux_r, np.nan)
        rel_error_r[mask_r] = (flux_r[mask_r] - flux_r_anal[mask_r]) / flux_r_anal[mask_r] * 100

        # Axial profile at mid-radius
        z_coords = Z_grid[:, mid_r_idx]
        flux_z = flux0_grid[:, mid_r_idx]
        flux_z_anal = flux_anal_grid[:, mid_r_idx]
        mask_z = flux_z_anal > threshold
        rel_error_z = np.full_like(flux_z, np.nan)
        rel_error_z[mask_z] = (flux_z[mask_z] - flux_z_anal[mask_z]) / flux_z_anal[mask_z] * 100

        # Store results
        results['nr'].append(nr)
        results['nz'].append(nz)
        results['keff'].append(model.keff())
        results['L2_error'].append(L2_error)
        results['Linf_error'].append(Linf_error)
        results['radial_profiles'].append(flux_r)
        results['axial_profiles'].append(flux_z)
        results['radial_errors'].append(rel_error_r)
        results['axial_errors'].append(rel_error_z)
        results['r_coords'].append(r_coords)
        results['z_coords'].append(z_coords)


print(f"\n{'='*80}")

#==============================================================================*
# Print summary table

print("\nMesh Convergence Study Summary:")

table_data_keff = [
    [nr_i] + [
        f"{results['keff'][i*len(nr_values) + j]:.6f}"
        for j in range(len(nz_values))
    ]
    for i, nr_i in enumerate(nr_values)
]
table_data_L2 = [
    [nr_i] + [
        f"{results['L2_error'][i*len(nr_values) + j]:.6e}"
        for j in range(len(nz_values))
    ]
    for i, nr_i in enumerate(nr_values)
]

print(tabulate(
    table_data_keff,
    headers=[r'k_eff (nr\nz)'] + nz_values,
    tablefmt='pipe'
))

print()

print(tabulate(
    table_data_L2,
    headers=[r'L2 error [%] (nr\nz)'] + nz_values,
    tablefmt='pipe'
))

print(f"\nAnalytical k_eff: {keff_anal:.6f}")


#==============================================================================*
# Summary plots

# Get analytical profiles for reference using finest mesh
finest_idx = [
    i
    for i in range(len(results['nr']))
    if results['nr'][i] == nr_values[-1] and results['nz'][i] == nz_values[-1]
][0]
r_ref = results['r_coords'][finest_idx]
z_ref = results['z_coords'][finest_idx]

# Normalize analytical profiles consistently with simulation
R_ref, Z_ref = np.meshgrid(r_ref, z_ref)
max_flux = np.max(analytical_flux(R_ref.flatten(), Z_ref.flatten()))
flux_r_anal_ref = analytical_flux(r_ref, fuelHeight/2) / max_flux
flux_z_anal_ref = analytical_flux(rOut/2, z_ref) / max_flux


#==============================================================================*
# Plot 1: Radial profiles for different resolutions (fixed nz)

indices = [i for i, nz in enumerate(results['nz']) if nz == nz_values[-1]]

fig1, (ax1a, ax1b) = plt.subplots(nrows=2, figsize=(5, 6), sharex=True, dpi=200)
ax1a.plot(r_ref, flux_r_anal_ref, 'k-', label='Analytical', linewidth=3, zorder=10)
for idx in indices:
    ax1a.plot(
        results['r_coords'][idx],
        results['radial_profiles'][idx],
        'o-',
        label=f"nr={results['nr'][idx]}",
        markersize=4,
        alpha=0.7
    )
ax1a.set_xlim(0, rOut)
ax1a.set_ylim(bottom=0)
ax1a.set_ylabel('Normalized flux [a.u.]')
ax1a.set_title(f'Radial profiles at Z = {z_ref[len(z_ref)//2]:.3f} m (nz={nz_values[-1]})')
ax1a.legend(loc='best')
ax1a.grid(True, alpha=0.3)

for idx in indices:
    ax1b.plot(
        results['r_coords'][idx],
        results['radial_errors'][idx],
        'o-',
        label=f"nr={results['nr'][idx]}",
        markersize=4
    )
ax1b.set_xlim(0, rOut)
ax1b.set_xlabel('Radial position [m]')
ax1b.set_ylabel('Relative error [%]')
ax1b.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
ax1b.legend(loc='best')
ax1b.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig("fig_results_radial_profiles_convergence.png")
plt.close()


#==============================================================================*
# Plot 2: Axial profiles for different resolutions (fixed nr)

indices = [i for i, nr in enumerate(results['nr']) if nr == nr_values[-1]]

fig2, (ax2a, ax2b) = plt.subplots(nrows=2, figsize=(5, 6), sharex=True, dpi=200)
ax2a.plot(z_ref, flux_z_anal_ref, 'k-', label='Analytical', linewidth=3, zorder=10)
for idx in indices:
    ax2a.plot(
        results['z_coords'][idx],
        results['axial_profiles'][idx],
        'o-',
        label=f"nz={results['nz'][idx]}",
        markersize=4,
        alpha=0.7
    )
ax2a.set_xlim(0, fuelHeight)
ax2a.set_ylabel('Normalized flux [a.u.]')
ax2a.set_title(f'Axial profiles at R = {r_ref[len(r_ref)//2]:.3f} m (nr={nr_values[-1]})')
ax2a.legend(loc='best')
ax2b.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
ax2a.grid(True, alpha=0.3)

for idx in indices:
    ax2b.plot(results['z_coords'][idx], results['axial_errors'][idx],
              'o-', label=f"nz={results['nz'][idx]}", markersize=4)
ax2b.set_xlim(0, fuelHeight)
ax2b.set_xlabel('Axial position [m]')
ax2b.set_ylabel('Relative error [%]')
ax2b.legend(loc='best')
ax2b.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig("fig_results_axial_profiles_convergence.png")
plt.close()


#==============================================================================*
# Plot 3: L2 error convergence (log-log plot)

# For each fixed nz, plot L2 vs delta_r
fig3, ax3 = plt.subplots(figsize=(5, 4), dpi=200)

# Plot for different fixed nz values
for fixed_nz_val in sorted(set(results['nz'])):
    indices = [i for i, nz in enumerate(results['nz']) if nz == fixed_nz_val]
    nr_subset = [results['nr'][i] for i in indices]
    L2_subset = [results['L2_error'][i] for i in indices]

    ax3.loglog(
        nr_subset,
        L2_subset,
        'o-',
        label=f'nz={fixed_nz_val}',
        markersize=8,
        linewidth=2
    )

# Add reference lines for different convergence orders
# Use most refined nz for reference
ref_nz = nz_values[-1]
nr_ref = np.array([
    results['nr'][i]
    for i in range(len(results['nr']))
    if results['nz'][i] == ref_nz]
)
if len(nr_ref) > 0:
    # Find L2 for the coarsest nr at refined nz
    ref_idx = [
        i
        for i in range(len(results['nr']))
        if results['nr'][i] == nr_ref[0] and results['nz'][i] == ref_nz
    ][0]
    L2_ref = results['L2_error'][ref_idx]
    nr_ref_val = nr_ref[0]

    # First order line
    ax3.loglog(
        nr_ref,
        L2_ref * (nr_ref_val / nr_ref)**1,
        'k--',
        label='1st order',
        linewidth=1.5,
        alpha=0.7
    )
    # Second order line
    ax3.loglog(
        nr_ref,
        L2_ref * (nr_ref_val / nr_ref)**2,
        'k:',
        label='2nd order',
        linewidth=1.5,
        alpha=0.7
    )

ax3.set_xlabel('Number of radial cells (nr) [-]')
ax3.set_ylabel('L2 relative error [%]')
ax3.set_title('L2 Error Convergence (Radial Direction)')
ax3.legend(loc='best')
ax3.grid(True, which='both', alpha=0.3)

fig3.tight_layout()
fig3.savefig("fig_results_L2_convergence_radial.png")
plt.close()


#==============================================================================*
# Plot 4: L2 error convergence in axial direction

fig4, ax4 = plt.subplots(figsize=(5, 4), dpi=200)

# For each fixed nr, plot L2 vs nz
for fixed_nr_val in sorted(set(results['nr'])):
    indices = [i for i, nr in enumerate(results['nr']) if nr == fixed_nr_val]
    nz_subset = [results['nz'][i] for i in indices]
    L2_subset = [results['L2_error'][i] for i in indices]

    ax4.loglog(
        nz_subset,
        L2_subset,
        's-',
        label=f'nr={fixed_nr_val}',
        markersize=8,
        linewidth=2
    )

# Add reference lines for different convergence orders
# Use most refined nr for reference
ref_nr = nr_values[-1]
nz_ref = np.array([
    results['nz'][i]
    for i in range(len(results['nz']))
    if results['nr'][i] == ref_nr
])
if len(nz_ref) > 0:
    # Find L2 for the coarsest nz at refined nr
    ref_idx = [
        i
        for i in range(len(results['nz']))
        if results['nz'][i] == nz_ref[0] and results['nr'][i] == ref_nr
    ][0]
    L2_ref = results['L2_error'][ref_idx]
    nz_ref_val = nz_ref[0]

    # First order line
    ax4.loglog(
        nz_ref,
        L2_ref * (nz_ref_val / nz_ref)**1,
        'k--',
        label='1st order',
        linewidth=1.5,
        alpha=0.7
    )
    # Second order line
    ax4.loglog(
        nz_ref,
        L2_ref * (nz_ref_val / nz_ref)**2,
        'k:',
        label='2nd order',
        linewidth=1.5,
        alpha=0.7
    )

ax4.set_xlabel('Number of axial cells (nz) [-]')
ax4.set_ylabel('L2 relative error [%]')
ax4.set_title('L2 Error Convergence (Axial Direction)')
ax4.legend(loc='best')
ax4.grid(True, which='both', alpha=0.3)

fig4.tight_layout()
fig4.savefig("fig_results_L2_convergence_axial.png")
plt.close()


#==============================================================================*
# Plot 5: k_eff convergence

# Reshape k_eff into 2D grid
nr_unique = sorted(set(results['nr']))
nz_unique = sorted(set(results['nz']))
keff_grid = np.zeros((len(nz_unique), len(nr_unique)))

for i, nz in enumerate(nz_unique):
    for j, nr in enumerate(nr_unique):
        keff_grid[i, j] = results['keff'][
            [k for k in range(len(results['nr']))
             if results['nr'][k] == nr and results['nz'][k] == nz][0]
        ]

fig5, ax5 = plt.subplots(figsize=(5, 4), dpi=200)

levels = np.r_[np.linspace(keff_anal-0.004, keff_anal, 10, endpoint=False),
               keff_anal,
               np.linspace(keff_anal, keff_anal+0.004, 10, endpoint=False)[1:]]

im = ax5.contourf(nr_unique, nz_unique, keff_grid,
                  levels=levels, cmap='RdBu_r')

ax5.contour(nr_unique, nz_unique, keff_grid,
            levels=levels, colors='black', alpha=0.3, linewidths=1)

ax5.contour(nr_unique, nz_unique, keff_grid,
            levels=[keff_anal], colors='black',
            linewidths=3, linestyles='--')

ax5.set(
    xlabel='Number of radial cells (nr) [-]',
    ylabel='Number of axial cells (nz) [-]',
    title='$k_{eff}$ (Analytical: '+f'{keff_anal:.6f})',
    xscale='log',
    yscale='log'
)

plt.colorbar(im, ax=ax5, label='$k_{eff}$')
fig5.tight_layout()
fig5.savefig("fig_results_keff_convergence.png")
plt.close()


#==============================================================================*
# Plot 6: 2D contour of L2 error as function of nr and nz

# Reshape L2 error into 2D grid
nr_unique = sorted(set(results['nr']))
nz_unique = sorted(set(results['nz']))
L2_grid = np.zeros((len(nz_unique), len(nr_unique)))

for i, nz_val in enumerate(nz_unique):
    for j, nr_val in enumerate(nr_unique):
        idx = [k for k in range(len(results['nr']))
               if results['nr'][k] == nr_val and results['nz'][k] == nz_val][0]
        L2_grid[i, j] = results['L2_error'][idx]

fig6, ax6 = plt.subplots(figsize=(5, 4), dpi=200)
im = ax6.contourf(nr_unique, nz_unique, L2_grid, levels=20, cmap='viridis')
ax6.contour(nr_unique, nz_unique, L2_grid, levels=20, colors='white', alpha=0.3, linewidths=1.5)
ax6.set_xlabel('Number of radial cells (nr) [-]')
ax6.set_ylabel('Number of axial cells (nz) [-]')
ax6.set_title('L2 Relative Error [%]')
ax6.set_xscale('log')
ax6.set_yscale('log')
plt.colorbar(im, ax=ax6, label='L2 error [%]')

fig6.tight_layout()
fig6.savefig("fig_results_L2_error_2D.png")
plt.close()

#==============================================================================*
