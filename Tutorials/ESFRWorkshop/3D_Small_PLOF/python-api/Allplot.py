"""
Post-processing script
"""
# Imports
from Allrun import *


model.plot_residuals(
    parameters=['fluxStar0'],
    title="Neutronics"
)
model.plot_residuals(
    parameters=['p_rgh', 'h'],
    title="Thermal-hydraulics"
)
model.plot_residuals(
    parameters=['fuelDisp', 'CRDisp', 'Dx', 'Dy', 'Dz'],
    title="Thermomechanics"
)


lastTime = model.get_time_steps()[-1]

print(f"keff = {model.keff(time=lastTime)}")

res = model.get_keff_from_log()

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(res['time'], res['keff'])
ax.set_xlabel("Time [s]")
ax.set_ylabel(r"k$_\text{eff}$")
fig.tight_layout()
fig.savefig("fig_results_keff.png")


# Plot flux along a line
points, flux0 = model.sample_over_line(
    region=nMesh.region,
    time=lastTime,
    fieldName='flux0',
    point1=(0, 0, -0.7),
    point2=(0, 0, 0.7)
)
z = points[:,2]

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(z, [e * 1e-4 for e in flux0])
ax.set_xlabel("Axial direction [m]")
ax.set_ylabel("Neutron flux [n/cm$^2$/s]")
fig.tight_layout()
fig.savefig("fig_results_flux0.png")



model.plot_mesh(
    region=thMesh,
    time=lastTime,
    fieldName='T',
    cmap='RdBu_r',
    unit='K'
)

model.plot_slice(
    region="neutroRegion",
    time=lastTime,
    fieldName='flux0',
    normal='z',
    cmap='Blues',
    unit='neutron/m2/s',
    show_edges=True
)
model.plot_slice(
    region="neutroRegion",
    time=lastTime,
    fieldName='powerDensity',
    normal='z',
    cmap='inferno',
    unit='W/m3',
    show_edges=True
)
for field in ['TCool', 'TClad', 'TFuel', 'TStruct', 'TStructMech']:
    model.plot_slice(
        region="neutroRegion",
        time=lastTime,
        fieldName=field,
        normal='x',
        cmap='RdBu_r',
        unit='K',
        show_edges=True
    )
    model.plot_slice(
        region="neutroRegion",
        time=lastTime,
        fieldName=field,
        normal='z',
        cmap='RdBu_r',
        unit='K',
        show_edges=True
    )

model.plot_animation(
    region=nMesh,
    fieldName='TCool',
    unit='K',
    cmap="RdBu_r",
    fps=1
)
