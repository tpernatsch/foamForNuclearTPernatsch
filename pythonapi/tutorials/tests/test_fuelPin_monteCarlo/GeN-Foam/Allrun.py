#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import foamForNuclear.boundaryConditions as bc
import openmc
import matplotlib.pyplot as plt


#==============================================================================*

ffn.allclean()


#==============================================================================*
# Mesh

nMesh = mesh.BlockMesh(region="neutroRegion")

fuelOR = 0.39218e-2
cladOR = 0.45720e-2
pitch = 2 * 0.63e-2
length = 0.001

nx = 50
ny = nx
nz = 1

nMesh.create_cylinder_along_z(
    'fuel',
    radius=fuelOR,
    lowZ=0, highZ=length,
    nx=nx, ny=ny, nz=nz,
    isAddBoundaryConditions=True
)

nMesh.create_ring_along_z(
    'cladding',
    innerRadius=fuelOR,
    outerRadius=cladOR,
    lowZ=0, highZ=length,
    nr=15, nt=nx, nz=nz,
    isAddBoundaryConditions=True
)

nMesh.create_cube_with_hole_along_z(
    'water',
    lowX=-pitch/2, highX=pitch/2,
    lowY=-pitch/2, highY=pitch/2,
    lowZ=0, highZ=length,
    radius=cladOR,
    nx=nx, ny=ny, nz=nz, nt=nx,
    isAddBoundaryConditions=True
)

nMesh.add_merge_patch_pairs()

nMesh.merge_patches_with_name('top', includeFacename=['Top'])
nMesh.merge_patches_with_name('bottom', includeFacename=['Bottom'])
nMesh.merge_patches_with_name('walls', includeFacename=['Wall'])


#==============================================================================*
# Fields

timeFolder0 = ffn.TimeFolder(time=0)

defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.Dimension(default='flux')
defaultFlux.internalField = 1e21
defaultFlux.set_boundary_condition("walls", bc.ZeroGradient())
defaultFlux.set_boundary_condition("top", bc.ZeroGradient())
defaultFlux.set_boundary_condition("bottom", bc.ZeroGradient())

defaultFlux2 = ffn.Field("defaultFlux2", region=nMesh.region)
defaultFlux2.dimensions = ffn.Dimension(default='flux')
defaultFlux2.internalField = 1e21
defaultFlux2.set_boundary_condition("walls", bc.ZeroGradient())
defaultFlux2.set_boundary_condition("top", bc.ZeroGradient())
defaultFlux2.set_boundary_condition("bottom", bc.ZeroGradient())

timeFolder0.append(defaultFlux)
timeFolder0.append(defaultFlux2)


#==============================================================================*
# Solver

neutronicsSolver = ffn.NeutronicsSolver(
    region=nMesh.region,
    # solver="diffusionNeutronics",
    solver="SP3Neutronics",
    # solver="SNNeutronics",
    mesh=nMesh,
    eigenvalueNeutronics=True,
    keff=1,
    power=30e3
)

# neutronicsSolver.quadratureSet = ffn.QuadratureSet(region=nMesh.region)
# neutronicsSolver.quadratureSet.set_default_4direction()
# neutronicsSolver.nuclearData.legendreMoments = 1


serpentState = ffn.NuclearDataState('serpent')
serpentState.read_from_serpent(
    outputFilename="../Serpent/main_res.m",
    universes=[
        ('uFuel', 'fuel', 1),
        ('uClad', 'cladding', 0),
        ('uWater', 'water', 0),
    ]
)

openmcState = ffn.NuclearDataState('reference')
openmcState.read_from_openmc(
    outputFilename="../OpenMC/statepoint.1000.h5",
    domains=[
        (openmc.Cell(name='Fuel'), 'fuel', 1),
        (openmc.Cell(name='Cladding'), 'cladding', 0),
        (openmc.Cell(name='Water'), 'water', 0),
    ],
    energyGroups=openmc.mgxs.EnergyGroups([0., 0.625, 20.0e6]),
    delayedGroups=list(range(1,7)),
    rootDomain=openmc.Universe(0, name='root universe')
)

neutronicsSolver.nuclearData.add_state(serpentState)
# neutronicsSolver.nuclearData.add_state(openmcState)


# Plot XS comparison
#-------------------

fig, ax = plt.subplots()
serpentState.plot_xs_spectrum(ax, 'fuel', 'sigmaRemoval', label=r'Serpent $\Sigma_r$')
serpentState.plot_xs_spectrum(ax, 'fuel', 'nuSigmaEff', label=r'Serpent $\nu \Sigma_f$')
openmcState.plot_xs_spectrum(ax, 'fuel', 'sigmaRemoval', label=r'OpenMC $\Sigma_r$', ls='--')
openmcState.plot_xs_spectrum(ax, 'fuel', 'nuSigmaEff', label=r'OpenMC $\nu \Sigma_f$', ls='--')
ax.set_xscale('log')
ax.set_xlim(1e-9, 20)
ax.set_xlabel('Energy [MeV]')
ax.legend()
fig.savefig("fig_xs.png")


#==============================================================================*
# Settings

solvers = ffn.Solvers([neutronicsSolver])

model = ffn.Model(solvers=solvers, timeFolders=[timeFolder0])

settings: ffn.ControlDict = model.settings

settings.application = 'GeN-Foam'
settings.endTime = 10
settings.deltaT = 1e-6
settings.writeInterval = settings.endTime
settings.writeControl = 'adjustableRunTime'
settings.adjustTimeStep = True
settings.maxPowerVariation = 0.01


#==============================================================================*
# Main

if __name__ == "__main__":
    print(model)

    # Add idxField to visualize cellZones, not used in calculations
    idxField = neutronicsSolver.create_zone_field(nMesh.cellZones)
    timeFolder0.append(idxField)

    # Export to OpenFOAM
    model.export_to_openfoam()


    # Preprocessing
    #--------------

    ffn.run_preprocessing(model=model)

    model.plot_mesh(region=nMesh, show_edges=True)
    model.plot_mesh(region=nMesh, show_edges=True, normal='z')
    model.plot_slice(
        region=nMesh,
        fieldName=idxField.name,
        show_edges=True,
        cmap="tab10",
        normal="z",
        limits=[0.5, 10.5]
    )


    # Run
    #----

    ffn.run(model=model)


    # Post-processing
    #----------------

    print(f"keff = {model.keff()}")

    model.plot_residuals(
        parameters=['fluxStar0', 'fluxStar1', 'angularFlux_1_0'],
        title='Neutronics'
    )

    model.plot_mesh(
        region=nMesh,
        time=settings.endTime,
        fieldName='oneGroupFlux',
        cmap='Blues',
        unit='n/m2/s',
    )
    model.plot_mesh(
        region=nMesh,
        time=settings.endTime,
        fieldName='powerDensity',
        cmap='inferno',
        unit='W/m3',
    )

    xSample, sampleOneGroup = model.sample_over_line(
        region=nMesh.region,
        time=settings.endTime,
        fieldName='oneGroupFlux',
        point1=(-pitch/2, 0, length/2),
        point2=(+pitch/2, 0, length/2),
    )
    _, sampleFluxFast = model.sample_over_line(
        region=nMesh.region,
        time=settings.endTime,
        fieldName='flux0',
        point1=(-pitch/2, 0, length/2),
        point2=(+pitch/2, 0, length/2),
    )
    _, sampleFluxTherm = model.sample_over_line(
        region=nMesh.region,
        time=settings.endTime,
        fieldName='flux1',
        point1=(-pitch/2, 0, length/2),
        point2=(+pitch/2, 0, length/2),
    )

    xSample = [x*1e2 for x in xSample[:,0]]

    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax2 = ax1.twinx()
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(("axes", 1.3))
    line1,  = ax1.plot(xSample, sampleOneGroup, label='Integrated', color='tab:green')
    line2,  = ax2.plot(xSample, sampleFluxFast, label='Fast', color='tab:blue')
    line3,  = ax3.plot(xSample, sampleFluxTherm, label='Thermal', color='tab:red')
    ax1.vlines(
        x=[x*1e2 for x in [-pitch/2, -cladOR, -fuelOR, fuelOR, cladOR, pitch/2]],
        ymin=0, ymax=max(sampleOneGroup)*1.1,
        color='tab:grey',
        alpha=0.5,
    )
    ax1.set_ylim((min(sampleOneGroup), max(sampleOneGroup)))
    ax1.set_xlabel("X [cm]")
    ax1.set_ylabel("Integrated neutron flux [neutron/m$^2$/s]")
    ax2.set_ylabel("Fast neutron flux [neutron/m$^2$/s]")
    ax3.set_ylabel("Thermal neutron flux [neutron/m$^2$/s]")
    ax1.legend(
        handles=[line1, line2, line3],
        bbox_to_anchor=(0, 1.02, 1, 0.2),
        loc="lower left",
        mode="expand",
        borderaxespad=0,
        ncol=3
    )
    fig.tight_layout()
    fig.savefig("fig_results_neutronFluxComparison.png")


#==============================================================================*
