"""

"""
#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh

import matplotlib.pyplot as plt


#==============================================================================*

def generate_mesh(nCells: int=10):
    nMesh = mesh.BlockMesh(region="neutroRegion")

    nMesh.createSphere(
        name='uranium',
        radius=0.087455,
        nCenter=nCells,
        nBorder=nCells,
        isAddBoundaryConditions=True
    )

    nMesh.mergeBoundaryFaces(["Wall"], "wall")

    return(nMesh)


def generate_model():
    #==========================================================================*
    # Mesh

    nMesh = generate_mesh()


    #==========================================================================*
    # Fields

    timeFolder0 = ffn.TimeFolder(0)

    defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
    defaultFlux.dimensions = ffn.Dimension(default='flux')
    defaultFlux.internalField = 1
    defaultFlux.set_boundary_condition("wall", bc.InletOutlet(phi="facePhi", inletValue=0, value=0))

    defaultFlux2 = ffn.Field("defaultFlux2", region=nMesh.region)
    defaultFlux2.dimensions = ffn.Dimension(default='flux')
    defaultFlux2.internalField = 1
    defaultFlux2.set_boundary_condition("wall", bc.FixedValue(value=0))


    timeFolder0.append(defaultFlux)
    timeFolder0.append(defaultFlux2)


    #==========================================================================*
    # Solvers

    neutronicsSolver = ffn.NeutronicsSolver(
        region=nMesh.region,
        solver="SNNeutronics",
        eigenvalueNeutronics=True,
        keff=1,
        power=1e+06,
        mesh=nMesh
    )

    neutronicsSolver.quadratureSet = ffn.QuadratureSet(
        region=nMesh.region,
        default=4
    )
    neutronicsSolver.quadratureSet.plot_directions()

    neutronicsSolver.neutronTransportOptions.integralPredictor = True
    neutronicsSolver.neutronTransportOptions.aitkenAcceleration = True
    neutronicsSolver.neutronTransportOptions.maxNeutronIterations = 50
    neutronicsSolver.neutronTransportOptions.neutronIterationResidual = 1e-6

    # Parallelization
    neutronicsSolver.decomposeParDict.numberOfSubdomains = 8
    neutronicsSolver.decomposeParDict.method = "scotch"

    # Nuclear data
    neutronicsSolver.nuclearData.import_from_openfoam("XS/nuclearData")

    refState = neutronicsSolver.nuclearData.get_state_by_name("reference")

    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    refState.plot_xs_spectrum(ax, zoneName='uranium', xsType='sigmaRemoval', label="$\Sigma_r$")
    refState.plot_xs_spectrum(ax, zoneName='uranium', xsType='nuSigmaEff', label=r"$\nu \Sigma_{f}$")
    ax.set_xlabel('Energy bins')
    ax.set_ylabel('Macroscopic XS [1/m]')
    ax.set_yscale('log')
    ax.legend()
    fig.tight_layout()
    fig.savefig("fig_xs_uranium_sigmaR_sigmaNuFiss.png")
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
    refState.plot_xs_spectrum(ax, zoneName='uranium', xsType='chiPrompt', label="$\chi_{prompt}$")
    refState.plot_xs_spectrum(ax, zoneName='uranium', xsType='chiDelayed', label="$\chi_{delayed}$")
    ax.set_xlabel('Energy bins')
    ax.set_ylabel('$\chi$')
    ax.set_yscale('log')
    ax.legend()
    fig.tight_layout()
    fig.savefig("fig_xs_uranium_chi.png")
    plt.close()

    uraniumZone = refState.get_zone_by_name("uranium")

    for i in range(6):
        fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
        if (i == 0):
            colormap = uraniumZone.plot_scattering_matrix(
                ax,
                xsType=f"scatteringMatrixP{i}",
                isLogScale=True,
                zMin=1e-9
            )
        else:
            colormap = uraniumZone.plot_scattering_matrix(
                ax,
                xsType=f"scatteringMatrixP{i}",
                isSymLogScale=True,
                linthresh=1e-6
            )
        ax.set_xlabel('Energy bins')
        ax.set_ylabel('Energy bins')
        fig.colorbar(colormap, ax=ax, label="Macroscopic XS [1/m]", location='right')
        fig.tight_layout()
        fig.savefig(f"fig_xs_uranium_P{i}.png")
        plt.close()

    #==========================================================================*
    # Settings

    model = ffn.Model()

    model.solvers.append(neutronicsSolver)
    model.timeFolders = [timeFolder0]

    settings: ffn.ControlDict = model.settings

    settings.application = 'GeN-Foam'
    settings.endTime = 20
    settings.deltaT = 1
    settings.writeControl = 'runTime'
    settings.writeInterval = settings.endTime
    settings.writePrecision = 6
    settings.timePrecision = 6
    settings.runTimeModifiable = False
    settings.adjustTimeStep = False

    return(model, neutronicsSolver, defaultFlux, defaultFlux2)


#==============================================================================*

def mainSolverComparison():
    model, neutronicsSolver, defaultFlux, defaultFlux2 = generate_model()

    #==========================================================================*
    # Run

    neutronicsSolvers = ["diffusionNeutronics", "SP3Neutronics", "SNNeutronics"]

    results = {}
    for solver in neutronicsSolvers:
        neutronicsSolver.solver = solver

        # Change BC
        if (solver != "SNNeutronics"):
            defaultFlux.set_boundary_condition("wall", bc.AlbedoSP3(value=0, alpha=0, forSecondMoment=False))
            defaultFlux2.set_boundary_condition("wall", bc.AlbedoSP3(value=0, alpha=0, forSecondMoment=True))

        # Allow low memory only for SN
        neutronicsSolver.nuclearData.isLowMemory = solver == "SNNeutronics"

        # Clear case and re-export
        ffn.allclean()
        model.export_to_openfoam()

        # Run
        ffn.run_preprocessing(model=model)

        ffn.run(model=model)

        if (model.is_parallel):
            ffn.run_reconstruction(model=model, isLatestTime=True)


        # On-the-fly Post-processing
        logResults = model.get_keff_from_log()

        results[solver] = {
            'time': logResults['time'],
            'keff': logResults['keff']
        }
        print(f"{solver}: keff = {logResults['keff'][-1]}")


    #==========================================================================*
    # Post-processing

    # Plot mesh
    model.plot_mesh(region=neutronicsSolver.mesh, show_edges=True)
    model.plot_mesh(region=neutronicsSolver.mesh, normal='y', show_edges=True)
    model.plot_slice(region=neutronicsSolver.mesh, show_edges=True)


    # Plot keff evolution
    rho = lambda k: (k-1)/k

    fig, ax = plt.subplots(dpi=200, figsize=(5, 4))
    ax2 = ax.twinx()

    for solver in neutronicsSolvers:
        time = results[solver]['time']
        keff = results[solver]['keff']

        drhodt = [
            1e5 * abs(rho(kf) - rho(ki)) / (tf - ti)
            for ti, tf, ki, kf in zip(time[:-1], time[1:], keff[:-1], keff[1:])
        ]

        ax.plot(time, keff, label=solver)
        ax2.plot(
            time[1:],
            drhodt,
            ls='--'
        )

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('k$_{eff}$')
    ax.set_ylim(0.8)
    ax.legend(
        bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
        ncols=2, mode="expand", borderaxespad=0.
    )
    ax2.set_ylabel(r'$\frac{d \rho}{dt}$ [pcm/s]')
    ax2.set_yscale('log')
    fig.tight_layout()
    fig.savefig(f"fig_results_keff.png")

    # Plot flux slice for the last solver
    for i in range(23):
        model.plot_slice(
            region=neutronicsSolver.mesh,
            time=model.settings.endTime,
            fieldName=f'flux{i}',
            cmap='Blues_r',
            unit='n/m2/s'
        )


#==============================================================================*

def mainMeshConvergence():
    # Create base case
    model, neutronicsSolver, _, _ = generate_model()

    print(model)

    results = {}
    for nCells in [3, 6, 10, 15, 20]:
        # Generate a new mesh
        neutronicsSolver.mesh = generate_mesh(nCells=nCells)

        model.caseFolder = f"Mesh{nCells}"

        # Export
        model.export_to_openfoam()

        # Run
        ffn.run_preprocessing(model=model)

        nMaxCell = neutronicsSolver.decomposeParDict.extract_info_from_log(caseFolder=model.caseFolder)["nMaxCells"]
        print(f"Max number of cells in 1 processor: {nMaxCell}")

        ffn.run(model=model)

        if (model.is_parallel):
            ffn.run_reconstruction(model=model, isLatestTime=True)


        # On-the-fly Post-processing
        logResults = model.get_keff_from_log()

        results[nCells] = {
            'time': logResults['time'],
            'keff': logResults['keff']
        }
        print(f"Ncells {nCells}: keff = {logResults['keff'][-1]}")


        # Plot keff over time for various mesh refinement
        fig, ax = plt.subplots(dpi=200, figsize=(5, 4))
        for nCells, values in results.items():
            time = values['time']
            keff = values['keff']
            ax.plot(time, keff, label=f"{nCells} cells")

        ax.set_xlabel('Time [s]')
        ax.set_ylabel('k$_{eff}$')
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"fig_results_keff.png")
        plt.close()


        # Plot keff over mesh refinement
        fig, ax = plt.subplots(dpi=200, figsize=(5, 4))
        ax.plot(list(results.keys()), [values['keff'][-1] for values in results.values()], marker=".")
        ax.set_xlabel('N cells')
        ax.set_ylabel('k$_{eff}$')
        ax.set_xlim(0)
        fig.tight_layout()
        fig.savefig(f"fig_results_keff_convergence.png")
        plt.close()


#==============================================================================*
# Main

if __name__ == "__main__":

    # mainSolverComparison()

    mainMeshConvergence()


#==============================================================================*
