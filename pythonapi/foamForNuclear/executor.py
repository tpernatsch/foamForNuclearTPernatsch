import os
import subprocess

import foamlib

from foamForNuclear.common import copyFolder
from foamForNuclear.mesh.blockMesh import BlockMesh
from foamForNuclear.mesh.polyMesh import PolyMesh
from foamForNuclear.mesh.unvMesh import UnvMesh
from foamForNuclear.model import Model


def _process_CLI_arguments(
        ffn_exec: str,
        mpi_args=None,
    ):
    args = [ffn_exec]

    if mpi_args is not None:
        args = mpi_args + args + ['-parallel']

    return(args)


def _run(args, applicationName, output: bool=False, cwd: str='.'):
    # Launch a subprocess
    if (output):
        p = subprocess.Popen(
            args,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=True
        )

        # Capture and re-print subprocess output in real-time
        lines = []
        while True:
            # If subprocess is finished, break loop
            line = p.stdout.readline()
            if not line and p.poll() is not None:
                break

            lines.append(line)
            if (output and "Time " in line):
                # If user requested output, print to screen
                print(line, end='')

        # Raise an exception if return status is non-zero
        if p.returncode != 0:
            # Get error message from output and simplify whitespace
            output = ''.join(lines)
            if 'ERROR: ' in output:
                _, _, error_msg = output.partition('ERROR: ')
            elif 'what()' in output:
                _, _, error_msg = output.partition('what(): ')
            else:
                error_msg = 'FFN aborted unexpectedly.'
            error_msg = ' '.join(error_msg.split())

            raise RuntimeError(error_msg)
    else:
        print(f"FFN run: {applicationName}...")

        log = open(f'log.{applicationName}', 'a')

        p = subprocess.run(
            args,
            cwd=cwd,
            stdout=log,
            stderr=log,
            universal_newlines=True,
            shell=True
        )

        print("FFN run completed")


def run_preprocessing(
        model: Model,
        verbose: bool=False
    ):
    """
    Execute all the OpenFOAM pre-processing utilities:
    - blockMesh
    - ideasUnvToFoam
    - topoSet
    - createPatch
    - setFields
    - createBaffles
    - checkMesh
    - decomposePar
    - GeN-Foam -initializeMappedFields

    Parameters
    ----------
    model : Model
    verbose : bool
        Add `echo` commands in `Allrun.prerun` (default `False`)
    """
    commands = ". ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions\n"

    addVerbose = lambda text: f"echo {text}\n" if verbose else ""

    # Mesh generation
    for solver in model.solvers:
        region = solver.region
        # No region in the case
        if (region == ""):
            # Mesh generation
            if (isinstance(solver.mesh, BlockMesh)):
                commands += f"runApplication blockMesh -merge-points\n"
            elif (isinstance(solver.mesh, UnvMesh)):
                commands += f"runApplication ideasUnvToFoam\n"
            elif (isinstance(solver.mesh, PolyMesh)):
                commands += f"runApplication -s {region} renumberMesh -region {region} -overwrite\n"

            # Mesh manipulation
            if (not solver.mesh.topoSetDict.is_empty):
                commands += f"runApplication topoSet -noZero\n"
            if (not solver.mesh.createPatchDict.is_empty):
                commands += f"runApplication createPatch -overwrite\n"
            if (not solver.setFieldsDict.is_empty):
                commands += f"runApplication setFields\n"
            if (not solver.mesh.createBafflesDict.is_empty):
                commands += f"runApplication createBaffles -overwrite\n"

            commands += f"runApplication checkMesh\n"

        # Multiple regions
        else:
            # Mesh generation
            if (isinstance(solver.mesh, BlockMesh)):
                commands += f"runApplication -s {region} blockMesh -region {region} -merge-points\n"
            elif (isinstance(solver.mesh, UnvMesh)):
                commands += f"runApplication -s {region} ideasUnvToFoam\n"
                commands += f"rm -rf constant/{region}/polyMesh\n"
                commands += f"mv constant/polyMesh constant/{region}\n"
            elif (isinstance(solver.mesh, PolyMesh)):
                commands += f"runApplication -s {region} renumberMesh -region {region} -overwrite -no-fields\n"

            # Mesh manipulation
            if (not solver.mesh.topoSetDict.is_empty):
                commands += f"runApplication -s {region} topoSet -region {region}\n"
            if (not solver.mesh.createPatchDict.is_empty):
                commands += f"runApplication -s {region} createPatch -region {region} -overwrite\n"
            if (not solver.setFieldsDict.is_empty):
                commands += f"runApplication -s {region} setFields -region {region}\n"
            if (not solver.mesh.createBafflesDict.is_empty):
                commands += f"runApplication -s {region} createBaffles -region {region} -overwrite\n"

            commands += f"runApplication -s {region} checkMesh -region {region}\n"

    # Mesh decomposition
    if (model.is_parallel):
        commands += "runApplication decomposePar -allRegions -copyUniform\n"

    # GeN-Foam initialize mapped fields
    if (
        model.settings.application == "GeN-Foam"
        and len(model.solvers) >= 2
    ):
        commands += "GeN-Foam -initializeMappedFields > log.GeN-Foam.initializeMappedFields\n"
        commands += addVerbose("Run GeN-Foam -initializeMappedFields")

    # Change directory to case folder
    cwd = os.getcwd()
    os.chdir(model.caseFolder)

    # Dump commands
    with open('Allrun.prerun', 'w') as f:
        f.write('#!/bin/sh\n')
        f.write(commands)

    # Dump a .foam file to visualize in ParaView
    with open('foam.foam', 'w') as f:
        f.write('')

    # Run command
    p = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    out, err = p.communicate(commands)

    # Return to the original folder
    os.chdir(cwd)


def run_reconstruction(
        model: Model,
        isAllRegions: bool=True,
        isLatestTime: bool=False
    ):
    """
    Execute the `reconstructPar` OpenFOAM utility

    Parameters
    ----------
    model : Model
    isAllRegions : bool
        If true, add the `-allRegions` flag (default `True`)
    isLatestTime : bool
        If true, add the `-latestTime` flag, else add the `-newTimes` flag,
        default `False`
    """
    commands = ". ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions\n"

    text = "reconstructPar"
    if (isAllRegions):
        text += " -allRegions"
    if (isLatestTime):
        text += " -latestTime"
    else:
        text += " -newTimes"

    commands += text + '\n'

    # Change directory to case folder
    cwd = os.getcwd()
    os.chdir(model.caseFolder)

    # Dump commands
    with open('Allrun.postrun', 'w') as f:
        f.write('#!/bin/sh\n')
        f.write(commands)

    # Run command
    p = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    out, err = p.communicate(commands)

    # Return to the original folder
    os.chdir(cwd)



def run(
        model: Model,
        is_preprocessing: bool=False,
    ):
    """
    Command to run a FFN simulation.

    Parameters
    ----------
    model : Model
        FFN Model containing the solvers, time folders, coupling objects, ...
    is_preprocessing : bool
        Execute preprocessing programs such as mesh generation before running
        the main simulation (default `False`).
    """
    if (is_preprocessing):
        run_preprocessing(model)

    if (not model.is_parallel):
        case = foamlib.FoamCase(model.caseFolder)
        case.run()
    else:
        # Need to do it manually, foamLib is using mpiexec and it is not writing
        # in log file
        app = model.settings.application
        nCores = model.get_number_processors

        # Change directory to case folder
        cwd = os.getcwd()
        os.chdir(model.caseFolder)

        # Run command
        case = foamlib.FoamCase()
        case.run(
            cmd=f"mpirun -np {nCores} {app} -parallel | tee log.{app}",
            check=False,
            log=False
        )

        # Return to the original folder
        os.chdir(cwd)


def reconstructParallelCase(
        model: Model,
        isLatestTimeStep: bool=False
    ):
    commands = ". $WM_PROJECT_DIR/bin/tools/RunFunctions\n"
    commands += f"runApplication reconstructPar -allRegions -newTimes"

    if (isLatestTimeStep):
        commands += " -latestTime"

    commands += "\n"

    # Change directory to case folder
    cwd = os.getcwd()
    os.chdir(model.caseFolder)

    p = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    out, err = p.communicate(commands)

    # Return to the original folder
    os.chdir(cwd)



def generateCaseAsFMU(caseName: str, pythonFileName: str):
    """
    Generate an OpenFOAM case as an FMU. The `pythonFileName` must contain a
    class object that inherite from the `OF2Fmu` object.

    Parameters
    ----------
    caseName : str
        Name of the case
    pythonFileName : str
        Name of the Python file that contains the FMU class wrapper inhereted
        from `OF2Fmu`.
    """
    if (pythonFileName[-3:] == '.py'):
        msg = "Please remove the '.py' at the end of 'pythonFileName'"
        raise ValueError(msg)

    commands =  ". $WM_PROJECT_DIR/bin/tools/RunFunctions\n"

    commands += f"fmu4foam build -f {pythonFileName}.py -of {caseName} --no-external-tool\n"

    process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    out, err = process.communicate(commands)


def allclean():
    """
    Remove all OpenFOAM files and folder, and `*.png`, `*.gif`, `*_log.txt`,
    `*.csv`, and `input_check.txt` files.
    """
    commands =  ". ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions\n"
    commands += "cleanCase0\n"
    commands += "rm -rf constant system processor* log.* Allrun.pre* *.png *.gif *_log.txt *.csv input_check.txt\n"

    process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    out, err = process.communicate(commands)


def cleanCase0():
    """
    Clean the OpenFOAM case, keeping the OpenFOAM case structure
    """
    commands =  ". ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions\n"
    commands += "cleanCase0\n"

    process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    out, err = process.communicate(commands)


def duplicateFolder(folderToCopy, newFolder):
    copyFolder(folderToCopy, newFolder)
