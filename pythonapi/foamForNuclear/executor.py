import os
import subprocess

import foamlib
from numpy import empty

from foamForNuclear.common import copyFolder, of_flavour
from foamForNuclear.mesh.blockMesh import BlockMesh, Face
from foamForNuclear.mesh.polyMesh import PolyMesh
from foamForNuclear.mesh.unvMesh import UnvMesh
from foamForNuclear.case import Case

import re
from pathlib import Path

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


# Read the boundary patches from boundary file
def read_boundary_patches(boundary_path: Path) -> dict[str, str]:
    txt = boundary_path.read_text(encoding="utf-8", errors="ignore")

    # isolate the patch list: first '(' ... last ')'
    i0 = txt.find('(')
    i1 = txt.rfind(')')
    if i0 == -1 or i1 == -1 or i1 <= i0:
        raise RuntimeError(f"Cannot locate patch list parentheses in {boundary_path}")

    block = txt[i0+1:i1]  # content inside parentheses

    patches: dict[str, str] = {}
    # match blocks:  name { ... type <word>; ... }
    for m in re.finditer(r"(?ms)^\s*([A-Za-z_]\w*)\s*\{(.*?)^\s*\}", block):
        name = m.group(1)
        body = m.group(2)

        mt = re.search(r"(?m)^\s*type\s+([A-Za-z_]\w*)\s*;", body)
        if mt:
            patches[name] = mt.group(1)

    return patches

# Give transform points command(s)
def _emit_transform_points(mesh, region: str) -> str:
    fl = of_flavour()
    if fl is None:
        # safest fallback: runtime detection in shell (or raise)
        raise RuntimeError("OpenFOAM environment not sourced (WM_PROJECT not set).")

    cmds = ""

    ops = getattr(mesh, "transforms", [])
    if not ops:
        return cmds

    # Detect OpenFOAM flavour at runtime
    # Done in shell, not Python, so it works inside Allrun
    for op, vec in ops:
        vec_str = f"{vec[0]} {vec[1]} {vec[2]}"

        prefix = f"runApplication -s {region} transformPoints -region {region} " if region else "runApplication -o transformPoints "

        if fl == "esi":
            cmds += prefix + f'-{op} "({vec_str})"\n'
        else:  # esi
            cmds += prefix + f'"{op}=({vec_str})"\n'

    return cmds



def run_preprocessing(
        case: Case,
        verbose: bool=False,
        isRenumberMesh: bool=False
    ):
    """
    Execute all the OpenFOAM pre-processing utilities:
    - blockMesh
    - ideasUnvToFoam
    - renumberMesh (optional)
    - topoSet
    - createPatch
    - setFields
    - createBaffles
    - checkMesh
    - decomposePar
    - GeN-Foam -initializeMappedFields

    Parameters
    ----------
    case : Case
    verbose : bool
        Add `echo` commands in `Allrun.prerun` (default `False`)
    """
    # Find OpenFOAM version
    of_version = of_flavour()

    # Prepare commands list
    commands = ". ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions\n"

    addVerbose = lambda text: f"echo {text}\n" if verbose else ""

    # Mesh generation
    for solver in case.solvers:
        region = solver.region
        # No region in the case
        if (region == ""):
            # Mesh generation
            if (isinstance(solver.mesh, BlockMesh)):
                if(of_version == "esi"):
                    commands += f"runApplication blockMesh -merge-points\n"
                elif(of_version == "foundation"):
                    commands += f"runApplication blockMesh\n"
            elif (isinstance(solver.mesh, UnvMesh)):
                commands += f"runApplication ideasUnvToFoam {solver.mesh.srcpath}\n"
                # Change patch type for special patches (user-provided)
                for face in solver.mesh.faces:
                    if(of_version == "esi"):
                        commands += (
                            f"runApplication -o foamDictionary constant/polyMesh/boundary "
                            f"-entry entry0.{face.name}.type -set {face.boundaryType}\n")
                                # set extra parameters if present (e.g. neighbourPatch/Region, owner, updateAMI, ...)

                        if getattr(face, "extraParameters", None):
                            for k, v in face.extraParameters.items():
                                commands += (
                                    f"runApplication -o foamDictionary constant/polyMesh/boundary "
                                    f"-entry entry0.{face.name}.{k} -set {v}\n"
                                )

                    elif(of_version == "foundation"):
                        commands += (
                            f"runApplication -o foamDictionary constant/polyMesh/boundary "
                            f"-entry entry0/{face.name}/type -set {face.boundaryType}\n")
                                # set extra parameters if present (e.g. neighbourPatch/Region, owner, updateAMI, ...)

                        if getattr(face, "extraParameters", None):
                            for k, v in face.extraParameters.items():
                                commands += (
                                    f"runApplication -o foamDictionary constant/polyMesh/boundary "
                                    f"-entry entry0/{face.name}/{k} -set {v}\n"
                                )

            elif (isinstance(solver.mesh, PolyMesh)):
                commands += f"runApplication -s {region} renumberMesh -region {region} -overwrite\n"

            # Mesh manipulation
            if (isRenumberMesh):
                commands += f"runApplication renumberMesh -region {region} -overwrite\n"
            if (not solver.mesh.topoSetDict.is_empty):
                commands += f"runApplication topoSet -noZero\n"
            if (not solver.mesh.createPatchDict.is_empty):
                commands += f"runApplication createPatch -overwrite\n"
            if (not solver.setFieldsDict.is_empty):
                commands += f"runApplication setFields\n"
            if (not solver.mesh.createBafflesDict.is_empty):
                commands += f"runApplication createBaffles -overwrite\n"
            if (solver.mesh.transforms):
                commands += _emit_transform_points(solver.mesh, region="")

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

            # Mesh manipulation
            if (isRenumberMesh):
                commands += f"runApplication -s {region} renumberMesh -region {region} -overwrite -no-fields\n"
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
    if (case.is_parallel):
        commands += "runApplication decomposePar -allRegions -copyUniform\n"

    # GeN-Foam initialize mapped fields only if more than 2 defined solvers
    if (
        case.settings.application == "GeN-Foam"
        and len([solver for solver in case.solvers if solver.solver != 'none']) >= 2
    ):
        commands += "GeN-Foam -initializeMappedFields > log.GeN-Foam.initializeMappedFields\n"
        commands += addVerbose("Run GeN-Foam -initializeMappedFields")

    # Change directory to case folder
    cwd = os.getcwd()
    os.chdir(case.caseFolder)

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

    # Collect new patches that have been created and might have not been known before
    region_meshes = {}
    for solver in case.solvers:
        patches_all = read_boundary_patches(Path(f'constant/{solver.region}/polyMesh/boundary'))

        faces_by_name = {f.name: f for f in solver.mesh.faces}

        for name, ptype in patches_all.items():
            if name in faces_by_name:
                # keep user intent OR sync from boundary file (choose one)
                # 1) keep user intent (recommended):
                continue

                # 2) OR sync from boundary file:
                # faces_by_name[name].boundaryType = ptype
            else:
                solver.mesh.faces.append(Face(name, boundaryType=ptype))

        region_meshes[solver.region] = solver.mesh

    # Check boundary field
    for timeFolder in case.timeFolders:
        timeFolder.export_to_openfoam(region_meshes=region_meshes)

    # Return to the original folder
    os.chdir(cwd)


def run_reconstruction(
        case: Case,
        isAllRegions: bool=True,
        isLatestTime: bool=False
    ):
    """
    Execute the `reconstructPar` OpenFOAM utility

    Parameters
    ----------
    case : Case
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
    os.chdir(case.caseFolder)

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
        case: Case,
        is_preprocessing: bool=False,
    ):
    """
    Command to run a FFN simulation.

    Parameters
    ----------
    case : Case
        FFN Case containing the solvers, time folders, coupling objects, ...
    is_preprocessing : bool
        Execute preprocessing programs such as mesh generation before running
        the main simulation (default `False`).
    """
    if (is_preprocessing):
        run_preprocessing(case)

    if (not case.is_parallel):
        app = case.settings.application
        case = foamlib.FoamCase(case.caseFolder)
        case.run(
            cmd=[app],
            check=True,
            log=True
        )
    else:
        # Need to do it manually, foamLib is using mpiexec and it is not writing
        # in log file
        app = case.settings.application
        nCores = case.get_number_processors

        # Change directory to case folder
        cwd = os.getcwd()
        os.chdir(case.caseFolder)

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
        case: Case,
        isLatestTimeStep: bool=False
    ):
    commands = ". $WM_PROJECT_DIR/bin/tools/RunFunctions\n"
    commands += f"runApplication reconstructPar -allRegions -newTimes"

    if (isLatestTimeStep):
        commands += " -latestTime"

    commands += "\n"

    # Change directory to case folder
    cwd = os.getcwd()
    os.chdir(case.caseFolder)

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


def allclean(case=None):
    """
    Clean an OpenFOAM case folder.

    If case is None:
        Clean the current working directory.

    If case is provided:
        Clean case.caseFolder.
        If the folder does NOT exist, exit without error.
    """

    # Determine target folder
    folder = case.caseFolder if case is not None else os.getcwd()

    # If folder does not exist, do nothing (safe)
    if not os.path.isdir(folder):
        return "Folder not found, nothing to clean", ""

    # Prepare cleaning commands
    commands  = ". ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions\n"
    # NOTE: Changed to cleanCase + rm -rf 0 to be compatible with Foundation v9
    commands += "cleanCase\n"
    commands += (
        "rm -rf 0\n"
        "rm -rf constant system processor* \n"
        "log.* Allrun.pre* *.png *.gif *_log.txt \n"
        "*.csv *.json input_check.txt\n"
    )

    # Run inside folder
    cwd = os.getcwd()
    os.chdir(folder)

    try:
        process = subprocess.Popen(
            "/bin/bash",
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        out, err = process.communicate(commands)
    finally:
        os.chdir(cwd)

    return out, err


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
