"""
Generate the doc for Sphinx RST files
"""

#==============================================================================*
# Imports

import os
import re
import sys

#==============================================================================*

def replaceInlineReference(text):
    """
    Replace "[@REF]" for inline reference
    """
    if ("[@" not in text):
        return(text)

    refValues = []
    for line in text.split("\n"):
        if ("[@" in line):
            lineSplit = line.split("[@")
            for section in lineSplit[1:]:
                refValues.append(section.split("]")[0])

    for refValue in refValues:
        text = text.replace(f"[@{refValue}]", f":ref:`{refValue} <{refValue}>`")

    return(text)


def replaceInlineMath(text):
    """
    Replace "\f$" and "$" for inline math
    """
    # Remove \f$
    mathInlineFlag = True
    while (r"\f$" in text):
        text = re.sub(
            r"\\f\$\s*" if mathInlineFlag else r"\s*\\f\$",
            ":math:`" if mathInlineFlag else "`",
            text,
            1
        )
        mathInlineFlag = not mathInlineFlag

    # Remove $
    refValues = []
    for line in text.split("\n"):
        if ("$" in line):
            lineSplit = line.split("$")
            if (len(lineSplit) % 2 == 0):
                continue
            for section in lineSplit[1::2]:
                refValues.append(section)

    for refValue in refValues:
        text = text.replace(f"${refValue}$", f":math:`{refValue}`")

    return(text)


def format_tables_for_rst(text: str) -> str:
    tableRow = []
    textSplit = text.split('\n')
    for i, line in enumerate(textSplit):
        if (r"\table" in line):
            nCols = len(textSplit[i+1].split('|'))
            tableRow.append((i+2+len(tableRow), nCols))

    for idx, nCols in tableRow:
        textSplit = textSplit[:idx] + [(nCols-1)*'-|' + '-'] + textSplit[idx:]

    new_text = '\n'.join(textSplit)
    new_text = new_text.replace(r'\table', '')
    new_text = new_text.replace(r'\endtable', '')

    return(new_text)


def format_table_for_rst(text: str) -> str:
    tableRow = []
    textSplit = text.split('\n')
    isRead = False
    new_text = text
    for i, line in enumerate(textSplit):
        if (r"\table" in line):
            nCols = textSplit[i+1].count('|')+1
            textSplit[i] = textSplit[i].replace(r'\table', f"\n.. list-table::\n    :widths:{nCols*' 50'}\n    :header-rows: 1\n", 1)

            isRead = True
            continue

        if (r"\endtable" in line):
            isRead = False

            for idx in tableRow:
                strippedLine = [element.strip() for element in textSplit[idx].strip().split("|")]
                textSplit[idx] = "    * - " + '\n      - '.join(strippedLine)

            textSplit[i] = textSplit[i].replace(r'\endtable', '', 1)
            new_text = '\n'.join(textSplit)

            tableRow = []
            continue

        if (isRead):
            tableRow.append(i)

    return(new_text)


def format_equation_for_rst(text: str) -> str:
    tableRow = []
    textSplit = text.split('\n')
    isRead = False
    new_text = text
    for i, line in enumerate(textSplit):
        if ((r"\f[" in line or r"$$" in line) and not isRead):
            textSplit[i] = ".. math::"
            isRead = True
            continue
        if ((r"\f]" in line or r"$$" in line) and isRead):
            isRead = False

            for idx in tableRow:
                textSplit[idx] = f"    {textSplit[idx].strip()}"

            textSplit[i] = ""
            new_text = '\n'.join(textSplit)

            tableRow = []
            continue

        if (isRead):
            tableRow.append(i)

    return(new_text)


def extract_descriptions_and_usage(file_path):
    """
    Function to extract Description and Usage from .H files

    Parses a .H file to extract the Description and Usage sections.
    """
    try:
        with open(file_path, 'r') as file:
            content = file.read()
    except:
        return("", "", "", "")

    description_pattern = re.compile(r"Description\s*\n\s*(.*?)(?=\n\s*(Usage|\\par Options|SourceFiles|Class|\\\*-----|\*/|\\vartable))", re.DOTALL)
    options_pattern = re.compile(r"\\par Options\s*\n\s*(.*?)(?=\n\s*(Usage|SourceFiles|Class|\\\*-----|\*/))", re.DOTALL)
    usage_pattern = re.compile(r"Usage\s*\n\s*(.*?)(?=\n\s*(SourceFiles|Class|\\\*-----|\*/))", re.DOTALL)
    vartable_pattern = re.compile(r"\\vartable.*?\\endvartable", re.DOTALL)

    description = description_pattern.search(content)
    options = options_pattern.search(content)
    usage = usage_pattern.search(content)
    vartable = vartable_pattern.search(content)

    description_text = description.group(1).strip() if description else "No description available."
    options_text = options.group(1).strip() if options else "No options available."
    usage_text = usage.group(1).strip() if usage else "No usage available."
    vartable_text = vartable.group(0).strip() if vartable else "No vartable available."

    return description_text, options_text, usage_text, vartable_text


def is_runTimeSelectable(file_path):
    """
    Function to check if a class is runTimeSelectable by looking in the
    corresponding .C file

    Check if the class is runTimeSelectable by searching for the
    'addToRunTimeSelectionTable' macro in the .C file.
    """
    c_file_path = file_path.replace('.H', '.C')  # Replace .H extension with .C to find the implementation file

    if os.path.exists(c_file_path):
        with open(c_file_path, 'r') as file:
            content = file.read()

        # Check for the 'addToRunTimeSelectionTable' macro in the .C file
        return 'addToRunTimeSelectionTable' in content
    return False  # If the .C file doesn't exist, assume it's not runTimeSelectable


def generate_class_rst_files(ffn_lib_dirs: list[str], rst_output_dir: str) -> dict:
    """
    Generates rst files for only runTimeSelectable classes in ffn, excluding lnInclude folder.
    Returns class entries for mkdocs.yml.
    """
    class_entries = {}

    for ffn_lib_dir in ffn_lib_dirs:
        for root, dirs, files in os.walk(ffn_lib_dir):
            # Skip any directory named 'lnInclude'
            dirs[:] = [d for d in dirs if (d != 'lnInclude' and d != 'Make')]

            # Preserve the relative path structure for the output directory
            relative_path = os.path.relpath(root, ffn_lib_dir)
            output_dir = os.path.join(rst_output_dir, relative_path)
            os.makedirs(output_dir, exist_ok=True)

            for file_name in files:
                if file_name.endswith(".H"):
                    file_path = os.path.join(root, file_name)

                    # Check if the corresponding .C file contains 'addToRunTimeSelectionTable'
                    # if is_runTimeSelectable(file_path):
                    class_name = file_name[:-2]  # Remove .H extension to get class name

                    # Extract Description and Usage
                    description, options, usage, vartable = extract_descriptions_and_usage(file_path)

                    # Transform description and usage for rst
                    description = transform_content(description)
                    options = transform_content(options, is_options_section=True)
                    usage = transform_content(usage)
                    vartable = transform_vartable(vartable)

                    # Create the content
                    rst_content = ""
                    rst_content += f".. _{class_name}:\n\n"

                    rst_content += f"{'':=<{len(class_name)}}\n"
                    rst_content += f"{class_name}\n"
                    rst_content += f"{'':=<{len(class_name)}}\n\n"

                    rst_content += f"Description\n"
                    rst_content += f"===========\n\n"
                    rst_content += f"{description}\n\n"

                    if (vartable != None):
                        rst_content += f"Variables\n"
                        rst_content += f"=========\n\n"
                        rst_content += f"{vartable}\n\n"

                    rst_content += f"Options\n"
                    rst_content += f"=======\n\n"
                    rst_content += f"{options}\n\n"

                    rst_content += f"Usage\n"
                    rst_content += f"=====\n\n"
                    rst_content += f"{usage}\n\n"

                    rst_content += f"Link to code\n"
                    rst_content += f"============\n\n"
                    # rst_content += f"- `{class_name}.H <https://gitlab.com/foamForNuclear/foamForNuclear/-/blob/main/{file_path}>`_\n"
                    # rst_content += f"- `{class_name}.C <https://gitlab.com/foamForNuclear/foamForNuclear/-/blob/main/{file_path.replace('.H', '.C')}>`_\n"
                    rst_content += f"- `Doxygen doc <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H.html>`_\n"
                    rst_content += f"- `{class_name}.H <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8H_source.html>`_\n"
                    rst_content += f"- `{class_name}.C <https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/{class_name}_8C_source.html>`_\n"

                    # Write to rst file in the corresponding output directory
                    rst_file_path = os.path.join(output_dir, f"{class_name}.rst")
                    with open(rst_file_path, 'w') as rst_file:
                        rst_file.write(rst_content)

                    # Create the entry for mkdocs.yml
                    relative_class_path = os.path.join(rst_output_dir, relative_path, f"{class_name}.rst")
                    class_entries[class_name] = relative_class_path

    return class_entries


def transform_vartable(content):
    content = re.sub(r'^\t|^ {4}', '', content, flags=re.MULTILINE)
    content = re.sub(r'\\vartable', '', content, flags=re.MULTILINE)
    content = re.sub(r'\\endvartable', '', content, flags=re.MULTILINE)
    if ("|" not in content):
        return(None)


    res  = ".. list-table:: Parameters\n"
    res += "    :widths: 50 50\n"
    res += "    :header-rows: 1\n\n"
    res += "    * - Parameter\n"
    res += "      - Description\n"
    for line in content.split('\n'):
        if ("|" in line):
            line = line.split("|")
            res += f"    * - :math:`{line[0].strip()}`\n"
            res += f"      - {line[1]}\n"

    res = replaceInlineMath(res)
    return(res)


def transform_content(content, is_options_section=False):
    """
    Function to transform specific syntax in the extracted content
    """
    # De-indent the entire content by removing one level of indentation at the start of each line
    content = re.sub(r'^\t|^ {4}', '', content, flags=re.MULTILINE)

    lines = content.splitlines()
    transformed_lines = []
    in_note_block = False
    in_warning_block = False
    in_info_block = False
    in_table = False

    for i, line in enumerate(lines):
        # Check for \note and start the note block
        if r"\note" in line:
            transformed_lines.append(".. note::\n\n")
            in_note_block = True
            continue

        # Check for \warning and start the warning block
        elif r"\warning" in line:
            transformed_lines.append(".. warning::\n\n")
            in_warning_block = True
            continue

        # Check for a table section in the Options section
        elif is_options_section and "Parameters in" in line:
            transformed_lines.append(f".. info:: \"{line.strip()}\"\n")
            transformed_lines.append("\t<table>")
            in_info_block = True
            in_info_block_first = True
            continue

        # Handle lines inside a note block
        if in_note_block:
            if line.strip() == "":  # End of note block on an empty line
                in_note_block = False
            else:
                transformed_lines.append(f"\t{line.strip()}")
                continue

        # Handle lines inside a warning block
        if in_warning_block:
            if line.strip() == "":  # End of warning block on an empty line
                in_warning_block = False
            else:
                transformed_lines.append(f"\t{line.strip()}")
                continue

        # Handle lines inside an info block for tables
        if in_info_block:
            # End the table on an empty line
            if (line.strip() == "" and in_info_block_first == False) or (i == len(lines)-1):
                transformed_lines.append(f"\t</table>\n")
                in_info_block = False
                continue

            # Convert bullet points in the "Parameters in" section to table rows
            match = re.match(r"- <b>`(.*?)`</b> - (.+)", line.strip())
            if match:
                if(in_info_block_first == False):
                    n_lines = len(transformed_lines)
                    transformed_lines[n_lines-1] += "</td>\n\t\t</tr>"
                in_info_block_first = False
                param_name = match.group(1)
                description = match.group(2).strip()
                transformed_lines.append(f"\t\t<tr>\n\t\t\t<td><strong>{param_name}</strong></td>\n\t\t\t<td>{description}")
            else:
                if(line.strip() == ""):
                    continue
                # Normal lines within the table are still included as-is
                n_lines = len(transformed_lines)
                transformed_lines[n_lines-1] += (f" {line.strip()}")
            continue

        # Normal lines outside of note/warning/info blocks
        transformed_lines.append(line)

    # Re-join all lines into a single content string
    transformed_content = "\n".join(transformed_lines)

    # Transform "\verbatim" blocks into a rst code block with C++ syntax highlighting
    transformed_content = re.sub(
        r"\\verbatim(.+?)\\endverbatim",
        lambda m: "\n.. code :: cpp\n" +
                #   "\n".join([line[1:] if line.startswith('\t') else line[4:] if line.startswith('    ') else line for line in m.group(1).splitlines()]) +
                  "\n    ".join([line for line in m.group(1).splitlines()]) +
                  "\n",
        transformed_content,
        flags=re.DOTALL
    )

    transformed_content = format_table_for_rst(transformed_content)
    transformed_content = format_equation_for_rst(transformed_content)

    transformed_content = replaceInlineMath(transformed_content)
    transformed_content = replaceInlineReference(transformed_content)

    return transformed_content


def generate_cppapi_index(
        cppapi_folder: str,
        class_entries: dict,
        sections: list[tuple]
    ) -> None:
    """
    Parameters
    ----------
    sections: list[tuple]
        tuple = (Key word path folder in generated folder, Section name, Section Depth, Add content)
    """

    filtered_class_entries = {
        key: item.replace(f"{cppapi_folder}/", "").replace(".rst", "")
        for key, item in class_entries.items()
        if "include" not in item
    }

    cppapi_index = ".. _cppapi:\n\n"

    cppapi_index += "===============\n"
    cppapi_index += "C++ Source Code\n"
    cppapi_index += "===============\n\n"

    for sectionKey, sectionName, sectionDepth, isAddContent in sections:
        if (sectionDepth == 0):
            cppapi_index += f"{'':-<{len(sectionName)}}\n"

        cppapi_index += f"{sectionName}\n"

        if (sectionDepth == 0):
            cppapi_index += f"{'':-<{len(sectionName)}}\n\n"
        if (sectionDepth == 1):
            cppapi_index += f"{'':-<{len(sectionName)}}\n\n"
        if (sectionDepth == 2):
            cppapi_index += f"{'':^<{len(sectionName)}}\n\n"

        if (isAddContent):
            # Add totree
            cppapi_index += ".. toctree::\n"
            cppapi_index += "   :hidden:\n"
            cppapi_index += "   :maxdepth: 1\n\n"
            for filename, filepath in filtered_class_entries.items():
                if (sectionKey in filepath):
                    cppapi_index += f"   {filepath}\n"

            cppapi_index += "\n"

            # Add classes in a table
            cppapi_index += ".. list-table::\n"
            cppapi_index += "    :widths: 50 50\n\n"
            for filename, filepath in filtered_class_entries.items():
                if (sectionKey in filepath):
                    cppapi_index += f"    * - :doc:`{filepath}`\n"

                    # Extract first line of hearder file
                    with open(class_entries[filename], 'r') as f:
                        for i, line in enumerate(f):
                            if (
                                i < 5
                                or line == "\n"
                                or "===" in line
                                or "Description" in line
                                or ".. _" in line
                            ):
                                continue

                            if ("Options" in line):
                                description = ""
                                break

                            description = line.replace("\n", "")
                            if ('.' in line):
                                description = line.split('.')[0]

                            break

                    cppapi_index += f"      - {description}\n"

        cppapi_index += "\n\n"

    # Write the file on disk
    with open(f"{cppapi_folder}/index.rst", 'w') as f:
        f.write(cppapi_index)


#==============================================================================*
# Main function to generate both rst files and index.rst

def main():
    # Paths to the necessary directories
    rst_output_dir = "documentation/sphinx/cppapi"  # Path to output rst files

    # Step 1: Generate rst files for all classes in offbeatLib
    print("Call generate_class_rst_files ...")
    os.makedirs(rst_output_dir, exist_ok=True)
    class_entries = generate_class_rst_files(
        [
            "src",
            "applications",
        ],
        rst_output_dir+"/generated"
    )

    # Step 2: Generate the main cppapi index
    print("Call generate_cppapi_index ...")
    generate_cppapi_index(
        cppapi_folder=rst_output_dir,
        class_entries=class_entries,
        sections=[
            # Key word, Section name, Section Depth, Add content
            ("modules/neutronics", "Neutronics", 0, True),
            ("modules/thermalHydraulics", "Thermal-hydraulics", 0, False),
            ("modules/thermalHydraulics/onePhase", "1-phase Solvers", 1, True),
            ("modules/thermalHydraulics/twoPhase", "2-phase Solvers", 1, True),
            ("porousMediaModels/phaseModels/structureModels", "Phase Models", 1, False),
            ("porousMediaModels/phaseModels/structureModels/powerModels", "Power Models", 2, True),
            ("porousMediaModels/phaseModels/structureModels/heatExchanger", "Heat Exchanger", 2, True),
            ("porousMediaModels/phaseModels/structureModels/pump", "Pump", 2, True),
            ("porousMediaModels/phaseModels/structureModels/powerOffCriterionModels", "Power Off Criterion", 2, True),
            ("porousMediaModels/phaseModels/structureModels", "Physics Models", 1, False),
            ("porousMediaModels/physicsModels/dragModels", "Drag Models", 2, True),
            ("porousMediaModels/physicsModels/heatTransferModels", "Heat transfer Models", 2, True),
            ("porousMediaModels/physicsModels/regimeMapModels", "Regime Map Models", 2, True),
            ("porousMediaModels/physicsModels/contactPartitionModels", "Contact Partition Models", 2, True),
            ("porousMediaModels/physicsModels/dispersionModels", "Dispersion Models", 2, True),
            ("porousMediaModels/physicsModels/fluidDiameterModels", "Fluid Diameter Models", 2, True),
            ("porousMediaModels/physicsModels/interfacialAreaModels", "Interfacial Area Models", 2, True),
            ("porousMediaModels/physicsModels/phaseChangeModels", "Phase Change Models", 2, True),
            ("porousMediaModels/physicsModels/phasePairs", "Phase Pairs", 2, True),
            ("porousMediaModels/physicsModels/templatedModels", "Templated Models", 2, True),
            ("porousMediaModels/physicsModels/turbulenceModels", "Turbulence Models", 2, True),
            ("porousMediaModels/physicsModels/Virtual Mass Models", "Virutal Mass Models", 2, True),
            ("offbeatLib", "Structural mechanics and fuel performance", 0, False),
            ("offbeatLib/physicsSubSolvers", "Physics Sub-Solvers", 1, True),
            ("offbeatLib/heatSource", "Heat Sources", 1, True),
            ("offbeatLib/fastFlux", "Fast Flux", 1, True),
            ("offbeatLib/materials/materialModel", "Materials", 1, True),
            ("offbeatLib/gapGasModel", "Gap Gas Models", 1, True),
            ("functionObjects", "Function Objects", 0, True),
            ("fvPatchFields", "Boundary Conditions", 0, False),
            ("fvPatchFields/neutronics", "Neutronics", 1, True),
            ("fvPatchFields/thermalHydraulics", "Thermal-hydraulics", 1, True),
            ("fvPatchFields/thermoMechanics", "Thermo-mechanics", 1, True),
            ("offbeatLib/fvPatchFields", "Thermo-mechanics", 1, True),
        ]
    )


    print("End")


#==============================================================================*
# Main

if __name__ == "__main__":
    main()

    print(f"python3 {sys.argv[0]} ... End")


#==============================================================================*
