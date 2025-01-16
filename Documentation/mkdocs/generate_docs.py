import os
import re

# Function to extract Description and Usage from .H files
def extract_descriptions_and_usage(file_path):
    """
    Parses a .H file to extract the Description and Usage sections.
    """
    with open(file_path, 'r') as file:
        content = file.read()

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

# Function to check if a class is runTimeSelectable by looking in the corresponding .C file
def is_runTimeSelectable(file_path):
    """
    Check if the class is runTimeSelectable by searching for the 'addToRunTimeSelectionTable' macro in the .C file.
    """
    c_file_path = file_path.replace('.H', '.C')  # Replace .H extension with .C to find the implementation file

    if os.path.exists(c_file_path):
        with open(c_file_path, 'r') as file:
            content = file.read()

        # Check for the 'addToRunTimeSelectionTable' macro in the .C file
        return 'addToRunTimeSelectionTable' in content
    return False  # If the .C file doesn't exist, assume it's not runTimeSelectable

# Function to create markdown files for only runTimeSelectable classes in ../offbeatLib/ and return class entries for mkdocs.yml
def generate_class_markdown_files(offbeat_lib_dir, markdown_output_dir):
    """
    Generates markdown files for only runTimeSelectable classes in offbeatLib, excluding lnInclude folder.
    Returns class entries for mkdocs.yml.
    """
    class_entries = {}

    for root, dirs, files in os.walk(offbeat_lib_dir):
        # Skip any directory named 'lnInclude'
        dirs[:] = [d for d in dirs if (d != 'lnInclude' and d != 'Make')]

        # Preserve the relative path structure for the output directory
        relative_path = os.path.relpath(root, offbeat_lib_dir)
        output_dir = os.path.join(markdown_output_dir, relative_path)
        os.makedirs(output_dir, exist_ok=True)

        for file_name in files:
            if file_name.endswith(".H"):
                file_path = os.path.join(root, file_name)

                # Check if the corresponding .C file contains 'addToRunTimeSelectionTable'
                if is_runTimeSelectable(file_path):
                    class_name = file_name[:-2]  # Remove .H extension to get class name

                    # Extract Description and Usage
                    description, options, usage, vartable = extract_descriptions_and_usage(file_path)

                    # Transform description and usage for markdown
                    description = transform_content(description)
                    options = transform_content(options, is_options_section=True)
                    usage = transform_content(usage)
                    vartable = transform_vartable(vartable)

                    # Create the markdown content
                    markdown_content = f"# **{class_name}**\n\n"
                    markdown_content += f"## **Description**\n\n{description}\n\n"
                    if (vartable != None):
                        markdown_content += f"## **Variables**\n\n{vartable}\n\n"
                    markdown_content += f"## **Options**\n\n{options}\n\n"
                    markdown_content += f"## **Usage**\n\n{usage}\n"

                    # Write to markdown file in the corresponding output directory
                    markdown_file_path = os.path.join(output_dir, f"{class_name}.md")
                    with open(markdown_file_path, 'w') as markdown_file:
                        markdown_file.write(markdown_content)

                    # Create the entry for mkdocs.yml
                    relative_class_path = os.path.join('docs/classes', relative_path, f"{class_name}.md")
                    class_entries[class_name] = relative_class_path

    return class_entries

def transform_vartable(content):
    content = re.sub(r'^\t|^ {4}', '', content, flags=re.MULTILINE)
    content = re.sub(r'\\vartable', '', content, flags=re.MULTILINE)
    content = re.sub(r'\\endvartable', '', content, flags=re.MULTILINE)
    if ("|" not in content):
        return(None)
    res  = "| Parameter | Description |\n"
    res += "|:-|:-|\n"
    for line in content.split('\n'):
        if ("|" in line):
            line = line.split("|")
            line[0] = f"${line[0].rstrip()}$"
            res += f"| {line[0]} | {line[1]} |\n"
    return(res)

# Function to transform specific syntax in the extracted content
def transform_content(content, is_options_section=False):
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
            transformed_lines.append("???+ note\n")
            in_note_block = True
            continue

        # Check for \warning and start the warning block
        elif r"\warning" in line:
            transformed_lines.append("???+ warning\n")
            in_warning_block = True
            continue

        # Check for a table section in the Options section
        elif is_options_section and "Parameters in" in line:
            transformed_lines.append(f"???+ info \"{line.strip()}\"\n")
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

    # Transform "\verbatim" blocks into a Markdown code block with C++ syntax highlighting
    transformed_content = re.sub(
        r"\\verbatim(.+?)\\endverbatim",
        lambda m: "```cpp\n" +
                  "\n".join([line[1:] if line.startswith('\t') else line[4:] if line.startswith('    ') else line for line in m.group(1).splitlines()]) +
                  "\n```",
        transformed_content,
        flags=re.DOTALL
    )

    # Replace "\f$" with "$" for inline math
    transformed_content = re.sub(r"\\f\$", "$", transformed_content)

    # Replace "#### Formulation" with "### **Formulation**"
    transformed_content = re.sub(r"#### Formulation", "### **Formulation**", transformed_content)

    return transformed_content

# Function to read class order from the Overview file
def read_order_from_overview(overview_file_path):
    """
    Reads the class order from the overview file (like user_guide/solverDict/thermal_solver.md).
    Assumes the classes are listed as "- ClassName" in Markdown.
    """
    order = []
    if os.path.exists(overview_file_path):
        with open(overview_file_path, 'r') as file:
            for line in file:
                match = re.match(r"-\s*(\w+)", line)
                if match:
                    order.append(match.group(1))
    return order

# Function to generate the dynamic nav structure for a specific folder, sorted by the order from the overview file
def generate_dynamic_nav_for_folder(class_entries, folder_path, overview_file=None):
    """
    Generates dynamic class navigation for the specified folder, respecting the order in the overview file if provided.
    Only checks the first level of folders, non-recursively.
    """
    # Read class order from overview file, if it exists
    class_order = []
    if overview_file and os.path.exists(overview_file):
        class_order = read_order_from_overview(overview_file)

    # Filter the classes that belong to the specified folder at the first level
    filtered_class_entries = {name: path.replace('docs/', '') for name, path in class_entries.items() if os.path.dirname(path) == os.path.dirname(folder_path)}

    # Sort class entries based on the order in the overview file, if it exists
    ordered_nav = []
    if class_order:
        for class_name in class_order:
            if class_name in filtered_class_entries:
                ordered_nav.append(f"      - {class_name}: {filtered_class_entries[class_name]}")

    # Add any classes that were not listed in the overview file at the end, or list all if no overview file
    unordered_classes = [f"      - {name}: {path}" for name, path in filtered_class_entries.items() if name not in class_order]
    ordered_nav.extend(unordered_classes)

    return "\n".join(ordered_nav)

# Function to scan release notes folder and generate entries for mkdocs.yml
def generate_release_notes_entries(release_notes_folder):
    """
    Scans the release_notes folder and generates entries for mkdocs.yml.
    Also updates the release_notes_index.md file.
    """
    release_notes_entries = []
    index_entries = []

    # Sort release notes by version, assuming files are named like 'v_X.X.md'
    for file_name in sorted(os.listdir(release_notes_folder)):
        if file_name.startswith("v") and file_name.endswith(".md"):
            version = file_name.replace(".md", "").replace("_",".")
            # Add the release note entry to mkdocs.yml
            release_notes_entries.append(f"    - Release Note {version}: release_notes/{file_name}")
            # Add the release note entry to the release_notes_index.md
            index_entries.append(f"- [Release Note {version}]({file_name})")

    return release_notes_entries, index_entries

# Function to update the release_notes_index.md file with all release note links
def update_release_notes_index(index_entries, index_file_path):
    """
    Writes the updated release notes entries to the release_notes_index.md file.
    """
    with open(index_file_path, 'w') as index_file:
        index_file.write("# <b>Release Notes Index</b>\n\n")
        index_file.write("\n".join(index_entries))
        index_file.write("\n")

# Functio to indent a nav block
def indent_navigation(nav_block, level):
    spaces = ' ' * level
    indented_nav = "\n".join([f"{spaces}{line}" for line in nav_block.splitlines()])
    return indented_nav

# Main function to generate both markdown files and mkdocs.yml
def main():
    # Paths to the necessary directories
    offbeat_lib_dir = "../../src/classes"  # Path to the OFFBEAT classes
    markdown_output_dir = "docs/classes"  # Path to output markdown files
    # release_notes_folder = "docs/release_notes"  # Path to release notes
    # release_notes_index_file = "docs/release_notes/index.md"  # Path to the release notes index file

    # Generate release notes entries for mkdocs.yml and the release_notes_index.md file
    # release_notes_entries, index_entries = generate_release_notes_entries(release_notes_folder)

    # Update the release_notes_index.md file
    # update_release_notes_index(index_entries, release_notes_index_file)

    # Step 1: Generate markdown files for all classes in offbeatLib
    os.makedirs(markdown_output_dir, exist_ok=True)
    class_entries = generate_class_markdown_files(offbeat_lib_dir, markdown_output_dir)

    # Step 2: Generate class entries for dynamic folders, respecting order from Overview files
    # Define the dictionary with folders and corresponding overview files
    nav_sections = {
        # Neutronics
        "diffusion": ("docs/classes/neutronics/diffusion/", "docs/user_guide/neutronics.md"),
        "adjoint_diffusion": ("docs/classes/neutronics/adjointDiffusion/", "docs/user_guide/neutronics.md"),
        "SP3": ("docs/classes/neutronics/SP3/", "docs/user_guide/neutronics.md"),
        "SN": ("docs/classes/neutronics/SN/", "docs/user_guide/neutronics.md"),
        "point_kinetics": ("docs/classes/neutronics/pointKinetics/", "docs/user_guide/neutronics.md"),
        # Thermal-hydraulics
        "one_phase": ("docs/classes/thermalHydraulics/solvers/onePhase/", "docs/user_guide/thermalHydraulicsSolvers.md"),
        "two_phase": ("docs/classes/thermalHydraulics/solvers/twoPhase/", "docs/user_guide/thermalHydraulicsSolvers.md"),
        "compressible_inter_foam": ("docs/classes/openFoamImportedSolvers/compressibleInterFoam/", "docs/user_guide/thermalHydraulicsSolvers.md"),
        # Thermal-hydraulics sub-scale structures
        "fixed_power": ("docs/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedPower/", "docs/user_guide/thermalHydraulicsSubScaleStructures.md"),
        "fixed_temperature": ("docs/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/fixedTemperature/", "docs/user_guide/thermalHydraulicsSubScaleStructures.md"),
        "heated_pin": ("docs/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/heatedPin/", "docs/user_guide/thermalHydraulicsSubScaleStructures.md"),
        "lumped_nuclear_structure": ("docs/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/lumpedNuclearStructure/", "docs/user_guide/thermalHydraulicsSubScaleStructures.md"),
        "nuclear_fuel_fmu": ("docs/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearFuelFMU/", "docs/user_guide/thermalHydraulicsSubScaleStructures.md"),
        "nuclear_fuel_pin": ("docs/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearFuelPin/", "docs/user_guide/thermalHydraulicsSubScaleStructures.md"),
        "nuclear_steady_state_pebble": ("docs/classes/thermalHydraulics/src/phaseModels/structureModels/powerModels/nuclearSteadyStatePebble/", "docs/user_guide/thermalHydraulicsSubScaleStructures.md"),
    }

    # Generate dynamic navigation for all sections
    dynamic_navs = {}
    for key, (folder, overview) in nav_sections.items():
        dynamic_navs[key] = generate_dynamic_nav_for_folder(class_entries, folder, overview)

    # Indent boundary condition navigation if necessary
    # for key in ["fixed_temperature_bcs", "convective_htc_bcs", "temperature_coupled_bcs", "radiative_convective_sink_bcs",
    #             "traction_displacement_bcs", "implicit_contact_bcs", "fixed_displacement_bcs", "fixed_displacement_zero_shear_bcs"]:
    #     dynamic_navs[key] = indent_navigation(dynamic_navs[key], 3)
    for key in [
        "one_phase", "two_phase", "compressible_inter_foam",
        "fixed_power", "fixed_temperature", "heated_pin", "lumped_nuclear_structure", "nuclear_fuel_fmu", "nuclear_fuel_pin", "nuclear_steady_state_pebble"
    ]:
        dynamic_navs[key] = indent_navigation(dynamic_navs[key], 2)

    # Indent nav blocks if necessary
    # dynamic_navs["multi_material_correction"] = indent_navigation(dynamic_navs["multi_material_correction"], 2)

    # Predefined structure of mkdocs.yml with placeholders for dynamic class nav inserted
    mkdocs_yml_template = """site_name: GeN-Foam Documentation
nav:
  - Introduction: index.md
  - Installation: installation.md
  - User's Guide:
    - user_guide/index.md
    - Neutronics:
      - Overview: user_guide/neutronics.md
{dynamic_diffusion_nav}
{dynamic_adjoint_diffusion_nav}
{dynamic_SP3_nav}
{dynamic_SN_nav}
{dynamic_point_kinetics_nav}
    - Thermal-hydraulics:
      - Overview: user_guide/thermalHydraulics.md
      - Solvers:
        - Overview: user_guide/thermalHydraulicsSolvers.md
{dynamic_one_phase_nav}
{dynamic_two_phase_nav}
{dynamic_compressible_inter_foam_nav}
      - Sub-scale structures:
        - Overview: user_guide/thermalHydraulicsSubScaleStructures.md
{dynamic_fixed_power_nav}
{dynamic_fixed_temperature_nav}
{dynamic_heated_pin_nav}
{dynamic_lumped_nuclear_structure_nav}
{dynamic_nuclear_fuel_fmu_nav}
{dynamic_nuclear_fuel_pin_nav}
{dynamic_nuclear_steady_state_pebble_nav}
    - Thermo-mechanics: user_guide/thermoMechanics.md
    - Coupling solvers and time stepping: user_guide/coupling.md
    - FMU: user_guide/fmu.md
    - Pre-processing: user_guide/pre_processing.md
    - Running GeN-Foam: user_guide/running_GeN-Foam.md
    - Post-processing: user_guide/post_processing.md
    - Important notes: user_guide/important_notes.md
    - Tips and tricks: user_guide/tips_and_tricks.md
  - Tools: tools/index.md
  - Tutorials: tutorials/index.md
  - Developer's Guide: developer_guide/index.md
  - Publications: publications.md
  - Bibliography: bibliography.md
  - Contributors: contributors.md

theme:
  name: material
  favicon: logo/GeN-Foam-logo.png
  logo: logo/GeN-Foam-logo.png
  icon:
    admonition:
      note: octicons/tag-16
      abstract: octicons/checklist-16
      info: octicons/info-16
      tip: octicons/squirrel-16
      success: octicons/check-16
      question: octicons/question-16
      warning: octicons/alert-16
      failure: octicons/x-circle-16
      danger: octicons/zap-16
      bug: octicons/bug-16
      example: octicons/beaker-16
      quote: octicons/quote-16
  features:
    - content.code.annotate
    - navigation.indexes
    - header.autohide
  palette:

    # Palette toggle for light mode
    - scheme: default
      toggle:
        icon: material/brightness-7
        name: Switch to dark mode

    # Palette toggle for dark mode
    - scheme: slate
      toggle:
        icon: material/brightness-4
        name: Switch to light mode

markdown_extensions:
  - footnotes
  - admonition
  - pymdownx.details
  - pymdownx.superfences:
      custom_fences:
        - name: mermaid
          class: mermaid
          format: !!python/name:pymdownx.superfences.fence_code_format
  - attr_list
  - pymdownx.highlight:
      anchor_linenums: true
      line_spans: __span
      pygments_lang_class: true
  - pymdownx.inlinehilite
  - pymdownx.snippets
  - pymdownx.superfences

  - pymdownx.emoji:
      emoji_index: !!python/name:material.extensions.emoji.twemoji
      emoji_generator: !!python/name:material.extensions.emoji.to_svg
  - md_in_html
  - pymdownx.snippets
  - pymdownx.arithmatex:
      generic: true

extra_javascript:
  - javascripts/mathjax.js
  - https://unpkg.com/mathjax@3/es5/tex-mml-chtml.js

extra_css:
  - css/extra.css

plugins:
  - search
  - bibtex:
      bib_file: "docs/references.bib"
"""
    # Insert the generated release notes entries into mkdocs.yml
    # release_notes_section = "\n".join(release_notes_entries)

    # Replace placeholders in mkdocs.yml
    mkdocs_yml = mkdocs_yml_template.format(
        # release_notes_section=release_notes_section,
        **{f"dynamic_{key}_nav": nav for key, nav in dynamic_navs.items()}
    )

    # Write the mkdocs.yml file
    with open("mkdocs.yml", 'w') as mkdocs_file:
        mkdocs_file.write(mkdocs_yml)

    print("mkdocs.yml generated successfully.")

if __name__ == "__main__":
    main()
