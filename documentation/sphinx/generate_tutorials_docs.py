#==============================================================================*
# Imports

import os
from re import finditer
import shutil
import sys


#==============================================================================*
# Functions

def camel_case_split(identifier):
    matches = finditer('.+?(?:(?<=[a-z])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])|$)', identifier)
    return(' '.join([m.group(0).capitalize() for m in matches]))


def getTitleFromMarkdown(text: str) -> str:
    for line in text.split('\n'):
        if (line[:2] == '# '):

            return(line[2:])

    return('')


def getTagsFromMarkdown(text: str, tutorialName: str) -> list[str]:
    for line in text.split('\n'):
        if ("Tags: " == line[:6]):
            tags = line[6:].replace('[![badge](', '').replace('.svg)]()', '').split()

            tags = [
                (
                    f"|{tutorialName}-{'-'.join(tag.split('badge/')[1].split('-')[:2])}-badge|",
                    f".. |{tutorialName}-{'-'.join(tag.split('badge/')[1].split('-')[:2])}-badge| image:: {tag}"
                )
                for tag in tags
            ]

            return(tags)

    return([])

#==============================================================================*


# Get the absolute path to the directory containing this script
script_dir = os.path.dirname(os.path.abspath(__file__))


# Construct the correct path to the tutorials directory
tutorials_dir = os.path.abspath("tutorials")

# Destination base directory for documentation
docs_base_dir =  os.path.join(script_dir, "usersguide")
docs_tutorials_dir = os.path.join(docs_base_dir, "tutorials")

# Ensure the tutorials documentation directory exists
os.makedirs(docs_tutorials_dir, exist_ok=True)

text_list_tutorials = """
List of tutorials
=================

.. toctree::
    :maxdepth: 2

"""

text_table_tutorials = """
Tutorials
=========

.. raw:: html

   <input type="text" class="tableFilter" placeholder="Search tutorials...">

.. list-table:: Tutorials
    :widths: 50 50 20
    :header-rows: 1
    :class: wy-table-responsive filterable-table

    * - Tutorial title
      - Tags and categories
      - Link to description
"""

text_tags_tutorials = ""


# Iterate over section folders in the tutorials source
for section in sorted(os.listdir(tutorials_dir)):
    section_path = os.path.join(tutorials_dir, section)
    if os.path.isdir(section_path):
        # Create corresponding section folder in documentation
        docs_section_path = os.path.join(docs_tutorials_dir, section)
        os.makedirs(docs_section_path, exist_ok=True)

        text_list_tutorials += f"    tutorials/{section}/index\n"

        # Create section index.rst in documentation
        section_text  = f"{camel_case_split(section)} Tutorials\n"
        section_text += f"{'=' * (len(section) + 11)}\n\n"
        section_text += ".. toctree::\n"
        section_text += "    :maxdepth: 1\n\n"

        # Iterate over tutorials in the section
        for tutorial in sorted(os.listdir(section_path)):
            tutorial_path = os.path.join(section_path, tutorial)

            if os.path.isdir(tutorial_path):
                readme_path = os.path.join(tutorial_path, "README.md")
                tutorial_rst_path = os.path.join(docs_section_path, f"{tutorial}.md")

                # Create tutorial .rst file in documentation
                with open(tutorial_rst_path, "w") as tutorial_rst:
                    # tutorial_rst.write(f"{tutorial}\n{'=' * len(tutorial)}\n\n")
                    if os.path.exists(readme_path):
                        with open(readme_path, "r") as readme_file:
                            readme_file_content = readme_file.read()

                        tutorial_rst.write(readme_file_content)

                        tutorial_name = "tutorials" + readme_path.split('tutorials')[-1].replace('/README.md', '')
                        tutorial_tags = getTagsFromMarkdown(readme_file_content, tutorial_name)
                        tutorial_title = getTitleFromMarkdown(readme_file_content)

                        tutorial_tags_label = '| ' + '\n        | '.join([label for label, _ in tutorial_tags])
                        text_tags_tutorials += '\n\n'.join([tag for _, tag in tutorial_tags]) + '\n\n'

                        text_table_tutorials += f"    * - {tutorial_title}\n"
                        text_table_tutorials += f"      - {tutorial_tags_label}\n"
                        text_table_tutorials += f"      - :doc:`Link <{tutorial_name}>`\n"

                # Add tutorial to section index
                section_text += f"    {tutorial}\n"

                # Add missing images
                tutorial_image_path = os.path.join(tutorial_path, "images")
                if os.path.isdir(tutorial_image_path):
                    doc_image_path = os.path.join(docs_section_path, "images")
                    os.makedirs(doc_image_path, exist_ok=True)

                    for image in os.listdir(tutorial_image_path):
                        image_path = os.path.join(tutorial_image_path, image)
                        shutil.copyfile(image_path, os.path.join(doc_image_path, image))


        # Write section index per tutorials class
        section_index_path = os.path.join(docs_section_path, "index.rst")
        with open(section_index_path, "w") as section_index:
            section_index.write(section_text)


#==============================================================================*

# Create the main tutorials.rst file
main_tutorials_path = os.path.join(docs_base_dir, "tutorials.rst")

import_readme_path = os.path.join(tutorials_dir,  "README.rst")
if os.path.exists(import_readme_path):
    with open(import_readme_path, "r") as readme_file:
        readme_file_content = readme_file.read()

main_text = f""".. _usersguide_tutorials:

{readme_file_content}

{text_table_tutorials}

{text_tags_tutorials}

"""
# {text_list_tutorials}

with open(main_tutorials_path, "w") as main_index:
    main_index.write(main_text)


#==============================================================================*

print(f"python3 {sys.argv[0]} ... End")


#==============================================================================*
