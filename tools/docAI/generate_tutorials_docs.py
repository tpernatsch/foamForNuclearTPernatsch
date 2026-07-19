#==============================================================================*
# Imports

import filecmp
import os
from re import finditer
import shutil
import sys


#==============================================================================*
# Functions

def write_if_changed(path: str, text: str) -> None:
    """
    Write text only if it differs from the current file content, so that
    unchanged files keep their mtime and Sphinx skips re-reading them.
    """
    if os.path.exists(path):
        with open(path, "r") as f:
            if f.read() == text:
                return
    with open(path, "w") as f:
        f.write(text)


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


# Construct the correct path to the tutorials directory
tutorials_dir = os.path.abspath("tutorials")

# Destination base directory for documentation (relative to repo root)
docs_base_dir = os.path.abspath("documentation/sphinx/usersguide")
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

# Track generated paths so stale pages (removed tutorials/sections) can be cleaned
expected_sections: set[str] = set()
expected_pages: set[str] = set()


# Iterate over section folders in the tutorials source
for section in sorted(os.listdir(tutorials_dir)):
    section_path = os.path.join(tutorials_dir, section)
    if os.path.isdir(section_path):
        # Create corresponding section folder in documentation
        docs_section_path = os.path.join(docs_tutorials_dir, section)
        os.makedirs(docs_section_path, exist_ok=True)
        expected_sections.add(os.path.abspath(docs_section_path))

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
                expected_pages.add(os.path.abspath(tutorial_rst_path))

                # Create tutorial .md file in documentation
                if os.path.exists(readme_path):
                    with open(readme_path, "r") as readme_file:
                        readme_file_content = readme_file.read()

                    write_if_changed(tutorial_rst_path, readme_file_content)

                    tutorial_name = "tutorials" + readme_path.split('tutorials')[-1].replace('/README.md', '')
                    tutorial_tags = getTagsFromMarkdown(readme_file_content, tutorial_name)
                    tutorial_title = getTitleFromMarkdown(readme_file_content)

                    tutorial_tags_label = '| ' + '\n        | '.join([label for label, _ in tutorial_tags])
                    text_tags_tutorials += '\n\n'.join([tag for _, tag in tutorial_tags]) + '\n\n'

                    text_table_tutorials += f"    * - {tutorial_title}\n"
                    text_table_tutorials += f"      - {tutorial_tags_label}\n"
                    text_table_tutorials += f"      - :doc:`Link <{tutorial_name}>`\n"
                else:
                    # No README: write a stub so the page still has a
                    # title (an empty page breaks the toctree link).
                    write_if_changed(
                        tutorial_rst_path,
                        f"# {tutorial}\n\n"
                        "**This tutorial does not have a description yet.**\n"
                    )

                # Add tutorial to section index
                section_text += f"    {tutorial}\n"

                # Add missing images
                tutorial_image_path = os.path.join(tutorial_path, "images")
                if os.path.isdir(tutorial_image_path):
                    doc_image_path = os.path.join(docs_section_path, "images")
                    os.makedirs(doc_image_path, exist_ok=True)

                    for image in os.listdir(tutorial_image_path):
                        image_path = os.path.join(tutorial_image_path, image)
                        dest_path = os.path.join(doc_image_path, image)
                        if not (os.path.exists(dest_path)
                                and filecmp.cmp(image_path, dest_path, shallow=False)):
                            shutil.copyfile(image_path, dest_path)


        # Write section index per tutorials class
        section_index_path = os.path.join(docs_section_path, "index.rst")
        write_if_changed(section_index_path, section_text)


#==============================================================================*

# Remove stale generated content for tutorials/sections that no longer exist
for entry in os.listdir(docs_tutorials_dir):
    entry_path = os.path.join(docs_tutorials_dir, entry)
    if os.path.isdir(entry_path) and os.path.abspath(entry_path) not in expected_sections:
        shutil.rmtree(entry_path)

for section_dir in expected_sections:
    for name in os.listdir(section_dir):
        if name.endswith(".md"):
            page_path = os.path.abspath(os.path.join(section_dir, name))
            if page_path not in expected_pages:
                os.remove(page_path)


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

write_if_changed(main_tutorials_path, main_text)


#==============================================================================*

print(f"python3 {sys.argv[0]} ... End")


#==============================================================================*
