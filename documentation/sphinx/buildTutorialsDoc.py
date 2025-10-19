import os

# Get the absolute path to the directory containing this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the path to the tutorials directory relative to the script
tutorials_dir = os.path.abspath(os.path.join(script_dir, "../../../tutorials")



# Destination base directory for documentation
docs_base_dir = "./usersguide"
docs_tutorials_dir = os.path.join(docs_base_dir, "tutorials")

# Ensure the tutorials documentation directory exists
os.makedirs(docs_tutorials_dir, exist_ok=True)

# Create the main tutorials.rst file
main_tutorials_path = os.path.join(docs_base_dir, "tutorials.rst")
with open(main_tutorials_path, "w") as main_index:

    main_index.write(".. _usersguide_tutorials:\n\n")
    main_index.write("Tutorials\n====================\n\n")
    main_index.write(".. toctree::\n   :maxdepth: 2\n   :caption: Tutorial Sections\n\n")

    # Iterate over section folders in the tutorials source
    for section in sorted(os.listdir(tutorials_dir)):
        section_path = os.path.join(tutorials_dir, section)
        if os.path.isdir(section_path):
            # Create corresponding section folder in documentation
            docs_section_path = os.path.join(docs_tutorials_dir, section)
            os.makedirs(docs_section_path, exist_ok=True)

            main_index.write(f"   tutorials/{section}/index\n")

            # Create section index.rst in documentation
            section_index_path = os.path.join(docs_section_path, "index.rst")
            with open(section_index_path, "w") as section_index:
                section_index.write(f"{section.capitalize()} Tutorials\n{'=' * (len(section) + 10)}\n\n")
                section_index.write(".. toctree::\n   :maxdepth: 1\n\n")

                # Iterate over tutorials in the section
                for tutorial in sorted(os.listdir(section_path)):
                    tutorial_path = os.path.join(section_path, tutorial)
                    
                    if os.path.isdir(tutorial_path):
                        readme_path = os.path.join(tutorial_path, "README.md")
                        tutorial_rst_path = os.path.join(docs_section_path, f"{tutorial}.rst")

                        # Create tutorial .rst file in documentation
                        with open(tutorial_rst_path, "w") as tutorial_rst:
                            tutorial_rst.write(f"{tutorial}\n{'=' * len(tutorial)}\n\n")
                            if os.path.exists(readme_path):
                                with open(readme_path, "r") as readme_file:
                                    tutorial_rst.write(readme_file.read())

                        # Add tutorial to section index
                        section_index.write(f"   {tutorial}\n")

