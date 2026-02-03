import os
import subprocess
import fileinput
import shutil
import re
import sys
import numpy as np

# Function to modify dictionary using foamDictionary
def modify_foam_dict(dict_rel_path, entry, value):
    subprocess.run(["foamDictionary", dict_rel_path, "-entry", entry, "-set", str(value)], check=True)            

# Function to run simulation
def run_simulation():
    # Make the .Allrun script executable and run it
    subprocess.run(["chmod", "+x", "./Allrun"], check=True)
    subprocess.run(["./Allrun"], check=True)

# Function to list and remove time directories
def clean_time_directories():
    # Run foamListTimes with the -rm option to remove time directories
    subprocess.run(["foamListTimes", "-rm"], check=True)    

def extract_internal_field(case_path, field_file):
    """
    Use foamDictionary to extract the internalField from an OpenFOAM field file.
    :param case_path: Path to the case directory.
    :param field_file: Name of the field file to extract data from.
    :return: List of floats representing the internalField data.
    """
    command = ["foamDictionary", case_path + field_file, "-entry", "internalField", "-value", "-writePrecision", "20"]
    result = subprocess.run(command, capture_output=True, text=True)
    data_str = result.stdout

    # Determine if the field is scalar or vector by checking for "List<vector>" in the data string
    if 'List<vector>' in data_str:
        # Field is vectorial
        # Remove the "nonuniform List<vector>" part and any other non-numeric, non-space characters except for newlines, '(', and ')'
        data_str = re.sub(r'nonuniform List<vector>\s*\d*\s*\(', '', data_str)
        data_str = re.sub(r'[^\d.\-e \n()]+', '', data_str)  # Keep numbers, '.', '-', 'e', ' ', new line, '(', and ')'

        # Convert the cleaned string to a list of tuples, each representing a vector
        data_list = [tuple(map(float, item.strip("()").split())) for item in data_str.split('\n') if item.strip()]
    else:
        # Field is scalar
        data_str = re.sub(r'nonuniform List<scalar>\s*\d*\s*\(', '', data_str)
        data_str = re.sub(r'[^\d.\-e\n]+', '', data_str)  # Keep numbers, '.', '-', 'e', and new line

        # Convert the cleaned string to a list of floats
        data_list = [float(item) for item in data_str.split()]

    return data_list

def compare_internal_fields(case_path, file1, file2, tolerance=1e-6):
    """
    Compare the internalField entries of two field files within a specified tolerance.
    :param case_path: Path to the case directory.
    :param file1: First field file name.
    :param file2: Second field file name.
    :param tolerance: Tolerance for numerical comparison.
    :return: Boolean indicating if the internalField entries are identical within the tolerance.
    """
    field1 = extract_internal_field(case_path, file1)
    field2 = extract_internal_field(case_path, file2)

    # The two fields must have the same length to be considered identical
    if len(field1) != len(field2):
        return False

    # Determine if we are dealing with scalar or vector fields based on the first element
    if isinstance(field1[0], tuple):
        # This is a vector field comparison
        for vec1, vec2 in zip(field1, field2):
            # Convert tuples to numpy arrays for vectorized operations
            vec1, vec2 = np.array(vec1), np.array(vec2)
            # Calculate the vector difference and then the magnitude of the vector
            if np.linalg.norm(vec1 - vec2)/max(np.linalg.norm(vec1), 1e-15) > tolerance:
                print(vec1)
                print(vec2)
                return False
    else:
        # This is a scalar field comparison
        for val1, val2 in zip(field1, field2):
            if abs(val1 - val2)/max(abs(val1), 1e-15) > tolerance:
                print(val1)
                print(val2)
                return False

    return True

def clean_directory_except_base(base_case_folder):
    """
    Remove all directories and specific file ('results') in the current working directory except for the base_case_folder.
    :param base_case_folder: Folder to exclude from deletion.
    """
    base_path = os.getcwd()
    results_file = 'results'  # Filename you want to delete

    # Ensure base_path exists and is a directory
    if os.path.exists(base_path) and os.path.isdir(base_path):
        for item in os.listdir(base_path):
            item_path = os.path.join(base_path, item)
            # Check if the current item is the 'results' file
            if item == results_file and os.path.isfile(item_path):
                os.remove(item_path)
                print(f"Removed {item_path}")
            # Ensure the current item is a directory and not the base_case_folder
            elif os.path.isdir(item_path) and item != base_case_folder:
                # Remove the folder
                shutil.rmtree(item_path)
                print(f"Removed {item_path}")
    else:
        print(f"Path {base_path} does not exist or is not a directory.")

def main():    

    # Check if any command-line argument is provided and if it's "clean"
    if len(sys.argv) > 1 and sys.argv[1] == "clean":
        print("Cleaning mode activated. Only cleaning directories...")
        base_case_dir = "baseCase"  # specify your base case directory
        clean_directory_except_base(base_case_dir)
        print("Cleaning done.")
        return  # Exit the script after cleaning, don't proceed further
    
    # Continue the script if no "clean" options is given
    base_case_dir = "baseCase"
    work_dir = "."
    base_path = os.getcwd()

    # Clean the directory, leaving only the base case
    clean_directory_except_base(base_case_dir)    

    # List of options to test
    options_to_test = [
        {"dict": ["constant/solverDict"], "entry": ["thermalSolver"], "value": ["solidConduction"]},
        {"dict": ["constant/solverDict"], "entry": ["gapGas"], "value": ["Frapcon"]},
        {"dict": ["constant/solverDict"], "entry": ["heatSource"], "value": ["timeDependentLhgr"]},
        {"dict": ["constant/solverDict"], "entry": ["fastFlux"], "value": ["timeDependentAxialProfile"]},
        {"dict": ["constant/solverDict"], "entry": ["heatSourceOptions/axialProfile/type"], "value": ["timeDependentTabulated"]}, 
        {"dict": ["0/T", "0/T"], "entry": ["boundaryField/fuelOuter/type", "boundaryField/cladInner/type"], "value": ["fuelRodGap", "fuelRodGap"]}, 
        {"dict": ["constant/solverDict"], "entry": ["burnup"], "value": ["fromPower"]},
        {"dict": ["constant/solverDict","constant/solverDict",] , "entry": ["neutronicsSolver","burnup"], "value": ["diffusion","Lassmann"]},    
        {"dict": ["constant/solverDict"], "entry": ["heatSourceOptions/radialProfile/type"], "value": ["fromBurnup"]},
        {"dict": ["constant/solverDict"], "entry": ["mechanicsSolver"], "value": ["smallStrain"]},
        {"dict": ["0/D", "0/D"], "entry": ["boundaryField/fuelOuter/type", "boundaryField/cladInner/type"], "value": ["gapContact", "gapContact"]}, 
        {"dict": ["constant/solverDict"], "entry": ["materials/fuel/isotropicCracking"], "value": ["on"]}, 
        {"dict": ["constant/solverDict"], "entry": ["fgr"], "value": ["SCIANTIX"]},
        {"dict": ["constant/solverDict"], "entry": ["fgrOptions/addToRegistry"], "value": ["on"]}, 
        {"dict": ["constant/solverDict"], "entry": ["materials/cladding/rheologyModelOptions/creepModel"], "value": ["misesPlasticCreep"]}, 
        # ... other options
    ]

    caseN = 0
    for option in options_to_test:

        dict_list = option["dict"]
        entry_list = option["entry"]
        value_list = option["value"]
        
        caseN +=1

        # Copy base case to a new working directory
        new_case_dir = os.path.join(work_dir, f"case_{caseN}_{entry_list[0].replace('/', '_')}_{value_list[0]}")
        shutil.copytree(base_case_dir, new_case_dir)

        print("\nEntering in " + str(new_case_dir))
        os.chdir(new_case_dir)

        print("\nChanging dictionaries:")
        # Modify the dictionary using foamDictionary
        # Check if all lists are of the same length to ensure one-to-one correspondence
        if len(dict_list) == len(entry_list) == len(value_list):
            for dict, entry, value in zip(dict_list, entry_list, value_list):
                modify_foam_dict(dict, entry, value)
        else:
            print("Error: Mismatched options length.")

        print("\n1-step simulation: running the Allrun script")
        # Run the .Allrun script (assumes it handles the simulation from 0 to 1 year)
        clean_time_directories()
        run_simulation()

        print("\n1-step simulation: saving results")
        # Save results
        shutil.copy("1576800/T", "results/T_1step")
        shutil.copy("1576800/D", "results/D_1step")
        shutil.copy("1576800/Bu", "results/Bu_1step")

        # Run the half-year simulation
        # You might need to modify something to tell .Allrun to stop at half a year
        print("\n2-steps simulation: running the Allrun script for 1st part")
        clean_time_directories()
        modify_foam_dict("system/controlDict", "endTime", "788400")
        run_simulation()

        # You might need to modify something to tell .Allrun to restart from latest time
        print("\n2-steps simulation: running the 2nd part")
        modify_foam_dict("system/controlDict", "endTime", "1576800")
        run_simulation()

        print("\n2-steps simulation: saving results")
        # Save results
        shutil.copy("1576800/T", "results/T_2step")
        shutil.copy("1576800/D", "results/D_2step")
        shutil.copy("1576800/Bu", "results/Bu_2step")

        # Perform the comparison
        fields_to_compare = [
            ("T_1step", "T_2step"),
            ("Bu_1step", "Bu_2step"),
            ("D_1step", "D_2step"),
            # add more pairs if needed
        ]

        results_file_path = os.path.join(base_path, "results")

        # It's good practice to use 'with open' as it handles file closure even if exceptions are raised
        with open(results_file_path, "a") as results_file:

            header_message = f"\nResults for "+str(new_case_dir)+":\n"
            print(header_message)
            results_file.write(header_message)

            for file1, file2 in fields_to_compare:
                if compare_internal_fields("results/", file1, file2):
                    message = f"{file1} and {file2} have identical internalField entries.\n"
                    print(message)
                    results_file.write(message)
                else:
                    message = f"{file1} and {file2} have different internalField entries.\n"
                    print(message)
                    results_file.write(message)    

        os.chdir(base_path)
        base_case_dir = new_case_dir


if __name__ == "__main__":
    main()
