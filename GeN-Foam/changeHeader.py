import os

folder = "/home/faculty/c/carlo.fiorina/data/GeN-Foam/"

for root, dirs, files in os.walk(folder):
    for filename in files:
        filepath = os.path.join(root, filename)

        try:
            with open(filepath, "r") as file:
                file_contents = file.read()

            file_contents = file_contents.replace("Built on OpenFOAM v2306", "Built on OpenFOAM v2306")

            with open(filepath, "w") as file:
                file.write(file_contents)

        except UnicodeDecodeError:
            print(f'Skipped file {filepath} due to UnicodeDecodeError')
