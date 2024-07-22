"""
Script to update OpenFOAM version in file header. Update recursively from
current folder location.

Usage in GeN-Foam root directory:
>>> python3 Tools/caseFormatting/changeHeader.py
"""

#===============================================================================*
# Imports

import os
import sys

#===============================================================================*
# User input

verbose: int = 1

old_version: str = 'v2312'
new_version: str = 'v2406'


#===============================================================================*

script_name: str = sys.argv[0]

folder = os.getcwd()

versions = [
    'v1806', 'v1812', 'v1906', 'v1912', 'v2006', 'v2012', 'v2106', 'v2112',
    'v2206', 'v2212', 'v2306', 'v2312', 'v2406'
]

stats: dict = {'totalFiles': 0}
for v in versions:
    stats[v] = 0


#===============================================================================*

for root, dirs, files in os.walk(folder):
    for filename in files:
        filepath = os.path.join(root, filename)

        # Exclude current script and FMU4FOAM package
        if (filename == script_name or 'FMU4FOAM' in filepath):
            continue

        try:
            with open(filepath, "r") as file:
                file_contents = file.read()

            # Replace specific string
            file_contents = file_contents.replace(
                f"Built on OpenFOAM {old_version}",
                f"Built on OpenFOAM {new_version}"
            )

            # Fill stats
            for version in versions:
                if (version in file_contents):
                    stats[version] += 1
                    break

            with open(filepath, "w") as file:
                file.write(file_contents)
                stats['totalFiles'] += 1


        except UnicodeDecodeError:
            if (verbose >= 2):
                print(f'Skipped file {filepath} due to UnicodeDecodeError')


#===============================================================================*

if (verbose >= 1):
    print('Stats\n-----\n')
    total = sum([stats[version] for version in versions])

    for version in versions:
        print(f"{version:6} {stats[version]:5} | {stats[version] / total * 100:6.2f}")

    print(f"{'':-^21}")
    print(f"Total  {total:5} |")


print("All files updated")

#===============================================================================*
