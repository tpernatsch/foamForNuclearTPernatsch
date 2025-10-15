"""
Script to update OpenFOAM version in file header. Update recursively from
current folder location. Use the WM_PROJECT_VERSION environment variable to
update all headers.

Usage in root directory:
>>> python3 tools/caseFormatting/changeHeader.py
"""

#===============================================================================*
# Imports

import os
import sys
from tabulate import tabulate


#===============================================================================*
# Usefull functions

def nextVersion(version: str) -> str:
    """
    Compute the next version, e.g v2412 -> v2506
    """
    year = int(version[1:3])
    semester = int(version[3:5])

    if (semester == 12):
        return(f"v{year+1}06")
    elif (semester == 6):
        return(f"v{year}12")

    msg = f"Version '{version}' is not a valide version format (vYYSS e.g v2506)"
    raise ValueError(msg)


def previousVersion(version: str) -> str:
    """
    Compute the next version, e.g v2412 -> v2406
    """
    year = int(version[1:3])
    semester = int(version[3:5])

    if (semester == 12):
        return(f"v{year}06")
    elif (semester == 6):
        return(f"v{year-1}12")

    msg = f"Version '{version}' is not a valide version format (vYYSS e.g v2506)"
    raise ValueError(msg)


def getYear(version: str) -> int:
    return("20"+version[1:3])


#===============================================================================*
# User input

verbose: int = 1

new_version: str = os.environ["WM_PROJECT_VERSION"]

new_copyright: str = f'Copyright 2011-2016 OpenFOAM Foundation, 2017-{getYear(new_version)} OpenCFD Ltd.'

strings_to_replace = []

folders_to_exclude = ['thirdParty', 'documentation']


#===============================================================================*

script_name: str = sys.argv[0]

folder = os.getcwd()

# Fill list of versions until the current OpenFOAM version
versions = ['v1806']
while (versions[-1] != new_version):
    strings_to_replace.append(
        {'old': f"Built on OpenFOAM {versions[-1]}", 'new': f"Built on OpenFOAM {new_version}"}
    )
    strings_to_replace.append(
        {'old': f'Copyright 2011-2016 OpenFOAM Foundation, 2017-{getYear(versions[-1])} OpenCFD Ltd.', 'new': new_copyright},
    )

    versions.append(nextVersion(versions[-1]))

stats: dict = {}
for v in versions:
    stats[v] = {'N Files': 0, 'Update to latest': 0}


#===============================================================================*

for root, dirs, files in os.walk(folder):
    for filename in files:
        filepath = os.path.join(root, filename)

        # Exclude current script and thirdParty package
        if (
            filename in script_name
            or any([folder_to_exclude in filepath for folder_to_exclude in folders_to_exclude])
        ):
            continue

        try:
            with open(filepath, "r") as file:
                file_contents = file.read()

            original_content = file_contents

            # Replace specific string
            for string_to_replace in strings_to_replace:
                file_contents = file_contents.replace(string_to_replace['old'], string_to_replace['new'])

            # Fill stats
            if (original_content != file_contents):
                for version in versions[:-1]:
                    if (version in original_content and new_version in file_contents):
                        stats[version]['Update to latest'] += 1
                        break

            for version in versions:
                if (version in file_contents):
                    stats[version]['N Files'] += 1
                    break

            # Break to avoid overwriting non modified files and ask make to
            # recompile files
            if (original_content == file_contents):
                continue

            # Dump file
            with open(filepath, "w") as file:
                file.write(file_contents)


        except UnicodeDecodeError:
            if (verbose >= 2):
                print(f'Skipped file {filepath} due to UnicodeDecodeError')

        except:
            print(f"Impossible to open {filepath}")


#===============================================================================*

if (verbose >= 1):
    data = [{'Version': version} | stat for version, stat in stats.items()]

    nTotalFiles = sum([stats[version]['N Files'] for version in versions])

    data.append({
        'Version': 'Total',
        'N Files': nTotalFiles,
        'Update to latest': sum([stats[version]['Update to latest'] for version in versions])
    })

    for i, line in enumerate(data):
        percent = f"({line['N Files'] / nTotalFiles * 100:.1f} %)"
        data[i]['N Files'] = f"{line['N Files']} {percent:>9}"

    table = tabulate(data, headers="keys", tablefmt="pipe", colalign=['left', 'right', 'right'])

    print(table)


print("All files updated")


#===============================================================================*
