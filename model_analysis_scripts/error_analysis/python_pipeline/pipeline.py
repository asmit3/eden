import sys
import os


# specify the pdb file has to be the full path pdb = files
# also mtz file for that , mtz = mtz path
# data-folder = some path, probably easier
# this folder has pdb, mtz, and cif files
# in that folder you search for the pdb file, mtz file, and multiple cif files
# finds these files and assigns them to the variables in the phil files.
"""
stores every line into an array, and write out a new phil file called ckanishk.phil
step 8a
input {
    pdb {
1X    file_name = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/Refine_4/MCR_REV_refine_4-coot-0.pdb"
    }

fchange the file_name here 

ignore the R free flags

search for file_name and look for a .pdb a extenion
find something that ends with a .pdb
the kicked.mtz will eventually be one of the 100 ouptut files

mainly just changing the paths
use protocol 2

step 8b

was told to delete sequence

"""

#fname = sys.argv[1]

#data_folder = fname 


def find_structure_files(directory, extensions=('.pdb', '.mtz', '.cif','.fasta')):
    file_paths = {ext: [] for ext in extensions}

    for root, _, files in os.walk(directory):
        for file in files:
            for ext in extensions:
                if file.lower().endswith(ext):
                    full_path = os.path.join(root, file)
                    file_paths[ext].append(full_path)

    return file_paths

directory_to_search = "/Users/yyklab/Desktop/Lab_Files/python_pipeline/data_folder"
structure_files = find_structure_files(directory_to_search)
"""
strucutre_files is an array with the absolute paths to each mtz, pdb, and cif file. 
example: 
{'.pdb': ['/Users/yyklab/Desktop/Lab_Files/python_pipeline/data_folder/example.pdb'], 
'.mtz': ['/Users/yyklab/Desktop/Lab_Files/python_pipeline/data_folder/example.mtz'], 
'.cif': ['/Users/yyklab/Desktop/Lab_Files/python_pipeline/data_folder/example2.cif', 
'/Users/yyklab/Desktop/Lab_Files/python_pipeline/data_folder/example1.cif']}
"""

# print results
for ext, paths in structure_files.items():
    print(f"{ext.upper()} files:")
    for path in paths:
        print(f"  {path}")


def update_eff_file(template_path, output_path, replacements):
    with open(template_path, 'r') as f:
        lines = f.readlines()
    print(lines)

    # extract replacements with safe defaults
    cif_paths = replacements.get('.cif', [])
    pdb_path = replacements.get('.pdb', [None])[0] or ""
    mtz_path = replacements.get('.mtz', [None])[0] or "" # code in multiple mtz
    fasta_path = replacements.get('.fasta', "")

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        # step 8a
        if line == "pdb {" and i + 1 < len(lines):
            lines[i + 1] = f'      file_name = "{pdb_path}"\n'

        # step 8b
        elif line == "xray_data {" and i + 1 < len(lines):
            lines[i + 1] = f'      file_name = "{mtz_path}"\n'

        # step 8c - replace multiple CIFs
        elif line == "monomers {":
            j = 1
            for k, cif_path in enumerate(cif_paths):
                idx = i + j
                if idx < len(lines) and "file_name" in lines[idx]:
                    lines[idx] = f'      file_name = "{cif_path}"\n'
                    j += 1
                else:
                    break
        
        # step 8d - fasta
        """
        DELETING SEQUENCE INSTEAD
        # step 8d - fasta
        elif line == "sequence {" and i + 1 < len(lines): # delete the sequence
            lines[i + 1] = f'      file_name = "{fasta_path}"\n'
        """
        elif line == "sequence {" and i + 1 < len(lines):
            del lines[i:i+2]

        # step 8f - output block
        elif line == "output {" and i + 2 < len(lines):
            lines[i + 1] = '    prefix = "MCR_XFEL_refine_endrapid"\n'
            lines[i + 2] = '    serial = 0\n'

        # step 8g - main refinement parameters
        elif line == "main {" and i + 4 < len(lines):
            lines[i + 1] = '    number_of_macro_cycles = 3\n'
            lines[i + 2] = '    nproc = 1\n'
            lines[i + 3] = '    random_seed = REPLACEME\n'
            lines[i + 4] = '    ordered_solvent = False\n'

        # step 8h - disable optimization weights
        elif line == "target_weights {" and i + 2 < len(lines):
            #115
            lines[i + 1] = '    optimize_xyz_weight = False\n'
            lines[i + 2] = '    optimize_adp_weight = False\n'

        # step 8i - remove GUI block
        elif line == "gui {" and i + 4 < len(lines):
            del lines[i:i + 5]
            continue  # lines shift when deleted!

        i += 1

    # step 9 - add refinement function
    lines.append("\nrefinement.modify_start_model {\n")
    lines.append("  modify {\n")
    lines.append("    sites {\n")
    lines.append('      atom_selection = "not ((resname WAT or resname HOH or resname OOO))"\n')
    lines.append("      shake = 0.25\n")
    lines.append("    }\n")
    lines.append("  }\n")
    lines.append("}\n")

    # write updated file
    with open(output_path, 'w') as f:
        f.writelines(lines)



                

    with open(output_path, 'w') as f:
        f.writelines(lines)

directory = "/path/to/your/folder"
template_file = "initial_MCR.eff"
output_file = "initial_MCR_modified.eff"

paths = find_structure_files(directory)
update_eff_file(template_file, output_file, structure_files)






