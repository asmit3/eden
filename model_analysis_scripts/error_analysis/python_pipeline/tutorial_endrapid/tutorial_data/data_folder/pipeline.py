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


MCR_REV_refine_5.mtz
use this to create 100 kicked mtz 
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

directory_to_search = "/Users/yyklab/Desktop/eden/model_analysis_scripts/error_analysis/python_pipeline/tutorial_endrapid/tutorial_data/data_folder"
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

    # extract replacements with safe defaults
    cif_paths = replacements.get('.cif', [])
    pdb_paths = replacements.get('.pdb', [])
    mtz_path = replacements.get('.mtz', [None])[0] or ""
#    fasta_path = replacements.get('.fasta', [None])[0] or ""
    if not pdb_paths:
        print("Warning: No .pdb files found.")
    if not mtz_path:
        print("Warning: No .mtz file found.")
    if not cif_paths:
        print("Warning: No .cif files found.")



    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # step 8a - replace pdb path
        if line.startswith("pdb {"):
            j = 1
            # Delete existing file_name lines
            while i + j < len(lines) and "file_name" in lines[i + j]:
                del lines[i + j]
            # Insert new file_name lines
            for k, pdb_path in enumerate(pdb_paths):
                lines.insert(i + 1 + k, f'      file_name = "{pdb_path}"\n')

        # step 8b - replace mtz path
        elif line.startswith("xray_data {") and i + 1 < len(lines):
            lines[i + 1] = f'      file_name = "{mtz_path}"\n'

        # step 8c - replace CIFs
        elif line.startswith("monomers {"):
            j = 1
        # remove all existing file_name lines
            while i + j < len(lines) and "file_name" in lines[i + j]:
                del lines[i + j]
            # insert new ones
            for k, cif_path in enumerate(cif_paths):
                lines.insert(i + 1 + k, f'      file_name = "{cif_path}"\n')


        # step 8d - delete sequence block
        elif line.startswith("sequence {"):
            j = 1
            while i + j < len(lines) and not lines[i + j].strip().startswith("}"):
                j += 1
            del lines[i:i + j + 1]
            continue

        # step 8f - output block
        elif line.startswith("output {") and i + 2 < len(lines):
            lines[i + 1] = '    prefix = "MCR_XFEL_refine_endrapid"\n'
            lines[i + 2] = '    serial = 0\n'

        # step 8g - refinement parameters
        elif line.startswith("main {") and i + 4 < len(lines):
            lines[i + 1] = '    number_of_macro_cycles = 3\n'
            lines[i + 2] = '    nproc = 1\n'
            lines[i + 3] = '    random_seed = REPLACEME\n'
            lines[i + 4] = '    ordered_solvent = False\n'

        # step 8h - disable target weight optimization
        elif line.startswith("target_weights {") and i + 2 < len(lines):
            lines[i + 1] = '    optimize_xyz_weight = False\n'
            lines[i + 2] = '    optimize_adp_weight = False\n'

        # step 8i - remove gui block
        elif line.startswith("gui {") and i + 4 < len(lines):
            del lines[i:i + 5]
            continue

        i += 1

    # step 9 - append modify_start_model block
    lines.append("\nrefinement.modify_start_model {\n")
    lines.append("  modify {\n")
    lines.append("    sites {\n")
    lines.append('      atom_selection = "not ((resname WAT or resname HOH or resname OOO))"\n')
    lines.append("      shake = 0.25\n")
    lines.append("    }\n")
    lines.append("  }\n")
    lines.append("}\n")

    with open(output_path, 'w') as f:
        f.writelines(lines)



template_file = "initial_MCR.eff"
output_file = "initial_MCR_modified.eff"

update_eff_file(template_file, output_file, structure_files)






