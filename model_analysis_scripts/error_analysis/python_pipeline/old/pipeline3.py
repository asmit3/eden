import sys
import os
import numpy as np
import matplotlib as plt
from iotbx import mtz
from scitbx.array_family import flex



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


MCR_REV_refine_5.mtz - use this to extract rfree flags
use this to create 100 kicked mtz 
"""

#fname = sys.argv[1]

#data_folder = fname 
seed = np.random.randint(0, 2**32 - 1)  

def create_mtz(mtz_file, output_dir):
    print("Using seed: ", seed)
    np.random.seed(seed)

    """
    UNCOMMENT LATER WITH INPUT AS WELL!
    try: 
        user_input = str(input("Enter a seed (or press Enter to use random seed): ")).strip()
    except SyntaxError as e: 
        print("No input detected. Using random seed:", seed)
    """
    
    # Initialize RNG with that seed
    np.random.seed(seed)
    print(mtz_file)
    mtz_obj = mtz.object(''.join(mtz_file))
    # /Users/kanishkkondaka/Desktop/phenix/0F.mtz 


    dataset = mtz_obj.as_miller_arrays_dict()

    key = ('crystal', 'Experimental-data-used-in-refinement', 'F-obs-filtered') # Fmodel 
    #miller_array = dataset[key]
    fobs_array = dataset[key]




    data_values = fobs_array.data()

    f_values = np.array(list(data_values))


    sigmas = fobs_array.sigmas()
    sigf_values = np.array(list(sigmas))
    

    fobs_simulations = np.random.normal(loc=f_values[None, :], scale=sigf_values[None, :], size=(100, len(f_values)))

    
    os.makedirs(output_dir, exist_ok=True)

    for i in range(len(fobs_simulations)): 
        x_flex = flex.double(fobs_simulations[i,:])
        new_miller_array = fobs_array.customized_copy(data=x_flex)

        mtz_dataset = new_miller_array.as_mtz_dataset(column_root_label="Fobs_perturbed")
        
        # Write file to 'endrapid_1/kicked1.mtz', etc.
        output_path = os.path.join(output_dir, f"kicked{i+1}.mtz")
        mtz_dataset.mtz_object().write(output_path)
    


def update_eff_file(original_mtz, template_path, output_path, replacements):
    with open(template_path, 'r') as f:
        lines = f.readlines()

    # extract replacements with safe defaults
    cif_paths = replacements.get('.cif', [])
    pdb_paths = replacements.get('.pdb', [])
    mtz_files = replacements.get('.mtz', [])
    mtz_path = original_mtz
    kicked_path = mtz_files[0] if len(mtz_files) > 0 else ""
#   fasta_path = replacements.get('.fasta', [None])[0] or ""
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
            lines[i + 1] = f'      file_name = "{kicked_path}"\n'
            lines[i + 2] = f'      labels = "F-obs-filtered,SIGF-obs-filtered"\n'

        # replacxe R-free flags with the mtz file
        elif line.startswith("r_free_flags {"): 
            print("here")
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



def find_structure_files(directory, extensions):
    file_paths = {ext: [] for ext in extensions}

    for root, _, files in os.walk(directory):
        for file in files:
            for ext in extensions:
                if file.lower().endswith(ext):
                    full_path = os.path.join(root, file)
                    file_paths[ext].append(full_path)

    return file_paths

def get_endrapid_paths(base_dir): # sorts numerically, previous version sorts lexiconally
    return sorted([
        os.path.join(base_dir, name)
        for name in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, name)) and name.startswith("endrapid_")
    ], key=lambda x: int(x.split("_")[-1]))


import os
def create_endrapid(base_dir): 
#    base_dir = os.getcwd()  # Or set to a specific path

    # Loop through 1 to 100 and create folders + copy files
    for i in range(1, 101):
        folder_name = f"endrapid_{i}"
        folder_path = os.path.join(base_dir, folder_name)
        os.makedirs(folder_path, exist_ok=True)


def main():
    directory_to_search = "/Users/yyklab/Desktop/eden/model_analysis_scripts/error_analysis/python_pipeline/tutorial_endrapid/tutorial_data" # CHANGE THIS!
    #/pscratch/sd/k/kkondaka/newfolder/python_endrapid_test_7_26/kanishk_endrapid/tutorial_endrapid/tutorial_data
    
    create_endrapid(directory_to_search)
    end_rapid_folders = get_endrapid_paths(directory_to_search)



    original_mtz = find_structure_files(directory_to_search, ('.mtz',))
    print(original_mtz)
    create_mtz(original_mtz['.mtz'], directory_to_search)

    structure_files = find_structure_files(directory_to_search, ('.pdb', '.mtz', '.cif'))

    structure_files['.mtz'] = [f for f in structure_files['.mtz'] if not f.endswith("kicked.mtz")]
    """
        print(structure_files)
        for ext, paths in structure_files.items():
            print(f"{ext.upper()} files:")
            for path in paths:
                print(f"  {path}")
     """
    for j in range(len(end_rapid_folders)):
        end_rapid_path = end_rapid_folders[j]
        template_file = "/Users/yyklab/Desktop/eden/model_analysis_scripts/error_analysis/python_pipeline/tutorial_endrapid/tutorial_data/initial_MCR.eff"
        #/pscratch/sd/k/kkondaka/newfolder/python_endrapid_test_7_26/kanishk_endrapid/tutorial_endrapid/tutorial_data
        output_file = os.path.join(end_rapid_path, "initial_MCR_modified.eff")
        update_eff_file(''.join(original_mtz['.mtz']), template_file, output_file, structure_files)

main()
"""
1. use these folders endrapid_1 -> kicked.mtz, write out a modified phil file
2. there a bunch of endrapid folders with like 100
3. phenix.refine params.phil, run batch jobs on the command line
"""



