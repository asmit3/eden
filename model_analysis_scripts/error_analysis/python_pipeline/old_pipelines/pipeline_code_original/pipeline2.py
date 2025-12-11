import sys
import os
import numpy as np
import matplotlib.pyplot as plt
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

number_of_macro cycles = 15



"""

#fname = sys.argv[1]

#data_folder = fname 


def create_mtz(mtz_file, output_dir, base_seed=None):
    """
    Generates perturbed MTZ files and returns a mapping of filenames to seeds.

    Args:
        mtz_file (str): Path to the input MTZ file.
        output_dir (str): Directory to save the output files.
        base_seed (int, optional): The starting seed for the random number generator.
                                   If None, a random seed will be generated.

    Returns:
        dict: A dictionary mapping each generated filename to the integer seed used.
              Example: {'kicked1.mtz': 1234, 'kicked2.mtz': 1235, ...}
    """
    if base_seed is None:
        base_seed = np.random.randint(0, 2**32 - 1)
    print("Using base seed: {0}".format(base_seed))

    # --- Load MTZ Data ---
    mtz_obj = mtz.object(''.join(mtz_file))
    dataset = mtz_obj.as_miller_arrays_dict()
    key = ('crystal', 'Experimental-data-used-in-refinement', 'F-obs-filtered')
    fobs_array = dataset[key]
    f_values = np.array(list(fobs_array.data()))
    sigf_values = np.array(list(fobs_array.sigmas()))
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Initialize an empty dictionary to store the seed mapping
    seed_map = {}

    num_simulations = 100
    for i in range(num_simulations):
        # Calculate a unique and deterministic seed for this file
        current_seed = base_seed + i
        np.random.seed(current_seed)

        # Generate a single perturbed dataset
        fobs_simulation = np.random.normal(loc=f_values, scale=sigf_values)

        # --- Create and write the new MTZ file ---
        x_flex = flex.double(fobs_simulation)
        new_miller_array = fobs_array.customized_copy(data=x_flex)
        mtz_dataset = new_miller_array.as_mtz_dataset(column_root_label="F-obs-filtered")

        filename = "kicked{0}.mtz".format(i+1)
        output_path = os.path.join(output_dir, filename)
        mtz_dataset.mtz_object().write(output_path)
        
        # 2. Add the filename and its seed to our dictionary
        seed_map[filename] = current_seed

    print("\n Successfully generated {0} MTZ files in '{1}'.".format(num_simulations, output_dir))
    print("   Returning the seed map dictionary.")

    # 3. Return the completed dictionary to the caller
    return seed_map


def update_eff_file(original_mtz, kicked_mtz, template_path, output_path, replacements, seed):
    with open(template_path, 'r') as f:
        lines = f.readlines()

    # extract replacements with safe defaults
    cif_paths = replacements.get('.cif', [])
    pdb_paths = replacements.get('.pdb', [])
    mtz_files = replacements.get('.mtz', [])
    mtz_path = original_mtz['.mtz'][0]
    
    kicked_path = kicked_mtz
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
            # delete existing file_name lines
            while i + j < len(lines) and "file_name" in lines[i + j]:
                del lines[i + j]
            # insert new file_name lines
            for k, pdb_path in enumerate(pdb_paths):
                lines.insert(i + 1 + k, '      file_name = "{0}"\n'.format(pdb_path))

        # step 8b - replace mtz path
        elif line.startswith("xray_data {") and i + 1 < len(lines):
            lines[i + 1] = '      file_name = "{0}"\n'.format(kicked_path)
            lines[i + 2] = '      labels = "F-obs-filtered,SIGF-obs-filtered"\n'
            

        # replacxe R-free flags with the mtz file
        elif line.startswith("r_free_flags {"): 
            print("here")
            lines[i + 1] = '      file_name = "{0}"\n'.format(mtz_path)


        # step 8c - replace CIFs
        elif line.startswith("monomers {"):
            j = 1
        # remove all existing file_name lines
            while i + j < len(lines) and "file_name" in lines[i + j]:
                del lines[i + j]
            # insert new ones
            for k, cif_path in enumerate(cif_paths):
                lines.insert(i + 1 + k, '      file_name = "{0}"\n'.format(cif_path))


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
            lines.insert(i + 1, '    number_of_macro_cycles = 15\n')
            lines.insert(i + 2, '    nproc = 1\n')
            lines.insert(i + 3, '    random_seed = {0}\n'.format(seed))
            lines.insert(i + 4,'    ordered_solvent = False\n')

        # step 8h - disable target weight optimization
        elif line.startswith("target_weights {") and i + 2 < len(lines):
            lines[i + 1] = '    optimize_xyz_weight = False\n'
            lines[i + 2] = '    optimize_adp_weight = False\n'

        # step 8i - remove gui block
        elif line.startswith("gui {") and i + 4 < len(lines):
            del lines[i:i + 8]
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
        folder_name = "endrapid_{0}".format(i)
        folder_path = os.path.join(base_dir, folder_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_path)

def get_kicked_number(path):
    filename = os.path.basename(path)
    if filename.startswith("kicked") and filename.endswith(".mtz"):
        number_part = filename[len("kicked"):-len(".mtz")]
        if number_part.isdigit():
            return int(number_part)
    return float('inf')  # Push non-matching files to the end


import shutil
def main():
    directory_to_search = "/pscratch/sd/k/kkondaka/newfolder/endrapid_test_8_28/tutorial_endrapid/tutorial_data"

    # 1. Ensure all 100 endrapid_* folders exist
    create_endrapid(directory_to_search)

    # 2. Clean up any old kicked*.mtz in the root before regenerating
    for f in os.listdir(directory_to_search):
        if f.startswith("kicked") and f.endswith(".mtz"):
            os.remove(os.path.join(directory_to_search, f))

    # 3. Get original mtz file (non-kicked) and generate fresh kicked*.mtz
    original_mtz = find_structure_files(directory_to_search, ['.mtz'])
    seed_map = create_mtz(original_mtz['.mtz'][0], directory_to_search)

    # 4. Get sorted list of kicked*.mtz from root directory
    kicked_files = [
        os.path.join(directory_to_search, f)
        for f in os.listdir(directory_to_search)
        if f.startswith("kicked") and f.endswith(".mtz")
    ]
    kicked_files.sort(key=get_kicked_number)

    # 5. Get all other structure files (.pdb, .mtz, .cif)
    structure_files = find_structure_files(directory_to_search, ('.pdb', '.mtz', '.cif'))

    # Remove any kicked*.mtz from structure_files['.mtz'] list
    structure_files['.mtz'] = [
        f for f in structure_files['.mtz']
        if not os.path.basename(f).startswith("kicked")
    ]

    # 6. Get sorted list of endrapid_* folders
    end_rapid_folders = get_endrapid_paths(directory_to_search)

    # 7. Template file path
    template_file = os.path.join(directory_to_search, "initial_MCR.eff")

    # 8. Move each kicked file into its corresponding folder & update .eff
    for i, end_rapid_path in enumerate(end_rapid_folders, start=1):
        kicked_src = kicked_files[i - 1]
        kicked_dst = os.path.join(end_rapid_path, os.path.basename(kicked_src))
        shutil.move(kicked_src, kicked_dst)

        output_file = os.path.join(end_rapid_path, "initial_MCR_modified.eff")
        seed = seed_map["kicked{0}.mtz".format(i)]

        print("Moving:", kicked_src, "->", kicked_dst)
        print("Seed:", seed)

        update_eff_file(
            original_mtz, kicked_dst,
            template_file, output_file,
            structure_files, seed
        )


main()

