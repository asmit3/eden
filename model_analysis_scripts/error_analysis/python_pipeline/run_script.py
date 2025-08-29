import subprocess
import os

directory_to_search = "/pscratch/sd/k/kkondaka/newfolder/endrapid_test_8_13_v2/tutorial_endrapid/tutorial_data"  # CHANGE THIS


def get_endrapid_paths(base_dir):  # sorts numerically
    return sorted([
        os.path.join(base_dir, name)
        for name in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, name)) and name.startswith("endrapid_")
    ], key=lambda x: int(x.split("_")[-1]))

    
def subprocess_run_py27(script_name):

    try:
        # Equivalent of subprocess.run with output capture
        output = subprocess.check_output(
            ["sbatch", script_name],
            stderr=subprocess.STDOUT
        )
        print("Submitted job:", output.strip())
    except subprocess.CalledProcessError as e:
        print("Job submission failed:")
        print(e.output)


def find_structure_files(directory, extensions):
    file_paths = {ext: [] for ext in extensions}

    for root, _, files in os.walk(directory):
        for file in files:
            for ext in extensions:
                if file.lower().endswith(ext):
                    full_path = os.path.join(root, file)
                    file_paths[ext].append(full_path)

    return file_paths


def run_script_main(directory_to_search):
    end_rapid_folders = get_endrapid_paths(directory_to_search)
    structure_files = find_structure_files(directory_to_search, ('.pdb', '.mtz', '.cif'))

    for i in range(len(end_rapid_folders)):
        script_name = "runlunus{0}.sh".format(i+1)
        script_name = os.path.join(end_rapid_folders[i], script_name)
        print(script_name)

        eff_file = os.path.join(end_rapid_folders[i], "initial_MCR_modified.eff")

        script_content = """#!/bin/bash -l
#SBATCH --job-name=phenix_pipeline
#SBATCH --qos=shared
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --constraint=cpu
#SBATCH -A m3289

source /pscratch/sd/k/kkondaka/phenix-1.20.1/phenix-1.20.1-4487/phenix_env.sh

phenix.refine {0} output.prefix=MCR_XFEL_refine_endrapid_{1}
""".format(eff_file,i+1)
        with open(script_name, "w") as f:
            f.write(script_content)
        subprocess_run_py27(script_name) 

run_script_main(directory_to_search)
