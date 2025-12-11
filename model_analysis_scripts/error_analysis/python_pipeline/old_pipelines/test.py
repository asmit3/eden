import sys
import os
import numpy as np
import matplotlib as plt
from iotbx import mtz
from scitbx.array_family import flex
import subprocess



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

seed = np.random.randint(0, 2**32 - 1)  

def create_mtz(mtz_file, output_dir):
    print("Using seed: ", seed)
       
    
    np.random.seed(seed)

    mtz_obj = mtz.object(mtz_file)


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
        
        output_path = os.path.join(output_dir, f"kicked{i+1}.mtz")
        mtz_dataset.mtz_object().write(output_path)


create_mtz("/pscratch/sd/k/kkondaka/newfolder/python_test_7_25/0F.mtz", "/pscratch/sd/k/kkondaka/newfolder/python_test_7_25")