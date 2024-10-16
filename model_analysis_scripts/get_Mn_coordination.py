from iotbx.data_manager import DataManager
from scitbx.matrix import col
import matplotlib.pyplot as plt
import copy
import sys
from scitbx.array_family import flex
import numpy as np

selection_to_use = "((chain a or chain A) and resid 601 and (name MN1 or name MN2 or name MN3 or name MN4))"
fn = sys.argv[1]
dm = DataManager()
dm.set_overwrite(True)
model = dm.get_model(fn)
ph = model.get_hierarchy()
sele_ph = ph.apply_atom_selection(selection_to_use)
sele_atoms = sele_ph.atoms()
sele_xyz = sele_atoms.extract_xyz()
sele_names = sele_atoms.extract_name()
all_xyz = ph.atoms().extract_xyz()
all_names = ph.atoms().extract_name()
all_elements = ph.atoms().extract_element()

cutoff_distance = 2.5 # Cutoff used for O = 2.5, C=3.2, Mn=2.9, 3.5, Ca=3.6
element_of_interest = 'O'
coordination_distances_element_specific = flex.double()
for i in range(sele_xyz.all()[0]):
  dxyz = (all_xyz - sele_xyz[i]).norms()
  Mn_coordination_selection = (dxyz<cutoff_distance)
  coordination_distances = dxyz.select(Mn_coordination_selection)
  distances_xyz = dxyz.select(Mn_coordination_selection)
  coordination_names = all_names.select(Mn_coordination_selection)
  coordination_elements = all_elements.select(Mn_coordination_selection)
  for j, element in enumerate(coordination_elements):
    # Uncomment below for Mn-O distances
    #if element.strip() == element_of_interest and coordination_names[j].strip() in ['O1', 'O2', 'O3', 'O4', 'O5']:

    # Uncomment below for Mn-O(N) distances
    if element.strip() == element_of_interest and coordination_names[j].strip() not in ['O1', 'O2', 'O3', 'O4', 'O5']:
    # Uncomment below for Mn-Mn, Mn-C, Mn-CA distances
    #if element.strip() == element_of_interest:
      coordination_distances_element_specific.append(coordination_distances[j])
      print (coordination_names[j], sele_names[i])

coordination_distances_element_specific=coordination_distances_element_specific.select(coordination_distances_element_specific > 0.2)

coordination_distances_element_specific = list(set(coordination_distances_element_specific))

mean_dist = np.mean(coordination_distances_element_specific)
mean_std = np.std(coordination_distances_element_specific)
print ('Mean Stdev of Mn-%s is %.2f +/- %.2f (N_samples = %d)'%(element_of_interest, mean_dist, mean_std, len(coordination_distances_element_specific)))



plt.hist(coordination_distances_element_specific)
plt.xlabel('Coordination distances for Mn-%s'%element_of_interest)
plt.show()

