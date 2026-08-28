from __future__ import absolute_import, division, print_function
from iotbx.data_manager import DataManager
from iotbx.pdb import hierarchy
import sys, os
from scitbx.array_family import flex
from scitbx.math import dihedral_angle
import glob

def filter_resname(resname):
  if resname[0] in ['C', 'N']: return 'IGN'
  if resname in ['HIS', 'HIP', 'HID', 'HIE']: return 'HIS'
  return resname

#http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
amino_acid_chi_dihedrals = {
    "ALA": {},
    "ARG": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD"],
        "chi3": ["CB", "CG", "CD", "NE"],
        "chi4": ["CG", "CD", "NE", "CZ"],
    },
    "ASN": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "OD1"],
    },
    "ASP": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "OD1"],
    },
    "CYS": {
        "chi1": ["N", "CA", "CB", "SG"],
    },
    "GLN": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD"],
        "chi3": ["CB", "CG", "CD", "OE1"],
    },
    "GLU": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD"],
        "chi3": ["CB", "CG", "CD", "OE1"],
    },
    "GLY": {},
    "HIS": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "ND1"],
    },
    "ILE": {
        "chi1": ["N", "CA", "CB", "CG1"],
        "chi2": ["CA", "CB", "CG1", "CD"],
    },
    "LEU": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD1"],
    },
    "LYS": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD"],
        "chi3": ["CB", "CG", "CD", "CE"],
        "chi4": ["CG", "CD", "CE", "NZ"],
    },
    "MET": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "SD"],
        "chi3": ["CB", "CG", "SD", "CE"],
    },
    "PHE": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD1"],
    },
    "PRO": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD"],
    },
    "SER": {
        "chi1": ["N", "CA", "CB", "OG"],
    },
    "THR": {
        "chi1": ["N", "CA", "CB", "OG1"],
    },
    "TRP": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD1"],
    },
    "TYR": {
        "chi1": ["N", "CA", "CB", "CG"],
        "chi2": ["CA", "CB", "CG", "CD1"],
    },
    "VAL": {
        "chi1": ["N", "CA", "CB", "CG1"],
    },
}

def dihedral_info_from_single_pdb(pdb_file_name):
  dm = DataManager()
  dm.set_overwrite(True)          
  model_filename = os.path.join(pdb_file_name)    
  model = dm.get_model(model_filename) 
  m = model.deep_copy() 
  pdb_hierarchy = m.get_hierarchy()    
  chi_angle_dict = {}  

  for chain in pdb_hierarchy.only_model().chains():
    for residue_group in chain.residue_groups():
      resname = residue_group.unique_resnames()[0].strip()
      resname = filter_resname(resname)
      if resname == 'IGN': continue #IGNORING C and N terminus
      assert len(list(residue_group.unique_resnames())) == 1, 'Multiple resnames for some reason'
      #print (resname)
      assert resname in amino_acid_chi_dihedrals.keys(), 'Resname is not in side chain dihedral definition'
      chi_angle_dict_key = '%s_%d'%(resname, residue_group.resseq_as_int())
      chi_angle_dict[chi_angle_dict_key] = [] 
      for atom_group in residue_group.atom_groups():
        atom_dict = {atom.name.strip(): atom for atom in atom_group.atoms()}
        names = atom_group.atoms().extract_name()
        names = [name.strip() for name in names]
        coords = atom_group.atoms().extract_xyz()
        atom_xyz_dict = {name:xyz for name,xyz in zip(names,coords)}
        for chi_angle in amino_acid_chi_dihedrals[resname].keys():
          sites = flex.vec3_double()
          chi_atoms = amino_acid_chi_dihedrals[resname][chi_angle]
          # First check that all 4 atoms exist in the residue
          for chi_atom in chi_atoms:
            assert chi_atom in names, 'Chi atom does not exist in PDB file'
            sites.append(atom_xyz_dict[chi_atom])
          #print (atom_xyz_dict)
          angle = dihedral_angle(sites=sites, deg=True)
          chi_angle_dict[chi_angle_dict_key].append(angle)
  return chi_angle_dict

def dihedral_info_all():
  pdb_file_names = glob.glob(os.path.join(base_path, '*.pdb'))
  all_chi_angle_dict = {}
  for i, pdb_file_name in enumerate(pdb_file_names):
    chi_angle_dict = dihedral_info_from_single_pdb(pdb_file_name)
    if i == 0:
      for uniq_res in chi_angle_dict.keys():
        all_chi_angle_dict[uniq_res] = []
        all_chi_angle_dict[uniq_res].append(chi_angle_dict[uniq_res])
    else:
      for uniq_res in all_chi_angle_dict.keys():
        all_chi_angle_dict[uniq_res].append(chi_angle_dict[uniq_res])
  return all_chi_angle_dict
      
def print_stats(all_chi_angle_dict):
  for uniq_res in all_chi_angle_dict.keys():
    if uniq_res == 'GLN_329':
      data = all_chi_angle_dict[uniq_res]
      num_chi_angles = len(data[0])
      
      chi_vals = {i:[] for i in range(num_chi_angles)}
      for i, val in enumerate(data):
        for n_chi in range(num_chi_angles):
          chi_vals[n_chi].append(val[n_chi])
      # Plot Dihedral angles
      import matplotlib.pyplot as plt
      
      fig, ax = plt.subplots(num_chi_angles, 1)
      ax = ax.ravel()
      for ii, a in enumerate(ax):
        a.hist(chi_vals[ii])
        a.set_xlabel('Angles')
        a.set_title('%s: Chi_%d'%(uniq_res, ii+1))
      plt.show()
       
def run():
  all_chi_angle_dict = dihedral_info_all()
  print_stats(all_chi_angle_dict)

if __name__ == '__main__':
  base_path='/Users/AB/'
  pdb_file_name='frame_000050.pdb'
  run()
