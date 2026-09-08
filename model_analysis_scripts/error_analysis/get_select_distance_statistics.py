# Go through all the PDB files from endrapid and get statistics for select pairs

import copy, os
from scitbx.matrix import col
from iotbx.data_manager import DataManager

def modify_info(atom, timepoint):
  # Special cases here
  info = copy.deepcopy(atom)
  if timepoint in ['0F', '1F']:
    altid=''
    info.append(altid)
    if 'OEZ' in atom:
      info[1] = 'OEC'
  if timepoint in ['2F']:
    altid='B'
    info.append(altid)
    if 'OEZ' in atom:
      info[1] = 'OEI'
  if timepoint in ['3F']:
    altid='B'
    info.append(altid)
    if 'OEZ' in atom:
      info[1] = 'OEC'
  if timepoint in ['2F50', '2F500',]:
    altid='B'
    info.append(altid)
  if timepoint in ['2F250', '2F730', '2F1200', '2F2000', '2F4000',]:
    altid='C'
    info.append(altid)
  return info
  


# Main function here

def get_distances_from_endrapid(atom1, atom2, timepoints):
  distances = {}
  for timepoint in timepoints:
    distance_list_for_timepoint = []
    # Modify atom1, atom2 based on some additional information I am providing here
    info1 = modify_info(atom1, timepoint)
    info2 = modify_info(atom2, timepoint)

    # Note that 0 should always be the deposited PDB file i.e non-endrapid PDB
    for i in range(0, 100, 1):
      f = os.path.join('/pscratch/sd/a/asmit/end_rapid_highres_S3/bootstrap_rev2/analysis/%s_endrapid_0624_12_%d.pdb'%(timepoint, i))
      dm = DataManager()
      dm.set_overwrite(True)
      m=dm.get_model(filename=f)
      pdb_hierarchy = m.get_hierarchy()
      for chain in pdb_hierarchy.only_model().chains():
        for residue_group in chain.residue_groups():
          for atom_group in residue_group.atom_groups():
            for atom in atom_group.atoms():
              chain_name = chain.id
              resn = atom_group.resname
              altid = atom_group.altloc
              resi = residue_group.resseq_as_int()
              atom_name = atom.name.strip()
              x,y,z = atom.xyz
              if chain_name == info1[3] and resi==int(info1[2]) and atom_name==info1[0] and resn==info1[1] and altid==info1[4]:
                r1 = col((x,y,z))
              if chain_name == info2[3] and resi==int(info2[2]) and atom_name==info2[0] and resn==info2[1] and altid==info2[4]:
                r2 = col((x,y,z))
      dist = (r1-r2) .length()
      distance_list_for_timepoint.append(dist)
    distances[timepoint] = distance_list_for_timepoint
  from libtbx.easy_pickle import dump
  pickle_fname = 'endrapid_alldists_'+atom1[0]+'_'+atom1[1]+'_'+atom1[2]+'_'+atom2[0]+'_'+atom2[1]+'_'+atom2[2]+'.pkl'
  dump(pickle_fname, distances)

if __name__=='__main__':
  import sys
  #atoms1 = [
  #          ['MN1', 'OEZ', '601', 'A'],
  #         ]
  #atoms2 = [
  #          ['MN4', 'OEZ', '601', 'A'],
  #         ]
  timepoints = ['2F']
  assert len(sys.argv)==9, 'Should have 9 input arguments otherwise something is wrong'
  get_distances_from_endrapid(sys.argv[1:5], sys.argv[5:9], timepoints) 
