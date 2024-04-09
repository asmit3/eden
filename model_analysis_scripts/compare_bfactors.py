from libtbx.phil import parse
from libtbx.utils import Sorry
import os, sys

phil_str = """
  pdb_1 = None
    .type = str
    .help = provide pdb file name 1
  pdb_2 = None
    .type = str
    .help = provide pdb file name 2
  chains = None
    .type = str
    .help = provide chains to use for plotting. Separated by comma. Eg. A,B,C
    .multiple=True
"""
phil_scope = parse(phil_str)

def params_from_phil(args, phil_scope=phil_scope):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      user_phil.append(parse(file_name=arg))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  return params


def run(args):
  from iotbx.data_manager import DataManager
  from scitbx.array_family import flex
  dm = DataManager()
  dm.set_overwrite(True)
  
  params = params_from_phil(args)
  fn1 = params.pdb_1
  fn2 = params.pdb_2
  chains_to_plot = []
  if len(params.chains) >0:
    for pch in params.chains:
      chains_to_plot.extend(pch.split(','))
  model_1 = dm.get_model(fn1)
  model_2 = dm.get_model(fn2)
  ph1 = model_1.get_hierarchy()
  ph2 = model_2.get_hierarchy()
  atom_selection_str = 'name CA'
  if len(chains_to_plot) > 0:
    atom_selection_str += ' and ('
    for ii, chain in enumerate(chains_to_plot):
      if ii != len(chains_to_plot)-1:
        atom_selection_str+='chain %s or '%chain
      else:
        atom_selection_str+='chain %s)'%chain
  print ('Phenix atom selection syntax:%s'%atom_selection_str)
  ph1_ca = ph1.apply_atom_selection(atom_selection_str)
  ph2_ca = ph2.apply_atom_selection(atom_selection_str)
  atoms_1 = ph1_ca.atoms()
  atoms_2 = ph2_ca.atoms()
  bfactors_1 = flex.double()
  bfactors_2 = flex.double()
  residue_info = flex.std_string()
  chain_info = flex.std_string()
  import matplotlib.colors as mcolors
  css4_dict = mcolors.CSS4_COLORS

  for chain in ph1_ca.chains():
    for chain_2 in ph2_ca.chains():
      if chain.id == chain_2.id:
        for residue_group in chain.residue_groups():
          for residue_group_2 in chain_2.residue_groups():
            if residue_group.resid() == residue_group_2.resid():
              #from IPython import embed; embed(); exit()
              if residue_group.have_conformers() or residue_group_2.have_conformers():continue
              #print (residue_group.id_str(), residue_group_2.id_str(), residue_group.atoms().extract_b().all(), residue_group_2.atoms().extract_b().all())
              assert residue_group.atoms().extract_b().all() == residue_group_2.atoms().extract_b().all(), 'Number of atoms in residue groups does not match'
              bfactors_1.append(residue_group.atoms().extract_b()[0])
              bfactors_2.append(residue_group_2.atoms().extract_b()[0])
              residue_info.append(residue_group.id_str())
              chain_info.append(chain.id)
 
  stats_1 = flex.mean_and_variance(bfactors_1)
  stats_2 = flex.mean_and_variance(bfactors_2)
  mean_1 = stats_1.mean()
  mean_2 = stats_2.mean()
  std_1 = stats_1.unweighted_sample_standard_deviation()
  std_2 = stats_2.unweighted_sample_standard_deviation()

  norm_b1 = (bfactors_1 - mean_1)/std_1
  norm_b2 = (bfactors_2 - mean_2)/std_2

  import matplotlib.pyplot as plt
  #plt.rcParams['axes.facecolor'] = 'lightgrey'
  fig = plt.figure()
  unique_chain_id = list(set(chain_info))
  css4_color_keys = list(css4_dict.keys())

  # Christina - please add a unique color for each chain below. Ue color schemes from here - 
  # https://matplotlib.org/stable/gallery/color/named_colors.html
  #
  #css4_dict = {'A':'red', 'B':'green', 'a': 'blue', 'b':'black'}

  color_dict = {}
  for ii, unique_chain in enumerate(unique_chain_id):
    color_dict[unique_chain] = css4_color_keys[ii]

  color_list = []
  for chain in chain_info:
    color_list.append(color_dict[chain])
  if len(chains_to_plot) == 0:
    plt.scatter(norm_b1, norm_b2, label='All chains') #,c=color_list)
  else:
    for chain in unique_chain_id:
      plt.scatter(norm_b1.select(chain_info==chain), norm_b2.select(chain_info==chain),c=color_dict[chain], label='chain %s'%chain)

  plt.xlabel('Dataset %s'%fn1)
  plt.ylabel('Dataset %s'%fn2)
  plt.legend()
  plt.title('Normalized B-factors by stdev, (B-Bavg)/sigmaB')
  plt.show()


if __name__ == '__main__':
  import sys
  message = """ 
Run script as cctbx.python compare_bfactors.py pdb_1=<pdb_file_1> pdb_2=<pdb_file_2> chains=<chain ids>
Make sure you have sourced phenix python on the command line prior to running script
chain ids if provided should be separated by commas
Script plots the normalized B-factor between pdb1 and pdb2 for the common set of Calpha atoms
Example usage:
 - cctbx.python compare_bfactors.py file1.pdb file2.pdb 
 - cctbx.python compare_bfactors.py file1.pdb file2.pdb chains=A,a 
            """
  print (message)
  
  pdb_f1, pdb_f2=sys.argv[1],sys.argv[2]
  run(sys.argv[1:])
