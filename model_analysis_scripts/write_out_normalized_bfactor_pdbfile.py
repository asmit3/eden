# Please make sure you have sourced phenix python before running this script

def run(arg):
  from iotbx.data_manager import DataManager
  from scitbx.array_family import flex
  dm = DataManager()
  dm.set_overwrite(True)
  fn = arg

  model = dm.get_model(fn)
  ph = model.get_hierarchy()
  atoms = ph.atoms()
  bfactors = atoms.extract_b()
  stats = flex.mean_and_variance(bfactors)
  mean = stats.mean()
  std = stats.unweighted_sample_standard_deviation()
  normalized_b = (bfactors - mean)/std 
  atoms.set_b(normalized_b)
  #
  # Limited selection to write out chain G/g statistics
  sele_str = 'chain G or chain g'
  isel=model.selection(sele_str).iselection()
  ph_sel = model.select(isel).get_hierarchy()
  atoms_sel = ph_sel.atoms()
  ph_sel.write_pdb_file('%s_only_chainGg.pdb'%(arg[:-4]))
  print ('Wrote out PDB file with normalized bfactor for selection: %s: %s_only_chainGg.pdb'%(sele_str, arg[:-4]))
  if False:
    bf=atoms_sel.extract_b()
    import matplotlib.pyplot as plt
    plt.hist(bf)
    plt.show()
  
  # Normalize bfactors here
  ph.write_pdb_file('%s_normalized_bfactors.pdb'%(arg[:-4]))
  print ('Wrote out PDB file with normalized bfactor: %s_normalized_bfactors.pdb'%(arg[:-4]))
  print ('Display in pymol for example: spectrum b, blue_white_red, minimum=-5, maximum=5')

if __name__ == '__main__':
  import sys
  message = """
Script to write out a PDB file with normalized b-factors. Run the script as following
cctbx.python write_out_normalized_bfactor_pdbfile.py <pdb_filename> 
"""
  print (message)
  run(sys.argv[1])
