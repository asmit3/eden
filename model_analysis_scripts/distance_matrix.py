from phenix.command_line import superpose_pdbs
from iotbx.data_manager import DataManager
from libtbx.program_template import ProgramTemplate
import scitbx.math
import numpy as np

class Program(ProgramTemplate):
  datatypes = ['model', 'phil']
  data_manager_options = ['model_skip_expand_with_mtrix']

  master_phil_str = """
    input {
      selection_str = None
        .type = str
        .help = Common selection string over which to calculate distance matrix. Otherwise assumed to be all\
                If it is all, will check that the same Calpha atoms are there If you are providing a selection syntax, try to use Calpha atoms if possible
  }
  """
  def validate(self):
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 2,
      exact_count = True)
    m = self.data_manager.get_model()

  def run(self):
    print ('Running the Distance Matrix program through ProgramTemplate')
    dmatrix = DistanceMatrix(self.params, self.data_manager)
    dmatrix.run()

class DistanceMatrix(object):
  def __init__(self, params, data_manager):
    self.model_1 = data_manager.get_model_names()[0]
    self.model_2 = data_manager.get_model_names()[1]
    self.modelname_1 = data_manager.get_model_names()[0] 
    self.modelname_2 = data_manager.get_model_names()[1] 
    if params.input.selection_str is not None:
      self.selection_str = params.input.selection_str
    else:
      self.selection_str = 'name CA'

  def superpose_models(self):
    var1 = 'input.pdb_file_name_fixed = %s'%self.model_1
    var2 = 'input.pdb_file_name_moving = %s'%self.model_2
    var3 = 'output.file_name = %s_superposed.pdb'%self.model_2
    superpose_pdbs.run([var1, var2, var3])
    self.model_2 = '%s_superposed.pdb'%self.model_2
      
  def read_in_models(self):
      dm1 = DataManager()                    #   Initialize the DataManager and call it dm
      dm1.set_overwrite(True)                #   tell the DataManager to overwrite files with the same name
      self.model_1 = dm1.get_model(self.model_1)
      #
      dm2 = DataManager()                    #   Initialize the DataManager and call it dm
      dm2.set_overwrite(True)                #   tell the DataManager to overwrite files with the same name
      self.model_2 = dm2.get_model(self.model_2)
      #self.cs_2 = self.model_2.crystal_symmetry()

  def get_atom_selection(self):
    isel1 = self.model_1.selection(self.selection_str).iselection()
    ph_sel_1 = self.model_1.select(isel1).get_hierarchy()
    self.atoms_sel_1 = ph_sel_1.atoms()
    self.coordinates_1=self.atoms_sel_1.extract_xyz()
    #
    isel2 = self.model_2.selection(self.selection_str).iselection()
    ph_sel_2 = self.model_2.select(isel2).get_hierarchy()
    self.atoms_sel_2 = ph_sel_2.atoms()
    self.coordinates_2=self.atoms_sel_2.extract_xyz()
    # Some checks are needed here
    print ('Length of number of atoms for distance matrix 1= %d'%(len(self.coordinates_1)))
    print ('Length of number of atoms for distance matrix 2= %d'%(len(self.coordinates_2)))
    assert len(self.coordinates_1) == len(self.coordinates_2), 'Number of atoms in selection has to be same for both models'
    assert ph_sel_1.is_ca_only() and ph_sel_2.is_ca_only(), 'All atoms in selection need to be C-alpha only'
    # Figure out boundaries of different chains
    self.unique_ids = []
    self.chain_ids = ph_sel_1.chain_ids()
    for ch in ph_sel_1.chains():
      chain_id = ch.id.strip()
      resid_group=ch.residue_groups()
      for res in resid_group:
        res_id = res.resid().strip()
        unique_id = '%s_%s'%(chain_id, res_id) 
        print (unique_id)
        self.unique_ids.append(unique_id)
  
  def get_distance_matrix(self):
    self.distance_matrix = scitbx.math.distance_difference_matrix(self.coordinates_1, self.coordinates_2)

  def get_xtick_values(self):
    # Find boundaries between chains first
    xvals = []
    xtick_vals = []
    chain_ids_seen = []
    chain_boundary_idx = []
    xtick_interval = 5
    for ii, unique_id in enumerate(self.unique_ids):
      chain_id = unique_id.split('_')[0]
      if chain_id not in chain_ids_seen:
        xvals.append(ii)
        xtick_vals.append(unique_id)
        chain_ids_seen.append(chain_id)
        chain_boundary_idx.append(ii)
      elif ii%xtick_interval == 0:
        xvals.append(ii)
        xtick_vals.append(unique_id)
    return (xvals, xtick_vals, chain_boundary_idx)

  def draw_plot(self):
    import matplotlib.pyplot as plt
    extremum_value = max(abs(self.distance_matrix))+0.05
    plt.imshow(self.distance_matrix.as_numpy_array(), cmap='bwr', origin='lower', vmin=-extremum_value, vmax=extremum_value) 
    cbar=plt.colorbar(label='Distance Matrix Value (Angstroms)', fraction=0.05,)
    xvals, xtick_custom_vals, chain_boundary_idx = self.get_xtick_values()
    plt.title('Calpha-Calpha distance matrix')
    plt.xlabel('Residue ID, model 1 (%s)'%self.modelname_1)
    plt.ylabel('Residue ID, model 2 (%s)'%self.modelname_2)
    plt.xticks(xvals, xtick_custom_vals, rotation=90) 
    plt.yticks(xvals, xtick_custom_vals, rotation=90) 
    linewidth=0.5
    for boundary_idx in chain_boundary_idx:
      #print (boundary_idx)
      plt.plot((0, boundary_idx), (boundary_idx, boundary_idx), 'k-', linewidth=linewidth)
      plt.plot((boundary_idx, boundary_idx), (0, boundary_idx), 'k-', linewidth=linewidth)
    plt.show()

  def run(self):
    self.superpose_models()
    self.read_in_models()
    self.get_atom_selection()
    self.get_distance_matrix()
    self.draw_plot()

if __name__ == '__main__':
  import sys
  from iotbx.cli_parser import run_program
  """ Distance Matrix Calculator Program
      Provide 2 model files (with same sequence) and give a selection string (optional)
      Example Usage:
      libtbx.python distance_matrix.py 3M1V.cif 7SUC.cif input.selection_str="chain A and resseq 80:100 and name CA"
  """
  run_program(program_class=Program)
