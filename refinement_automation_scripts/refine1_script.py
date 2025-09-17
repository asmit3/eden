import os
import sys
from iotbx import pdb
from iotbx.cli_parser import run_program
from phenix.programs import phenix_refine
from phenix.programs.phenix_refine import custom_process_arguments
from phenix.command_line import superpose_pdbs
from iotbx.data_manager import DataManager
from libtbx.program_template import ProgramTemplate

class Program(ProgramTemplate):
  datatypes = ['model', 'phil', 'miller_array', 'restraint']
  data_manager_options = ['model_skip_expand_with_mtrix']

  master_phil_str = """
  input_files {
    eff_file = None
      .type = path
      .help = Path to eff file for phenix.refine
    template_pdb = None
      .type = path
      .help = Template PDB to superpose and extract OEC (A/a 601) after Step 1
  }
  num_chained_steps = 1
    .type = int
    .help = Number of chained steps to run in the workflow
  labels_name = None
    .type = str
    .help = Label of input X-ray data to be used for refinement
  r_free_flags_name = None
    .type = str
    .help = Label of the R-free flag array to be used for refinement
  r_free_flags_test_value = 1
    .type = int
    .help = Integer value identifying test reflections in the R-free flag array
  add_atoms_mode = *all split
  .type = choice
  .help = Add all OEC atoms at once or split into metals (Mn & Ca) and then oxygen
  upper_high_resolution = None
    .type = float
    .help = Upper limit of high_resolution to use. Eg. for paired-refinement between 2.0-1.5A, this would be 1.5A
  lower_high_resolution = None
    .type = float
    .help = Lower limit of high_resolution to use. Eg. for paired-refinement between 2.0-1.5A, this would be 2.0A
  resolution_bin_size = None
    .type = float
    .help = Resolution bin size to use. For example, 0.01 or 0.1A would be good picks
  phenix_refine_params_file = None
    .type = str
    .help = Supply all parameters or basic other refinement stuff in this parameter file. Please do not specify high resolution as that should be done using the phil scope of this current program.
  folder_prefix = None
    .type = str
    .help = Prefix for the paired refinement folders. Default will be dataset_res1_res2 with subfolders dataset_res1 and dataset_res2. Res1 and res2 will be numbers 
  nproc = 1
    .type = int
    .help = Number of processors to use to run the refinement jobs
  output_folder = None
    .type = str
    .help = Output folder where results will be saved. Will be relative to the base_path 
  randomize_model = True
    .type = bool
    .help = Randomize the x,y,z coordinates by an rms (right now 0.25) and also adp. Random seed is hardcoded for now 
  randomize_outerbin = True
    .type = bool
    .help = Randomize the outerbin of each paired dataset as a control. Example, if we do a 2.0 and 2.10 A paired refinement, we will randomize the 2.0-2.1A bin for an additional refinement for the 2.0A case to see how the R-factors look 
  binning_strategy = *resolution volume
    .type = choice
    .help = Select between binning mode of equal resolution bins or equal volume (i.e same # of reflection bins).
  volume_nbins = 10
    .type = int
    .help = Number of bins between upper_high_res and lower_high_res to use when using the binning_strategy of volume
  output_prefix = "refine"
    .type =  str
    .help = Prefix for output files from each refinement step
  """
  
  """Add a waters template pdb above once Step 4 is implemented: waters_template_pdb = None
      .type = path
      .help = Template PDB for waters (Chain G/g) replacement in Step 4"""
  def validate(self):
    self.data_manager.has_models(raise_sorry = True)
    self.data_manager.has_miller_arrays(raise_sorry=True)
    print ('Inputs OK', file=self.logger)

  def run(self):
    num_chained_steps = self.params.num_chained_steps
    refiner = Refiner(self.params, self.data_manager)
    refiner.run_chain(steps=num_chained_steps)

class Refiner(object):
  def __init__(self, params, data_manager):
    self.base_path = os.getcwd()
    self.params = params
    self.data_manager = data_manager
    self.arguments_for_phenixrefine = []
    
    self.sf_filenames = [os.path.join(self.base_path, f) for f in data_manager.get_miller_array_names()]
    for sf in self.sf_filenames:
      self.arguments_for_phenixrefine.append(sf)

    model_names = data_manager.get_model_names()
    template_basename = None
    if self.params.input_files.template_pdb:
      template_basename = os.path.basename(self.params.input_files.template_pdb)

    working_model = None
    for name in model_names:
      if template_basename is None or os.path.basename(name) != template_basename:
        working_model = name
        break
    self.model_filename = os.path.join(self.base_path, working_model)
    self.arguments_for_phenixrefine.append(self.model_filename)

    for restraint_file in data_manager.get_restraint_names():
      self.arguments_for_phenixrefine.append(os.path.join(self.base_path, restraint_file))
    
    if params.input_files.eff_file  is not None:
      self.arguments_for_phenixrefine.append(os.path.join(self.base_path, params.input_files.eff_file))
    if params.labels_name is not None:
      self.arguments_for_phenixrefine.append(f"labels.name={params.labels_name}")
    if params.r_free_flags_name is not None:
      self.arguments_for_phenixrefine.append(f"xray_data.r_free_flags.label={params.r_free_flags_name}")
      self.arguments_for_phenixrefine.append(f"xray_data.r_free_flags.test_flag_value={params.r_free_flags_test_value}")

    self.arguments_for_phenixrefine.append(f"output.prefix={params.output_prefix}")
    #REMOVE THIS LINE AT THE END:
    self.arguments_for_phenixrefine.append("refinement.main.number_of_macro_cycles=1")
    self.arguments_for_phenixrefine.append("overwrite=true")

  def run_refine(self):
    run_program(
        program_class=phenix_refine.Program,
        args=self.arguments_for_phenixrefine,
        custom_process_arguments=custom_process_arguments,
        logger=None
    )

  def get_latest_pdb(self):
    pdb_files = [f for f in os.listdir() if f.startswith("refine_") and f.endswith(".pdb")]
    pdb_files.sort(key=lambda f: os.path.getmtime(f))
    return pdb_files[-1]
  
  def add_OEC_atoms(self, step1_file, aligned_template_file, output_pdb, element_selection="all"):
    dm = DataManager()
    dm.set_overwrite(True)
    model_step1 = dm.get_model(step1_file)
    model_template = dm.get_model(aligned_template_file)
    ph1 = model_step1.get_hierarchy()
    ph2 = model_template.get_hierarchy()

    allowed_elements = {
      "all": {"MN", "CA", "O"},
      "metals": {"MN", "CA"},
      "oxygen": {"O"}
    }
    selected_elements = allowed_elements[element_selection]
    
    for chain in ph2.chains():
      if chain.id.strip() in ["A", "a"]:
        for rg in chain.residue_groups():
          if rg.resid().strip() == "601":
            rg_copy = rg.detached_copy()
            for atom_group in rg_copy.atom_groups():
              delete_atoms = [atom for atom in atom_group.atoms() if atom.element.upper() not in selected_elements]
              for atom in delete_atoms:
                atom_group.remove_atom(atom)
            for target_chain in ph1.chains():
              if target_chain.id == chain.id:
                target_chain.append_residue_group(rg_copy)
                break

    
    dm.write_model_file(model_step1, output_pdb)

  def add_solvent(self, input_pdb, output_pdb):
    args = []
    replaced_pdb = False
    for arg in self.arguments_for_phenixrefine:
      if arg.endswith(".pdb"):
        if not replaced_pdb:
          args.append(input_pdb)
          replaced_pdb = True
      else:
        args.append(arg)

    args.extend([
        f"output.prefix={os.path.splitext(output_pdb)[0]}",
        "refinement.main.number_of_macro_cycles=1",
        "refinement.refine.strategy=individual_sites",
        "main.ordered_solvent=true"
    ])

    run_program(
        program_class=phenix_refine.Program,
        args=args,
        custom_process_arguments=custom_process_arguments,
        logger=None
    )

  """def change_waters(self, main_pdb, template_pdb, output_pdb, threshold = 1):
    dm = DataManager()
    dm.set_overwrite(True)

    model_main = dm.get_model(main_pdb)
    model_template = dm.get_model(template_pdb)
    main = model_main.get_hierarchy()
    template = model_template.get_hierarchy()
    
    chain_S = [c for c in main.chains() if c.id.strip() == "S"][0]
    chain_Gg = [c for c in template.chains() if c.id.strip() in ["G", "g"]]
    
    waters_to_remove = []

  for chain in chain_Gg:
      for rg_template in chain.residue_groups():
        for atom_group_template in rg_template.atom_groups():
          for atom_template in atom_group_template.atoms():
            matched = False
            for rg_S in chain_S.residue_groups():
              for ag_S in rg_S.atom_groups():
                for atom_S in ag_S.atoms():
                #threshold: use 1 Angstrom for now
                  if atom_S.element.strip().upper() == "O" and (atom_S.xyz == atom_template.xyz):
                    atom_template.set_xyz(atom_S.xyz)
                    atom_template.set_b(atom_S.b)
                    atom_template.set_occ(atom_S.occ)
                    waters_to_remove.append(atom_S)
                    matched = True
                    break
                if matched: break
              if matched: break
    for atom in waters_to_remove:
      atom.parent().remove_atom(atom)
    for chain in chain_Gg:
      main.append_chain(chain.detached_copy())

    dm.write_model_file(model_main, output_pdb)"""
    
  def run_chain(self, steps=3):
    for i in range(steps):
      if i == 1 and self.params.input_files.template_pdb:
        step1_output = os.path.join(self.base_path, self.get_latest_pdb())
        template_pdb = os.path.join(self.base_path, self.params.input_files.template_pdb)
        superposed_path = os.path.join(self.base_path, "superposed.pdb")
        step2_input = os.path.join(self.base_path, "step2_input.pdb")

        var1 = f"input.pdb_file_name_fixed={step1_output}"
        var2 = f"input.pdb_file_name_moving={template_pdb}"
        var3 = f"output.file_name={superposed_path}"
        print(var1)
        print(var2)
        print(var3)
        superpose_pdbs.run([var1, var2, var3])
        
        if self.params.add_atoms_mode == "all":
          self.add_OEC_atoms(step1_output, superposed_path, step2_input, element_selection="all")
          self.model_filename = step2_input
          for j, arg in enumerate(self.arguments_for_phenixrefine):
            if arg.endswith(".pdb"):
              self.arguments_for_phenixrefine[j] = self.model_filename
              break
          self.run_refine()

        elif self.params.add_atoms_mode == "split":
          #Step 2a: Adding manganese and calcium
          step2a_input = os.path.join(self.base_path, "step2a_input.pdb")
          self.add_OEC_atoms(step1_output, superposed_path, step2a_input, element_selection="metals")
          self.model_filename = step2a_input
          for j, arg in enumerate(self.arguments_for_phenixrefine):
            if arg.endswith(".pdb"):
              self.arguments_for_phenixrefine[j] = self.model_filename
              break
          self.run_refine()

          #Step 2b: Adding oxygens on top of refined file from Step 2a
          step2b_input = os.path.join(self.base_path, "step2b_input.pdb")
          pdb_after2a = os.path.join(self.base_path, self.get_latest_pdb())
          self.add_OEC_atoms(pdb_after2a, superposed_path, step2b_input, element_selection="oxygen")
          self.model_filename = step2b_input
          for j, arg in enumerate(self.arguments_for_phenixrefine):
            if arg.endswith(".pdb"):
              self.arguments_for_phenixrefine[j] = self.model_filename
              break
          self.run_refine()

      elif i == 2:
        step2_output = os.path.join(self.base_path, self.get_latest_pdb())
        step3_input = os.path.join(self.base_path, "step3_input.pdb")
        self.add_solvent(step2_output, step3_input)
        self.model_filename = step3_input
        for j, arg in enumerate(self.arguments_for_phenixrefine):
          if arg.endswith(".pdb"):
            self.arguments_for_phenixrefine[j] = self.model_filename
            break
        self.run_refine()

        """elif i == 3 and self.params.input_files.waters_template_pdb:
        #step3_output = os.path.join(self.base_path, self.get_latest_pdb())
        #step4_input = os.path.join(self.base_path, "step4_input.pdb")
        aligned_template = os.path.join(self.base_path, "aligned_template.pdb")
        superpose_pdbs.run(
        fixed_pdb_file=step3_output,
        moving_pdb_file=self.template_pdb,
        output_pdb_file=aligned_template,
        log=sys.stdout,
        )
        #waters_template = os.path.join(self.base_path, self.params.input_files.waters_template_pdb)
        #self.change_waters(step3_output, waters_template, step4_input)
        #self.model_filename = step4_input
        #for j, arg in enumerate(self.arguments_for_phenixrefine):
          #if arg.endswith(".pdb"):
            #self.arguments_for_phenixrefine[j] = self.model_filename
            #break
        #self.run_refine()"""

      elif i == 0:
        self.run_refine()

      latest_pdb = os.path.join(self.base_path, self.get_latest_pdb())
      self.data_manager.process_model_file(latest_pdb)
      self.model_filename = latest_pdb
      for j, arg in enumerate(self.arguments_for_phenixrefine):
        if arg.endswith(".pdb"):
          self.arguments_for_phenixrefine[j] = self.model_filename
          break

def run_main():
  run_program(program_class=Program)
                
if __name__ == "__main__":
  run_main()