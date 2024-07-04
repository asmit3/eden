from iotbx.cli_parser import run_program
from phenix.programs import phenix_refine
from phenix.programs.phenix_refine import custom_process_arguments
from libtbx.phil import parse
import numpy as np
import os
import copy
import glob
import mmtbx, iotbx
from iotbx import reflection_file_reader
from scitbx.array_family import flex

phil_str = """
  model_filename = None
    .type = str
    .help = Name of the model file. Please provide the path
  sf_filename = None
    .type = str
    .help = Name of mtz/sf-cif filename
  upper_high_resolution = None
    .type = float
    .help = Upper limit of high_resolution to use. Eg. for paired-refinement between 2.0-1.5A, this would be 1.5A
  lower_high_resolution = None
    .type = float
    .help = Lower limit of high_resolution to use. Eg. for paired-refinement between 2.0-1.5A, this would be 2.0A
  resolution_bin_size = None
    .type = float
    .help = Resolution bin size to use. For example, 0.01 or 0.1A would be good picks
  restraint_files = None
    .type = str
    .help = Supply restraint files for refinement
    .multiple = True 
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
    .help = randomize the x,y,z coordinates by an rms (right now 0.25) and also adp. Random seed is hardcoded for now. 
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


class paired_refiner(object):

  def __init__(self, args):
    self.base_path = os.getcwd()
    params = params_from_phil(args)
    assert params.model_filename is not None, 'Please supply model filename'
    assert params.sf_filename is not None, 'Please supply structure factor file'
    #assert params.phenix_refine_params_file, 'Please supply basic phenix_refine parameter file'
    self.model_filename=os.path.join(self.base_path, params.model_filename)
    self.sf_filename=os.path.join(self.base_path, params.sf_filename)
    self.arguments_for_phenixrefine = [self.model_filename, self.sf_filename,]
    for restraint_file in params.restraint_files:
      self.arguments_for_phenixrefine.append(os.path.join(self.base_path, restraint_file))
    if params.phenix_refine_params_file is not None:
      self.arguments_for_phenixrefine.append(os.path.join(self.base_path, params.phenix_refine_params_file))
    self.randomize_model=params.randomize_model
    if self.randomize_model:
      randomization_str = ['refinement.modify_start_model.modify.sites.shake=0.25', 'refinement.modify_start_model.modify.sites.atom_selection=all', 'refinement.modify_start_model.modify.adp.randomize=True']
      self.arguments_for_phenixrefine.extend(randomization_str)

    self.upper_high_resolution = params.upper_high_resolution
    self.lower_high_resolution = params.lower_high_resolution
    self.resolution_bin_size = params.resolution_bin_size
    self.nproc = params.nproc
    if params.folder_prefix is not None:
      self.folder_prefix = params.folder_prefix
    else:
      self.folder_prefix = 'dataset' 
    if params.output_folder is None:
      self.output_folder = '.'
    else:
      self.output_folder = params.output_folder

    # Now establish the resolutions at which refinements should be done
    self.bins = np.arange(self.upper_high_resolution, self.lower_high_resolution+self.resolution_bin_size, self.resolution_bin_size)
    self.n_subfolders=3 # Hardcoding number of subfolders t1,t2,t3 for now

  def make_folders(self, do_not_create=False):
    number_of_folders = len(self.bins) - 1
    self.all_folders = [] # List of all folders where refinement is done
    self.all_parent_folders = []
    for i in range(number_of_folders):
      self.parent_folder = os.path.join(self.base_path, self.output_folder,  '%s_pair_%.2f_%.2f'%(self.folder_prefix, self.bins[i], self.bins[i+1]))
      self.child_folder = os.path.join(self.parent_folder, '%s_%.2f'%(self.folder_prefix, self.bins[i]))
      self.child_folder2 = os.path.join(self.parent_folder, '%s_%.2f'%(self.folder_prefix, self.bins[i+1]))
      if not do_not_create:
        os.makedirs(self.parent_folder, exist_ok=False)
        os.makedirs(self.child_folder, exist_ok=False) 
        os.makedirs(self.child_folder2, exist_ok=False) 
        for j in range(self.n_subfolders):
          os.makedirs(os.path.join(self.child_folder, 't%d'%(j+1)))
          os.makedirs(os.path.join(self.child_folder2, 't%d'%(j+1)))
      self.all_parent_folders.append(self.parent_folder)
      self.all_folders.append(self.child_folder)
      self.all_folders.append(self.child_folder2)

  def run_refinement(self, folder, trial_folder, args):
    os.chdir(os.path.join(self.output_folder, folder, trial_folder))
    print ('Starting Paired refinement in folder %s'%folder)
    r = run_program(
      program_class = phenix_refine.Program,
      args          = args,
      custom_process_arguments = custom_process_arguments,
      logger        = None # use default
      )
    print(dir(r))
    rwork = round(r.fmodel.r_work()*100, 1)
    rfree = round(r.fmodel.r_free()*100, 1)
    print(rwork,rfree)
    print ('Finished Paired refinement in folder %s'%folder)


  def run_paired_refinement(self):
    # Create args for multicore run
    argstuples = []
    for ii, folder in enumerate(self.all_folders):
      os.chdir(os.path.join(self.output_folder, folder))
      high_resolution = float(folder.split('_')[-1])
      args = copy.deepcopy(self.arguments_for_phenixrefine)
      #self.run_refinement(args+['xray_data.high_resolution=%.2f'%high_resolution])
      for jj in range(self.n_subfolders):
        trial_folder='t%d'%(jj+1)
        if self.randomize_model:
          args += ['main.random_seed=%d'%((ii+1)*self.n_subfolders+jj+22)] # Offset of 22 just for fun
        argstuples.append((folder,trial_folder,  args+['xray_data.high_resolution=%.2f'%high_resolution]))
    from libtbx import easy_mp
    for args, res, errstr in easy_mp.multi_core_run(self.run_refinement, argstuples, self.nproc):
      print ("arguments: %s \nresult: %s \nerror: %s\n" %(args, res, errstr))
    print ('Done with all paired refinement')
    os.chdir(self.base_path)

  def get_rfactor_statistics(self):
    table = []
    for folder in self.all_parent_folders:
      print ('Getting R-factor from %s'%folder)
      res1=float(folder.split('_')[-2])
      res2=float(folder.split('_')[-1])
      rwork1, rfree1 = self.get_rfactors_from_single_refinement(lower_high_resolution=res2, folder = os.path.join(folder, self.folder_prefix+'_%.2f'%res1))
      rwork2, rfree2 = self.get_rfactors_from_single_refinement(lower_high_resolution=res2, folder = os.path.join(folder, self.folder_prefix+'_%.2f'%res2))

      text = ['%.2f - %.2f'%(res1, res2), rwork1, rwork2, rwork1-rwork2,rfree1, rfree2, rfree1-rfree2]
      table.append(text)
    from tabulate import tabulate
    print (tabulate(table, headers=['Resolution bins', 'Rwork(higher)', 'Rwork(lower)', 'deltaR(high-low)','Rfree(higher)', 'Rfree(lower)', 'deltaRfree(high-low)', ], tablefmt='orgtbl'))

  def get_rfactors_from_single_refinement(self, lower_high_resolution=None, folder=None,):
      # First read in data from lower resolution refinement
      all_rwork = flex.double()
      all_rfree = flex.double()
      for ii in range(self.n_subfolders):
        trial_folder = 't%d'%(ii+1)
        os.chdir(os.path.join(self.output_folder, folder, trial_folder))
        pdb_files = glob.glob( '*.pdb') # XXX This assumes only pdb file format
        assert len(pdb_files) == 1
        phenix_refine_prefix = pdb_files[0].split('.pdb')[0]
        mtz_file = phenix_refine_prefix+'.mtz' 
        #print ('PDB and MTZ file = ', pdb_files[0], mtz_file)
        # Now read in fmodel
        pdb_inp = iotbx.pdb.input(file_name = pdb_files[0])
        model = mmtbx.model.manager(model_input = pdb_inp)
        miller_arrays = reflection_file_reader.any_reflection_file(file_name =
        mtz_file).as_miller_arrays()
        for m_array in miller_arrays:
          if m_array.info().label_string() == 'F-obs-filtered,SIGF-obs-filtered':
            f_obs = m_array#.resolution_filter(d_min=1.95)
            #f_obs = m_array.as_amplitude_array()#.resolution_filter(d_min=1.95) <-- Only needed for I-obs
          if m_array.info().label_string() == 'R-free-flags':
            r_free_flags = m_array
        f_obs, r_free_flags = f_obs.common_sets(r_free_flags)
        r_free_flags = r_free_flags.array(data = r_free_flags.data()==0)
        #print(r_free_flags.data().count(True), r_free_flags.data().count(False))
        fmodel = mmtbx.f_model.manager(f_obs=f_obs, r_free_flags=r_free_flags, xray_structure=model.get_xray_structure())
        fmodel = fmodel.resolution_filter(d_min=lower_high_resolution) # Cutting data here to right resolution
        fmodel.update_all_scales()
        rwork = fmodel.r_work()
        rfree = fmodel.r_free()
        all_rwork.append(rwork)
        all_rfree.append(rfree)
      os.chdir(self.base_path)
      avg_rwork = flex.mean(all_rwork)
      avg_rfree = flex.mean(all_rfree)
      return (avg_rwork, avg_rfree)

if __name__ == '__main__':
  import sys
  pr = paired_refiner(sys.argv[1:]) 
  pr.make_folders(do_not_create=False)
  pr.run_paired_refinement()
  pr.get_rfactor_statistics()
