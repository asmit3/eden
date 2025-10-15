from __future__ import division, print_function
import sys, time, math, os
from iotbx import reflection_file_reader
from libtbx.utils import to_str
from iotbx.data_manager import DataManager
from scitbx.array_family import flex
import numpy as np
import matplotlib.pyplot as plt
import iotbx.pdb
import mmtbx.model
from scitbx import regular_grid_on_unit_sphere
from mpl_toolkits import mplot3d
from libtbx.phil import parse
from libtbx.utils import Sorry
import time


phil_str = """ 
  model_filename = None
    .type = str
    .help = Name of the model file of the resting state or reference state
  mtz_filename = None
    .multiple = True
    .type = str
    .help = Name of the mtz file with difference coefficients. Multiple values allowed
  residue_numbers = None
    .type = ints
    .help = The residue numbers where the map is going to be calculated
  chain_names = None
    .type = str
    .help = Provide a list of chain names separated by a comma that accompany the residue numbers. This needs\
            to be given if a chain other than chain A is used. One-to-one correspondence\
            with the residue numbers is required.
  conformer_ids = None
    .type = str
    .help = Provide a list of conformer ids separated by a comma that accompany residue numbers. If a certain residue has no \
            alternate conformer or you want to select all components then just leave it as blank for that residue. Eg 3rd residue here has no alternate conformer -->conformer_ids=A,A,,A,B \
            If you do not provide an entry then all atoms (of all conformers) will be used for all the residues provided 
  sphere_radius_angstrom = 2.0
    .type = float
    .help = Sphere radius of probe around atom within which the positive/negative amplitudes will be summed
  radius_increment_angstrom = 0.25
    .type = float
    .help = increment the radius of the fictitous sphere starting from 0 up to sphere_radius in increments of radius_increment_angstrom
  sigma_level_display_line = 3.0
    .type = float
    .help = the sigma level dashed line to be displayed to help guide the eye in the plot
  save_png_fname = None
    .type = str
    .help = Name of filename to save image. Will be preprended to .png extension
  use_extreme_map_value = False
    .type = bool
    .help = If you want to print out the min/max of the map values in the sphere instead of mean values
"""


message = """ Difference density tracer
              -------------------------------------------------------------
              Based on the idea published here: https://aca.scitation.org/doi/10.1063/1.5126921
              Example Usage:
              phenix.python difference_density_tracer.py model_filename=5b6v-H_protein_ret_h2o.pdb mtz_filename=1725us_31.mtz residue_numbers=82,85,89,182,212,216
              See entire phil string for more input options
              -------------------------------------------------------------
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

class map_value_calculator(object):
  def __init__(self, args):
    params = params_from_phil(args)
    
    self.model_fn = params.model_filename                    #   Name of model file
    #from IPython import embed; embed(); exit()
    self.map_coeffs_filenames = params.mtz_filename
    self.residue_numbers = params.residue_numbers
    if params.chain_names is not None:
      self.chain_names = params.chain_names.split(',')
      assert len(self.chain_names) == len(self.residue_numbers), 'Mismatch in length of provided chain_names and residue_numbers'
    else:
      self.chain_names = ['A']*len(self.residue_numbers)

    if params.conformer_ids is not None:
      self.conformer_ids = params.conformer_ids.split(',')
      assert len(self.conformer_ids) == len(self.residue_numbers), 'Mismatch in length of provided conformer_ids and residue_numbers'
    else:
      self.conformer_ids = None


    self.sphere_radius = params.sphere_radius_angstrom
    self.radius_increment = params.radius_increment_angstrom
    self.sigma_level_display_line = params.sigma_level_display_line
       
    self.stepsize=0.01
    self.mtzlabel='FoFo,PHFc'
    self.use_extreme_map_value=params.use_extreme_map_value
    #self.plot_title='%s, %s'%(self.mtzlabel, self.map_coeffs_fn)
    #self.plot_save_fname='%s.png'%(params.save_png_fname)
    

  def read_in_model_mtz(self, fastmode=True):
    print ('Reading in model and mtz files')
    self.ma = {}
    if not fastmode:
      dm = DataManager()                    #   Initialize the DataManager and call it dm
      dm.set_overwrite(True)                #   tell the DataManager to overwrite files with the same name
    #
      self.model = dm.get_model(self.model_fn)        #   Deliver model object with model info
      self.unit_cell = self.model.crystal_symmetry().unit_cell()
      for map_coeffs_fn in self.map_coeffs_filenames:
        miller_arrays = reflection_file_reader.any_reflection_file(
                          file_name=to_str(map_coeffs_fn)).as_miller_arrays()
        ma = dm.get_miller_arrays(labels   = [self.mtzlabel],
                                  filename = map_coeffs_fn)
        self.ma[map_coeffs_fn]=ma[0] # Take the first entry please to maintain consistency with fast mode

    else:
      pdb_inp = iotbx.pdb.input(file_name = self.model_fn)
      self.model = mmtbx.model.manager(model_input = pdb_inp)
      self.unit_cell = self.model.crystal_symmetry().unit_cell()
      for map_coeffs_fn in self.map_coeffs_filenames:
        miller_arrays = reflection_file_reader.any_reflection_file(file_name =
        map_coeffs_fn).as_miller_arrays()
        for m_array in miller_arrays:
          if m_array.info().label_string() == self.mtzlabel:
            self.ma[map_coeffs_fn] = m_array

  def get_map_3d(self): 
    print('Doing FFT')
    self.map_3d = {}
    for map_coeffs_fn in self.map_coeffs_filenames:
      fft_map = self.ma[map_coeffs_fn].fft_map(resolution_factor=0.25)
      fft_map.apply_sigma_scaling()
      self.map_3d[map_coeffs_fn] = fft_map.real_map_unpadded()

  def get_map_value_at_point(self, x,y,z, map_filename=None):
    #print ('getting map value at point')
    site_cart = (x,y,z)
    site_frac = self.unit_cell.fractionalize(site_cart)
    map_value = self.map_3d[map_filename].eight_point_interpolation(site_frac)
    return map_value

  def draw_sphere_of_specified_radius(self, radius=2.0, center_of_sphere=(-18.135, -16.404, -6.739), plot=False):
    x0,y0,z0=center_of_sphere
    spherical_grid=regular_grid_on_unit_sphere.rosca(m=6,hemisphere=False)
    x,y,z = spherical_grid.parts()
    x1=x+x0
    y1=y+y0
    z1=z+z0
    l = x1-x0+0.01
    m = y1-y0+0.01
    n = z1-z0+0.01
    r1=1.0 # Initial radius
    r2 = radius
    s=np.sqrt((r2*r2)/(l*l+m*m+n*n)) # THis is the stepsize needed to go from x0 to x2 in the expanded sphere of radius r2
    # Get coordinates here now
    x2 = l*s+x0
    y2 = y0+m*(x2-x0)/l
    z2 = z0+n*(x2-x0)/l
  
  
    if plot:
      fig=plt.figure()
      ax=plt.axes(projection='3d')
      #ax.scatter3D(x1,y1,z1, color='black')
      ax.scatter3D(x2,y2,z2)
      plt.show()
    return (x2,y2,z2)

  def get_map_value_around_spherical_grid(self, center_of_sphere=(-18.135, -16.404, -6.739), threshold_sigma=1.0, map_filename=None,):
    map_value = 0.0
    positive_map_value = flex.double()
    negative_map_value = flex.double()
    list_radius =  list(np.arange(0, self.sphere_radius+self.radius_increment, self.radius_increment))# List of probe radii to use; will calculate map value at all these points
    for radius in list_radius:
      xs,ys,zs=self.draw_sphere_of_specified_radius(radius=radius, center_of_sphere=center_of_sphere)
      for x,y,z in zip(xs,ys,zs):
        map_value = self.get_map_value_at_point(x,y,z, map_filename=map_filename)
        if map_value > 0 and abs(map_value) > threshold_sigma :
          positive_map_value.append(map_value)
        if map_value < 0 and abs(map_value) > threshold_sigma:
          negative_map_value.append(map_value)
    avg_positive_map_value = 0.0
    avg_negative_map_value = 0.0
    if len(positive_map_value) > 0:
      if self.use_extreme_map_value:
        avg_positive_map_value = flex.max(positive_map_value) #  
      else:
        avg_positive_map_value = flex.mean(positive_map_value) #  
    if len(negative_map_value) > 0:
      if self.use_extreme_map_value:
        avg_negative_map_value = flex.min(negative_map_value) #
      else:
        avg_negative_map_value = flex.mean(negative_map_value) #
    return (avg_positive_map_value, avg_negative_map_value) 

  def get_map_values_around_residue(self):
    self.read_in_model_mtz()
    self.get_map_3d()

    #list_of_residue_numbers = self.residue_numbers #[82, 85, 89, 182, 212, 216]
    nrows = len(self.map_coeffs_filenames)
    ncol = 1
    if nrows == 1:
      f,ax2 = plt.subplots(nrows)
      ax = [ax2]
    else:
      f,ax = plt.subplots(nrows)
    f.set_figheight(10)
    f.set_figwidth(10)
    for imap,map_coeffs_fn in enumerate(self.map_coeffs_filenames):
      pos_vals=[]
      neg_vals=[]
      labels=[]
      residue_demarcation = []
      for ii, residue in enumerate(self.residue_numbers):
        start = time.time()
        if self.conformer_ids is not None:
          if self.conformer_ids[ii] != '':
            sele_str='chain %s and resseq %d and altid %s'%(self.chain_names[ii],residue, self.conformer_ids[ii])
          else:
            sele_str='chain %s and resseq %d'%(self.chain_names[ii],residue,)
        else:
          sele_str='chain %s and resseq %d'%(self.chain_names[ii],residue,)
        isel = self.model.selection(sele_str).iselection()
        ph_sel = self.model.select(isel).get_hierarchy()
        atoms_sel = ph_sel.atoms()
        # demarcation stuff
        if ii == 0:
          first_pt=0
        else:
          first_pt = residue_demarcation[-1]
        second_pt = first_pt+len(atoms_sel)
        residue_demarcation.extend([first_pt, second_pt])

        for jj, atom in enumerate(atoms_sel):
          center=atom.xyz
          name=atom.name.strip()
          pos, neg = self.get_map_value_around_spherical_grid(center_of_sphere=center, map_filename=map_coeffs_fn)
          print ('Mean map values at %s are (signed) %.2f, %.2f'%(atom.id_str(), pos, neg))
          pos_vals.append(pos)
          neg_vals.append(neg)
          if jj == math.floor(len(atoms_sel)/2.0):
            labels.append(name)
            bbox_str = '%d_%s'%(residue, self.chain_names[ii])
            if imap == nrows - 1:
              ax[imap].text(len(labels), 0.0, '%s'%bbox_str, bbox={'facecolor':'white', 'alpha':1, 'edgecolor':'black', 'pad':2}, ha='center', va='center')
          else:
            labels.append(name)
        end=time.time()
        print ("Time elapsed for residue %d:"%residue, end - start)

      ax[imap].plot(pos_vals, '-o',color='g', label=map_coeffs_fn)
      ax[imap].plot(neg_vals, '-o',color='r')
      critical_value=3.0
      if max(pos_vals) > critical_value and min(neg_vals) < -critical_value:
        ax[imap].axvspan(first_pt, second_pt, facecolor='b', alpha=0.3)
      # Demarcate in x
      for pt in residue_demarcation:
        ax[imap].axvline(pt, color='k')
      # Demarcate in y
      ax[imap].axhline(0.0, color='k')
      ax[imap].axhline(self.sigma_level_display_line, color='grey', linestyle="--")
      ax[imap].axhline(-self.sigma_level_display_line, color='grey', linestyle="--")
      # Extra axhlines
      #plt.axhline(2.5, color='grey', linestyle="--", alpha=0.3)
      #plt.axhline(-2.5, color='grey', linestyle="--", alpha=0.3)
      #
      if imap == nrows - 1:
        ax[imap].set_xticks(range(len(labels)), labels, fontsize=8)
      else:
        ax[imap].set_xticks([])
      #plt.title(self.plot_title)
      y_min, y_max=ax[imap].get_ylim()
      y_limit = math.ceil(max(abs(y_min), abs(y_max)))
      ax[imap].set_ylim(-y_limit, y_limit)
      ax[imap].legend()
    
    if False:
      plt.savefig(self.plot_save_fname, dpi=300)
    else:
      plt.tight_layout()
      plt.show()

if (__name__ == "__main__"):

  t0 = time.time()
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print (message)
    print (phil_str)
    exit()
  mvc = map_value_calculator(sys.argv[1:])

  mvc.get_map_values_around_residue()

  print("Total time: %-8.4f"%(time.time()-t0))
