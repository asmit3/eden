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
from libtbx.utils import Sorry
from scitbx import regular_grid_on_unit_sphere
from scitbx.matrix import col
from multiprocessing import Pool

class map_value_calculator(object):
  def __init__(self, args):
    self.state=args[0]
    # Sort out mode here whether fcalc or fobs
    if len(args)==2:
      self.mode = 'fobs' # Only 2 options, fobs or fcalc 
    if len(args) > 2:
      if args[2] == 'fobs':
        self.mode = 'fobs' # Only 2 options, fobs or fcalc 
      if args[2] == 'fcalc':
        self.mode = 'fcalc'
      if args[2] == 'fofc':
        self.mode = 'fofc'
      if args[2] == 'anom':
        self.mode = 'anom'
      else:
        Sorry('2nd argument needs to be either fobs, fcalc or just leave it blank')
    # Now go !
    self.dataset = args[1]
    self.model_fn = os.path.join(args[0]+'.pdb')                    #   Name of model file
    if self.mode == 'fobs':
      self.map_coeffs_fn = args[1]+'.mtz'         # Might or might not have fcalc 
      #self.mtzlabel='FOFCWT,PHFOFCWT'
      self.mtzlabel='2FOFCWT_no_fill,PH2FOFCWT_no_fill'
    self.stepsize=0.01
    if self.mode == 'fcalc':
      self.map_coeffs_fn = args[1]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='F-model,PHIF-model'
    if self.mode == 'fofc':
      self.map_coeffs_fn = args[1]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='FOFCWT,PHFOFCWT'
    if self.mode == 'anom':
      self.map_coeffs_fn = args[1]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='ANOM,PANOM'
      

  def read_in_model_mtz(self, fastmode=True):
    print ('Reading in model and mtz files')
    if not fastmode:
      dm = DataManager()                    #   Initialize the DataManager and call it dm
      dm.set_overwrite(True)                #   tell the DataManager to overwrite files with the same name
    #
      self.model = dm.get_model(self.model_fn)        #   Deliver model object with model info
      self.unit_cell = self.model.crystal_symmetry().unit_cell()
      miller_arrays = reflection_file_reader.any_reflection_file(
                        file_name=to_str(self.map_coeffs_fn)).as_miller_arrays()
      ma = dm.get_miller_arrays(labels   = [self.mtzlabel],
                                filename = self.map_coeffs_fn)
      self.ma=ma[0] # Take the first entry please to maintain consistency with fast mode

    else:
      pdb_inp = iotbx.pdb.input(file_name = self.model_fn)
      self.model = mmtbx.model.manager(model_input = pdb_inp)
      self.unit_cell = self.model.crystal_symmetry().unit_cell()
      miller_arrays = reflection_file_reader.any_reflection_file(file_name =
      self.map_coeffs_fn).as_miller_arrays()
      for m_array in miller_arrays:
        if m_array.info().label_string() == self.mtzlabel:
          self.ma = m_array

  def get_map_3d(self, volume_scale=False): 
    print('Doing FFT')
    fft_map = self.ma.fft_map(resolution_factor=0.10)
    if volume_scale:
      fft_map.apply_volume_scaling()
    else:
      fft_map.apply_sigma_scaling()
    self.map_3d = fft_map.real_map_unpadded()
    # stuff below done only to get a sense of vacuum level
    #map_1d = self.map_3d.as_1d()
    #import matplotlib.pyplot as plt
    #x=plt.hist(map_1d, bins=100)
    #print ("Max of peak height = ", x[1][np.argmax(x[0])], np.mean(map_1d))
    #plt.show()
    #from IPython import embed; embed(); exit()


  def get_map_value_at_point(self, x,y,z):
    #print ('getting map value at point')
    site_cart = (x,y,z)
    site_frac = self.unit_cell.fractionalize(site_cart)
    map_value = self.map_3d.eight_point_interpolation(site_frac)
    #map_value = self.map_3d.tricubic_interpolation(site_frac)
    return map_value

  def get_single_atom_coordinates(self, sele_str=None):
    print ('Getting single atom coordinates')
    if sele_str is None:
      sele_str = 'chain A and resname OEZ and name O5'
    isel = self.model.selection(sele_str).iselection()
    ph_sel = self.model.select(isel).get_hierarchy()
    atoms_sel = ph_sel.atoms()
    x,y,z=atoms_sel[0].xyz
    return (x,y,z)

  def get_Mn_Mn_distances_from_anom_map(self):
    self.read_in_model_mtz()
    self.get_map_3d()
    # Run this function for MN1, MN3, MN4, get peak position and then distances
    r_site = {}
    r_centroid = {}
    r_peakmax = {}
    mn1_selection_str = 'chain A and resid 601 and resname OEC and  name MN1'
    mn3_selection_str = 'chain A and resid 601 and resname OEC and  name MN3'
    mn4_selection_str = 'chain A and resid 601 and resname OEC and  name MN4'
    r_site['MN1'], r_centroid['MN1'], r_peakmax['MN1'] = self.get_anom_map_value_around_sphere(mn1_selection_str, 'MN1')
    r_site['MN3'], r_centroid['MN3'], r_peakmax['MN3'] = self.get_anom_map_value_around_sphere(mn3_selection_str, 'MN3')
    r_site['MN4'], r_centroid['MN4'], r_peakmax['MN4'] = self.get_anom_map_value_around_sphere(mn4_selection_str, 'MN4')
    print ('-------------------------------------------------------------------------------------------------------')
    # 1-4 distance
    atom1 = 'MN1'
    atom2 = 'MN4'
    dist_peakmax = (r_peakmax[atom1] - r_peakmax[atom2]).length()
    dist_centroid = (r_centroid[atom1] - r_centroid[atom2]).length()
    dist_model = (r_site[atom1] - r_site[atom2]).length()
    print ('%s-%s: PeakMax Distance=%.2f, Centroid Distance=%.2f, Model Distance=%.2f'%(atom1, atom2, dist_peakmax, dist_centroid, dist_model))
    # 1-3 
    atom1 = 'MN1'
    atom2 = 'MN3'
    dist_peakmax = (r_peakmax[atom1] - r_peakmax[atom2]).length()
    dist_centroid = (r_centroid[atom1] - r_centroid[atom2]).length()
    dist_model = (r_site[atom1] - r_site[atom2]).length()
    print ('%s-%s: PeakMax Distance=%.2f, Centroid Distance=%.2f, Model Distance=%.2f'%(atom1, atom2, dist_peakmax, dist_centroid, dist_model))
    # 3-4
    atom1 = 'MN3'
    atom2 = 'MN4'
    dist_peakmax = (r_peakmax[atom1] - r_peakmax[atom2]).length()
    dist_centroid = (r_centroid[atom1] - r_centroid[atom2]).length()
    dist_model = (r_site[atom1] - r_site[atom2]).length()
    print ('%s-%s: PeakMax Distance=%.2f, Centroid Distance=%.2f, Model Distance=%.2f'%(atom1, atom2, dist_peakmax, dist_centroid, dist_model))

    print ('-------------------------------------------------------------------------------------------------------')

  def f_map_calculator(self, newr):
    map_values = flex.double()
    coordinates = flex.vec3_double()
    local_map_value = self.get_map_value_at_point(newr[0], newr[1], newr[2])
    if local_map_value > self.threshold_map_value:
      map_values.append(local_map_value)
      coordinates.append(newr)
    return (map_values, coordinates)

  def get_anom_map_value_around_sphere(self,selection_str=None, atom_name=None):

    xs,ys,zs=self.get_single_atom_coordinates(sele_str=selection_str)
    map_value_at_selection_xyz = self.get_map_value_at_point(xs,ys,zs)
    print ('Map value at selection %s, coords= %s, map_value = %.2f'%(selection_str, str((xs,ys,zs)), map_value_at_selection_xyz))
    rs = col((xs, ys, zs))

    map_values = flex.double()
    coordinates = flex.vec3_double()
    #amplitudes=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    amplitudes=list(np.arange(0.01, 1.1, 0.1))
    # Scale number of points on a sphere with radius^2
    points_on_smallest_sphere = 10
    points_on_spheres = []
    for amplitude in amplitudes:
      # Doing area scaling instead of volume to save time
      points_on_spheres.append(int(points_on_smallest_sphere*(amplitude*amplitude/(amplitudes[0]*amplitudes[0]))))
      # Constant scaling, do not use as it gives misleading results for centroid
      #points_on_spheres.append(10000)
    self.threshold_map_value = 5.00 # Number needs to be tuned by user
    # Pick a random number within 0.2A of the central coordinate and sample the space
    for ii,amplitude in enumerate(amplitudes):
      for n_pt in range(points_on_spheres[ii]):
        dr = amplitude*col(flex.random_double_point_on_sphere())
        newr = dr+col((xs,ys,zs))
        local_map_value = self.get_map_value_at_point(newr[0], newr[1], newr[2])
        if local_map_value > self.threshold_map_value:
          map_values.append(local_map_value)
          coordinates.append(newr)
    
    x,y,z = coordinates.parts()
    # Centroid
    xc = flex.mean_weighted(x, map_values)
    yc = flex.mean_weighted(y, map_values)
    zc = flex.mean_weighted(z, map_values)
    rc = col((xc, yc, zc))
    # Now max peak position
    rm = col(coordinates[flex.max_index(map_values)])
    # Now print out the distances
    print ('%s: Distance between peak max and peak centroid = %.2f'%(atom_name, (rm-rc).length()))
    print ('%s: Distance between peak max and model position = %.2f'%(atom_name, (rm-rs).length()))
    print ('%s: Distance between peak centroid and model position = %.2f'%(atom_name, (rc-rs).length()))
    if False:
      import matplotlib.pyplot as plt
      plt.figure()
      rlength_s = (coordinates - rs).norms()
      rlength_m = (coordinates - rm).norms()
      plt.scatter(rlength_s, map_values,)
      plt.scatter(rlength_m, map_values,)

    # Now getting all the statistics around peakmax to get centroid position around peakmax
    # For debugging only
    if False:
      print ('-----------------STATS around peakmax, not model position -------------------------')
      map_values = flex.double()
      coordinates = flex.vec3_double()
      for amplitude in amplitudes:
        for n_pt in range(10000):
          dr = amplitude*col(flex.random_double_point_on_sphere())
          newr = dr+rm
          local_map_value = self.get_map_value_at_point(newr[0], newr[1], newr[2])
          if local_map_value > threshold_map_value:
            map_values.append(local_map_value)
            coordinates.append(newr)
          
      x,y,z = coordinates.parts()
      xc = flex.mean_weighted(x, map_values)
      yc = flex.mean_weighted(y, map_values)
      zc = flex.mean_weighted(z, map_values)
      rc = col((xc, yc, zc))
      # Now max peak position
      rm = col(coordinates[flex.max_index(map_values)])
      if False:
        plt.figure()
        rlength_s = (coordinates - rs).norms()
        rlength_m = (coordinates - rm).norms()
        plt.scatter(rlength_s, map_values,)
        plt.scatter(rlength_m, map_values,)

    return (rs, rc, rm)

if (__name__ == "__main__"):

  t0 = time.time()
  datasets = ['protein']
  timepoints = datasets 
  
  
  for ii, timepoint in enumerate(timepoints):
    print ('=================================================== Stats for dataset %s ======================================================'%(datasets[ii]))
    args = [timepoint, datasets[ii], 'anom',]
    mvc = map_value_calculator(args)
    mvc.get_Mn_Mn_distances_from_anom_map()
  print("Total time: %-8.4f"%(time.time()-t0))
