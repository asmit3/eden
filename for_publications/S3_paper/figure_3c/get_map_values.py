from __future__ import division, print_function
import sys, time, math
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
import copy


class map_value_calculator(object):
  def __init__(self, args):
    self.state=args[0]
    if int(args[2]) == 1:
      self.swap_order=True # Swap the order of atoms to plot MN1 on origin, not O6 in the case of MN1-O6
    if int(args[2]) == 0:
      self.swap_order=False # Swap the order of atoms to plot MN1 on origin, not O6 in the case of MN1-O6

    if int(args[3]) == 1:
      self.use_f000 = True
    if int(args[3]) == 0:
      self.use_f000 = False
    if int(args[4]) == 1:
      self.use_abs_scale = True
    if int(args[4]) == 0:
      self.use_abs_scale = False

    # Sort out mode here whether fcalc or fobs
    if len(args)==1:
      self.mode = 'fobs' # Only 2 options, fobs or fcalc 
    if len(args) > 1:
      if args[1] == 'fobs':
        self.mode = 'fobs' # Only 2 options, fobs or fcalc 
      if args[1] == 'fcalc':
        self.mode = 'fcalc'
      if args[1] == 'fofc':
        self.mode = 'fofc'
      if args[1] == 'ANOM':
        self.mode = 'ANOM'
      else:
        Sorry('2nd argument needs to be either fobs, fcalc or just leave it blank')
    # Now go !
    self.model_fn = args[0]+'.pdb'                    #   Name of model file
    if self.mode == 'fobs':
      self.map_coeffs_fn = args[0]+'.mtz'         # Might or might not have fcalc 
      #self.mtzlabel='FOFCWT,PHFOFCWT'
      #self.mtzlabel='2FOFCWT_no_fill,PH2FOFCWT_no_fill'
      self.mtzlabel='2FOFCWT,PH2FOFCWT'
    self.stepsize=0.01
    if self.mode == 'fcalc':
      self.map_coeffs_fn = args[0]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='F-model,PHIF-model'
    if self.mode == 'fofc':
      self.map_coeffs_fn = args[0]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='FOFCWT,PHFOFCWT'
    if self.mode == 'ANOM':
      self.map_coeffs_fn = args[0]+'.mtz'         # Might or might not have fcalc 
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

  def get_map_3d(self,): 
    print('Doing FFT')
    if self.use_abs_scale:
      volume_scale=True
    else:
      volume_scale=False

    if self.use_f000:
      fft_map = self.ma.fft_map(resolution_factor=0.10, f_000=3056554.731)
    if not self.use_f000:
      fft_map = self.ma.fft_map(resolution_factor=0.10,)
    #fft_map_copy = copy.deepcopy(fft_map)

    if volume_scale:
      fft_map.apply_volume_scaling()
    else:
      fft_map.apply_sigma_scaling()
    self.map_3d = fft_map.real_map_unpadded()
    #map_3d_sigma = fft_map_copy.apply_sigma_scaling().real_map_unpadded()
    #from IPython import embed; embed(); exit()
    # stuff below done only to get a sense of vacuum level
    #map_1d = self.map_3d.as_1d()
    #import matplotlib.pyplot as plt
    #x=plt.hist(map_1d, bins=100)
    #print ("Max of peak height = ", x[1][np.argmax(x[0])], np.mean(map_1d))
    #plt.show()

  def get_map_value_at_point(self, x,y,z):
    #print ('getting map value at point')
    site_cart = (x,y,z)
    site_frac = self.unit_cell.fractionalize(site_cart)
    map_value = self.map_3d.tricubic_interpolation(site_frac)
    #map_value = self.map_3d.eight_point_interpolation(site_frac)
    return map_value

  def get_single_atom_coordinates(self, sele_str=None):
    print ('Getting single atom coordinates')
    sele_str='chain A and resname OEZ and name O6'
    isel = self.model.selection(sele_str).iselection()
    ph_sel = self.model.select(isel).get_hierarchy()
    atoms_sel = ph_sel.atoms()
    x,y,z=atoms_sel[0].xyz
    return (x,y,z)

  def run_single_atom_map_value(self):
    self.read_in_model_mtz()
    self.get_map_3d()
    x,y,z=self.get_single_atom_coordinates()
    map_value=self.get_map_value_at_point(x,y,z)
    print ('Map value at %.2f, %.2f, %.2f is %.2f'%(x,y,z,map_value))

  def get_map_value_along_bond(self, sele_str=None, filename_template=None):
    self.read_in_model_mtz()
    self.get_map_3d()

    #sele_str='chain A and resname OEZ and (name MN1 or name MN4)'
    isel = self.model.selection(sele_str).iselection()
    ph_sel = self.model.select(isel).get_hierarchy()
    atoms_sel = ph_sel.atoms()
    map_values = flex.double()
    x1,y1,z1=atoms_sel[0].xyz
    x2,y2,z2=atoms_sel[1].xyz
    name1=atoms_sel[0].name.strip()
    name2=atoms_sel[1].name.strip()
    #print (x1,x2)
    r_dist=(col((x1,y1,z1))-col((x2,y2,z2))).length()
    l = x2-x1
    m = y2-y1
    n = z2-z1
    stepsize=self.stepsize
    steps = np.arange(-1, 2, stepsize)
    for step in steps:
      x = l*step+x1
      y = y1+m*(x-x1)/l
      z = z1+n*(x-x1)/l
      map_values.append(self.get_map_value_at_point(x,y,z))
    # Filename template would be like lm_O5O6 or sm_O5O6 depending on the monomer
    with open('%s_%s_map_values.csv'%(self.state, filename_template), 'w') as fout:
      fout.write('atoms, %s, %s\n'%(name1, name2))
      for step, map_value in zip(steps, map_values):
        fout.write("%.3f, %.2f\n"%(step, map_value))
    return map_values

  def get_map_value_along_bond_noNormalization_fixed_points(self, sele_str=None, filename_template=None):
    """ Do not normalize distance along bond. Instead write out map value at certai distances from center atom"""

    #sele_str='chain A and resname OEZ and (name MN1 or name MN4)'
    isel = self.model.selection(sele_str).iselection()
    ph_sel = self.model.select(isel).get_hierarchy()
    atoms_sel = ph_sel.atoms()
    map_values = flex.double()
    if self.swap_order:
      x1,y1,z1=atoms_sel[1].xyz
      x2,y2,z2=atoms_sel[0].xyz
      name1=atoms_sel[1].name.strip()
      name2=atoms_sel[0].name.strip()
    else:
      x1,y1,z1=atoms_sel[0].xyz
      x2,y2,z2=atoms_sel[1].xyz
      name1=atoms_sel[0].name.strip()
      name2=atoms_sel[1].name.strip()
    #print (x1,x2)
    r_dist=(col((x1,y1,z1))-col((x2,y2,z2))).length()

    # Please define the points at which to interpolate map values
    gridpoints = np.arange(-4, 4, 0.02) 

    l = x2-x1
    m = y2-y1
    n = z2-z1
    stepsize=self.stepsize
    # use parametric form x=x0+tl; y=y0+tm;z=z0+tn
    with open('%s_%s_points.pdb'%(self.state, filename_template,), 'w') as fout:
      for dist in gridpoints:
        t = dist*math.sqrt((dist*dist)/(l*l+m*m+n*n))/abs(dist)
        x = x1+t*l
        y = y1+t*m
        z = z1+t*n
        map_value = self.get_map_value_at_point(x,y,z) 
        map_values.append(map_value)
        fake_occu = (map_value+2.0)/25.0
        fake_adp = map_value + 2.0 # add an offset
        fout.write('HETATM53017  O   OOO Q1460     %6.3f %6.3f %6.3f   %1.2f %02.2f           O\n'%(x,y,z, fake_occu, fake_adp))
    # Filename template would be like lm_O5O6 or sm_O5O6 depending on the monomer
    with open('%s_%s_map_values_noNormalization.csv'%(self.state, filename_template), 'w') as fout:
      fout.write('atoms, %s, %s\n'%(name1, name2))
      for dist, map_value in zip(gridpoints, map_values):
        fout.write("%.3f, %.2f\n"%(dist, map_value))

  def get_map_value_summed_over_sphere(self, sele_str=None, filename_template=None):
    from scitbx import regular_grid_on_unit_sphere
    self.read_in_model_mtz()
    self.get_map_3d()

    #sele_str='chain A and resname OEZ and (name MN1 or name MN4)'
    isel = self.model.selection(sele_str).iselection()
    ph_sel = self.model.select(isel).get_hierarchy()
    atoms_sel = ph_sel.atoms()
    map_values = flex.double()
    x1,y1,z1=atoms_sel[0].xyz
    x2,y2,z2=atoms_sel[1].xyz
    name1=atoms_sel[0].name.strip()
    name2=atoms_sel[1].name.strip()
    #print (x1,x2)
    l = x2-x1
    m = y2-y1
    n = z2-z1
    stepsize=self.stepsize
    steps = np.arange(-1, 2, stepsize)
    for step in steps:
      xs = l*step+x1
      ys = y1+m*(xs-x1)/l
      zs = z1+n*(xs-x1)/l
      spherical_grid=regular_grid_on_unit_sphere.rosca(m=5,hemisphere=False)
      local_map_value = 0.0
      for xn,yn,zn in spherical_grid: 
        x=xn+xs
        y=yn+ys
        z=zn+zs
        local_map_value += self.get_map_value_at_point(x,y,z) 
      map_values.append(local_map_value)
    # Filename template would be like lm_O5O6 or sm_O5O6 depending on the monomer
    with open('%s_%s_map_values.csv'%(self.state, filename_template), 'w') as fout:
      fout.write('atoms, %s, %s\n'%(name1, name2))
      for step, map_value in zip(steps, map_values):
        fout.write("%.3f, %.2f\n"%(step, map_value))
    return map_values

  def get_map_value_around_sphere(self, sele_str=None, filename_template=None):
    """ This function might be incomplete"""
    self.read_in_model_mtz()
    self.get_map_3d()

    #sele_str='chain A and resname OEZ and (name MN1 or name MN4)'
    xs,ys,zs=self.get_single_atom_coordinates()
    spherical_grid=regular_grid_on_unit_sphere.rosca(m=5,hemisphere=False)
    map_values = flex.double()
    for xn,yn,zn in spherical_grid:
      x=xn+xs
      y=yn+ys
      z=zn+zs
      local_map_value = self.get_map_value_at_point(x,y,z)
      map_values.append(local_map_value)
    return map_values

if (__name__ == "__main__"):

  t0 = time.time()

  # Loop over multiple selection strings here

  chainA='A'
  monomer='lm'
  if int(sys.argv[3]) == 1:
    swap_order = True
  if int(sys.argv[3]) == 0:
    swap_order = False
  if swap_order:
    sele_strings = [
                    'chain %s and resname OEZ and (name MN1 or name O1)'%(chainA),
                    'chain %s and resname OEZ and (name MN1 or name O3)'%(chainA),
                    'chain %s and resname OEZ and (name MN1 or name O6)'%(chainA),
                    'chain %s and resname OEZ and (name MN2 or name O1)'%(chainA),
                    'chain %s and resname OEZ and (name MN2 or name O2)'%(chainA),
                    'chain %s and resname OEZ and (name MN2 or name O3)'%(chainA),
                    'chain %s and resname OEZ and (name MN3 or name O2)'%(chainA),
                    'chain %s and resname OEZ and (name MN3 or name O3)'%(chainA),
                    'chain %s and resname OEZ and (name MN3 or name O4)'%(chainA),
                    'chain %s and resname OEZ and (name MN3 or name O5)'%(chainA),
                    'chain %s and resname OEZ and (name MN4 or name O4)'%(chainA),
                    'chain %s and resname OEZ and (name MN4 or name O5)'%(chainA),
                    'chain %s and resname OEZ and (name MN4 or name W1)'%(chainA),
                    'chain %s and resname OEZ and (name MN4 or name W2)'%(chainA),
                    'chain %s and resname OEZ and (name CA1 or name O1)'%(chainA),
                    'chain %s and resname OEZ and (name CA1 or name O2)'%(chainA),
                    'chain %s and resname OEZ and (name CA1 or name O5)'%(chainA),
                    'chain %s and resname OEZ and (name CA1 or name W3)'%(chainA),
                    'chain %s and resname OEZ and (name CA1 or name W4)'%(chainA),
                    'chain %s and resname OEZ and (name CA1 or name O6)'%(chainA),
  ]
    file_name_templates = [
                    '%s_MN1O1'%(monomer),
                    '%s_MN1O3'%(monomer),
                    '%s_MN1O6'%(monomer),
                    '%s_MN2O1'%(monomer),
                    '%s_MN2O2'%(monomer),
                    '%s_MN2O3'%(monomer),
                    '%s_MN3O2'%(monomer),
                    '%s_MN3O3'%(monomer),
                    '%s_MN3O4'%(monomer),
                    '%s_MN3O5'%(monomer),
                    '%s_MN4O4'%(monomer),
                    '%s_MN4O5'%(monomer),
                    '%s_MN4W1'%(monomer),
                    '%s_MN4W2'%(monomer),
                    '%s_CA1O1'%(monomer),
                    '%s_CA1O2'%(monomer),
                    '%s_CA1O5'%(monomer),
                    '%s_CA1W3'%(monomer),
                    '%s_CA1W4'%(monomer),
                    '%s_CA1O6'%(monomer),
  ]
  #
  else:
    sele_strings = [ 
                    'chain %s and resname OEZ and (name O5 or name O6)'%(chainA),
  ]
    file_name_templates = [
                    '%s_O5O6'%(monomer),
  ]
  #
  # Command line inputs are 2F, fobs, swap_order_bool, use_f000_bool, use_abs_scale_bool
  
  mvc = map_value_calculator(sys.argv[1:])
  mvc.read_in_model_mtz()
  mvc.get_map_3d()

  if True:
    for sele_str, filename_template in zip(sele_strings, file_name_templates):
      print ('Working on filename template %s'%(filename_template))
      mvc.get_map_value_along_bond_noNormalization_fixed_points(sele_str=sele_str, filename_template=filename_template)

  print("Total time: %-8.4f"%(time.time()-t0))
