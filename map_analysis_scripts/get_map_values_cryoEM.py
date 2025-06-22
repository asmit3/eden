from __future__ import division, print_function
import sys, time, math
from iotbx import reflection_file_reader
from libtbx.utils import to_str
from iotbx.data_manager import DataManager
from scitbx.array_family import flex
from libtbx.utils import Sorry
from scitbx.matrix import col
from libtbx.program_template import ProgramTemplate
import numpy as np


class Program(ProgramTemplate):

  description = """
     phenix.python get_map_values_cryoEM.py emd_50019.map 9evx.cif atom_selection_str='chain A and resname OEX and (name O5 or name MN4)' filename_template='lm_Mn4O5' plot=True
"""

  datatypes = ['model', 'real_map', 'phil',]
  data_manager_options = ['model_skip_expand_with_mtrix']
  master_phil_str = """
  input_files {
    include scope iotbx.map_model_manager.map_model_phil_str
  }
  atom_selection_str = None
    .type = str
    .help = CCTBX Selection string of the 2 atoms along which to compute map values. Should strictly be 2 atoms
  filename_template= None
    .type = str
    .help = Filename template to write out map values. E.g. filename_template= 'lm_MN1MN4'
  plot = False
    .type = bool
    .help = Flag to plot the map values or not
  """
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    self.data_manager.has_real_maps(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    print ('Inputs OK', file=self.logger)

  def run(self):
    mvc_em = map_value_calculator_em(self.params, self.data_manager)
    mvc_em.get_map_value_along_bond_noNormalization_fixed_points()

class map_value_calculator_em(object):
  def __init__(self, params, data_manager):
    self.swap_order=True # Swap the order of atoms to plot MN1 on origin, not O6 in the case of MN1-O6
    self.model = data_manager.get_model() 
    if data_manager.get_model().crystal_symmetry() is None:
      data_manager.get_model().add_crystal_symmetry_if_necessary()
    self.crystal_symmetry = data_manager.get_real_map().crystal_symmetry()
    self.unit_cell = self.crystal_symmetry.unit_cell()
    self.map_3d = data_manager.get_real_map().map_data()
    self.sele_str = params.atom_selection_str
    self.filename_template = params.filename_template
    self.plot = params.plot

  def get_map_value_at_point(self, x,y,z):
    site_cart = (x,y,z)
    site_frac = self.unit_cell.fractionalize(site_cart)
    #map_value = self.map_3d.eight_point_interpolation(site_frac)
    map_value = self.map_3d.tricubic_interpolation(site_frac) # Seems to give smoother maps ?
    return map_value

  def get_map_value_along_bond_noNormalization_fixed_points(self,):
    isel = self.model.selection(self.sele_str).iselection()
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
    r_dist=(col((x1,y1,z1))-col((x2,y2,z2))).length()

    # Please define the points at which to interpolate map values
    gridpoints = np.arange(-4, 6, 0.02) 
    l = x2-x1
    m = y2-y1
    n = z2-z1
    # use parametric form x=x0+tl; y=y0+tm;z=z0+tn
    for dist in gridpoints:
      t = dist*math.sqrt((dist*dist)/(l*l+m*m+n*n))/abs(dist)
      x = x1+t*l
      y = y1+t*m
      z = z1+t*n
      map_values.append(self.get_map_value_at_point(x,y,z))
    with open('%s_map_values_noNormalization.csv'%(self.filename_template), 'w') as fout:
      fout.write('atoms, %s, %s\n'%(name1, name2))
      for dist, map_value in zip(gridpoints, map_values):
        fout.write("%.3f, %.2f\n"%(dist, map_value))
    if self.plot:
      import matplotlib.pyplot as plt
      plt.scatter(gridpoints, map_values)
      plt.xlabel('Distance, 0 point =%s, 2nd atom =%s'%(name1, name2))
      plt.ylabel('Map Value')
      plt.show()

if (__name__ == "__main__"):
  from iotbx.cli_parser import run_program

  t0 = time.time()
  run_program(program_class=Program)

  print("Total time: %-8.4f"%(time.time()-t0))
