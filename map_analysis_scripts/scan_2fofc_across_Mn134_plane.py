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

# Goal is to scan 2Fo-Fc map over the Mn1-Mn3-Mn4 plane to check for peaks

class map_value_calculator(object):
  def __init__(self, args):
    self.state=args[0]
    # Sort out mode here whether fcalc or fobs
    if len(args)==1:
      self.mode = 'fobs' # Only 2 options, fobs or fcalc 
    if len(args) > 1:
      if args[1] == 'fobs':
        self.mode = 'fobs' # Only 2 options, fobs or fcalc 
      if args[1] == 'fcalc':
        self.mode = 'fcalc'
      if args[1] == 'phenixfmodel':
        self.mode = 'phenixfmodel'
      if args[1] == 'fofc':
        self.mode = 'fofc'
      if args[1] == 'fobs_with_fill':
        self.mode = 'fobs_with_fill'
      else:
        Sorry('2nd argument needs to be either fobs, fcalc or just leave it blank')
    # Now go !
    self.model_fn = args[0]+'.pdb'                    #   Name of model file
    if self.mode == 'fobs':
      self.map_coeffs_fn = args[0]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='2FOFCWT_no_fill,PH2FOFCWT_no_fill'
    self.stepsize=0.01
    if self.mode == 'fobs_with_fill':
      self.map_coeffs_fn = args[0]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='2FOFCWT,PH2FOFCWT'
    if self.mode == 'fcalc':
      self.map_coeffs_fn = args[0]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='F-model,PHIF-model'
    if self.mode == 'phenixfmodel':
      self.map_coeffs_fn = args[0]+'_fcalc.mtz'         # Might or might not have fcalc 
      self.mtzlabel='FMODEL,PHIFMODEL'
    if self.mode == 'fofc':
      self.map_coeffs_fn = args[0]+'.mtz'         # Might or might not have fcalc 
      self.mtzlabel='FOFCWT,PHFOFCWT'

  def read_in_model_mtz(self, fastmode=True):
    print ('Reading in model and mtz files')
    if True:
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
    fft_map = self.ma.fft_map(resolution_factor=0.25)
    if volume_scale:
      fft_map.apply_volume_scaling()
    else:
      fft_map.apply_sigma_scaling()
    self.map_3d = fft_map.real_map_unpadded()

  def get_map_value_at_point(self, x,y,z):
    site_cart = (x,y,z)
    site_frac = self.unit_cell.fractionalize(site_cart)
    map_value = self.map_3d.eight_point_interpolation(site_frac)
    return map_value

  def get_plane_grid_points(self,nhat, r1, r2, r3, offset=0.0, write_grid_points_to_pdb=False):
    # r1,r2,r3 define the 3 points defining the plane
    all_grid_points = flex.vec3_double()
    all_steps = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0] #,1.1,1.2, 1.3,1.4,1.5,1.6]
    all_steps = list(np.arange(0.00, 1.0, 0.01))
    all_steps = [float(ss) for ss in all_steps]
    nsteps = len(all_steps)
    self.nsteps=nsteps

    r1_plane = offset*nhat+r1 
    r2_plane = offset*nhat+r2 
    r3_plane = offset*nhat+r3 

    self.rotation_angles = np.arange(-30,40, 5.0) 
    self.zero_degree_index = list(self.rotation_angles).index(0.0)
    for ang in self.rotation_angles:
      r1new=(r1_plane-r3_plane).rotate_around_origin(nhat, ang, deg=True)+r3_plane
      for step_size in all_steps:
        grid_point = r3_plane+step_size*(r1new-r3_plane)
        all_grid_points.append(col(grid_point))
    if write_grid_points_to_pdb:
      template = 'HETATM53017  O   HOH S1460     -51.833 -45.864 294.872  1.00 30.00           O'
      with open('grid_points.pdb', 'w') as fout:
        for grid_point in all_grid_points:
          coord_line = template[0:30]+'%8.3f'%grid_point[0]+'%8.3f'%grid_point[1]+'%8.3f'%grid_point[2]+template[54:]+'\n'
          fout.write(coord_line)
    return all_grid_points

  def run(self, sele_str='chain A and resname OEC and (name MN1 or name MN3 or name MN4)'):
    self.read_in_model_mtz()
    self.get_map_3d()
    isel = self.model.selection(sele_str).iselection()
    ph_sel = self.model.select(isel).get_hierarchy()
    atoms_sel = ph_sel.atoms()
    r1,r2,r3=atoms_sel.extract_xyz()
    r1=col(r1)
    r2=col(r2)
    r3=col(r3)
    nhat=(r1-r2).cross(r1-r3)
    nhat=nhat.normalize()
    all_map_values_dict = {}
    map_value_floor = -33.0


    #for offset in [-0.5, 0.5, 0.0]:
    for offset in [0.0]:
      all_grid_points = self.get_plane_grid_points(nhat, r1,r2,r3, offset=offset)
      all_map_values_dict[offset] = []
      for grid_point in all_grid_points:
        x,y,z = grid_point
        map_value = self.get_map_value_at_point(x,y,z)
        if map_value < map_value_floor:
          map_value = 0.0
        all_map_values_dict[offset].append(map_value)

    all_map_values = flex.double(len(all_map_values_dict[0.0]))
    for offset in all_map_values_dict.keys():
          all_map_values += flex.double(all_map_values_dict[offset])
    all_map_values=all_map_values/len(all_map_values_dict.keys())


    show_plot=True
    if show_plot:
      # First transform all the 3d grid points onto 2d plane. Define x-axis as Mn4-Mn1
      xp = (r1-r3).normalize()
      r1orth = (r1-r3).rotate_around_origin(nhat, 90, deg=True)+r3
      yp = (r1orth-r3).normalize()
      all_x2d = []
      all_y2d = []
      Mn1_arc_x2d = []
      Mn1_arc_y2d = []
      rMn1_O = all_grid_points - r1
      dMn1_O = rMn1_O.norms()
      for jj, grid_point in enumerate(all_grid_points):
        x_2d = (col(grid_point)-r3).dot(xp) 
        y_2d = (col(grid_point)-r3).dot(yp)
        all_x2d.append(x_2d)
        all_y2d.append(y_2d)
        if dMn1_O[jj] < 2.3 and dMn1_O[jj] > 1.7 :
          Mn1_arc_x2d.append(x_2d)
          Mn1_arc_y2d.append(y_2d)
      map_min=min(all_map_values)
      map_max=max(all_map_values)
      colors = [(x-map_min)/(map_max-map_min) for x in all_map_values]
      #colors = [(x,x,x) for x in colors]
      #plt.scatter(all_x2d, all_y2d,c=colors)
      #plt.show()
        
      fig = plt.figure()
      ax = plt.axes(projection='3d')
      # Will be in sequence of self.nsteps points
      for ii in range(0, len(all_x2d)//self.nsteps):
        if ii==self.zero_degree_index:
          ax.plot3D(all_x2d[ii*self.nsteps:(ii+1)*self.nsteps], all_y2d[ii*self.nsteps:(ii+1)*self.nsteps], all_map_values[ii*self.nsteps:(ii+1)*self.nsteps],'--')
        else:
          ax.plot3D(all_x2d[ii*self.nsteps:(ii+1)*self.nsteps], all_y2d[ii*self.nsteps:(ii+1)*self.nsteps], all_map_values[ii*self.nsteps:(ii+1)*self.nsteps],)
      ax.scatter3D(0,0,0, color='purple', s=40)
      ax.scatter3D(Mn1_arc_x2d, Mn1_arc_y2d, [0]*len(Mn1_arc_x2d), c='slategrey',alpha=0.1)
      #ax.text(0,0,0, 'MN4')
      ax.scatter3D((r1-r3).length(),0,0, color='purple', s=40)
      #ax.text((r1-r3).length(),0,0, 'MN1')
      # For point r2
      x_2d = (r2-r3).dot(xp) 
      y_2d = (r2-r3).dot(yp)
      ax.scatter3D(x_2d, y_2d, 0.0, color='purple', s=50)
      #ax.text(x_2d, y_2d, 0, 'MN3')
      isel = self.model.selection('chain A and resname OEC and (name O5)').iselection()
      ph_sel = self.model.select(isel).get_hierarchy()
      atoms_sel = ph_sel.atoms()
      coords = atoms_sel.extract_xyz()
      if len(coords) == 1:
        rO5= coords[0]
        rO5 = col(rO5)
        x_2d = (rO5-r3).dot(xp)
        y_2d = (rO5-r3).dot(yp)
        ax.scatter3D(x_2d, y_2d, 0.0, color='red', s=50)
      ax.set_title('2D traces for state %s. Map = %s, plane offset=%.2f'%(self.state, self.mtzlabel, offset),)  
      #ax.set(xlabel='x_2d', ylabel='y_2d', zlabel='map value (sigma)')
      ax.view_init(elev=0., azim=-91)
      plt.legend()
      plt.tight_layout()
      #plt.savefig('figure_map_trace_%s.pdf'%(self.state), format='pdf')
      plt.show()

if (__name__ == "__main__"):

  t0 = time.time()
  mvc = map_value_calculator(sys.argv[1:])
  mvc.run()
