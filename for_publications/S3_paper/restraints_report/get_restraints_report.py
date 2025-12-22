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
from tabulate import tabulate

def gauss(x, amplitude, mu, sigma):
  return amplitude*np.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))

class map_value_calculator(object):
  def __init__(self, mtz, pdb, geo):
    self.model_fn = pdb                    #   Name of model file
    self.map_coeffs_fn = mtz         # Might or might not have fcalc 
    self.geo_file = geo
    self.mtzlabel='FOFCWT,PHFOFCWT'
    self.datasets = ['distance pairs', 'Restraint set (harmonic)', 'FBHR']

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
    fft_map = self.ma.fft_map(resolution_factor=0.25,)
    if volume_scale:
      fft_map.apply_volume_scaling()
    else:
      fft_map.apply_sigma_scaling()
    self.map_3d = fft_map.real_map_unpadded()

  def get_map_value_at_point(self, x,y,z):
    #print ('getting map value at point')
    site_cart = (x,y,z)
    site_frac = self.unit_cell.fractionalize(site_cart)
    map_value = self.map_3d.tricubic_interpolation(site_frac)
    #map_value = self.map_3d.eight_point_interpolation(site_frac)
    return map_value

  def get_single_atom_coordinates(self, sele_str=None):
    print ('Getting single atom coordinates')
    sele_str='chain A and resname OEZ and name O5'
    isel = self.model.selection(sele_str).iselection()
    ph_sel = self.model.select(isel).get_hierarchy()
    atoms_sel = ph_sel.atoms()
    x,y,z=atoms_sel[0].xyz
    return (x,y,z)

  def get_central_point(self, sele_str = None):
    print ('Getting central point of OEC')
    sele_str = 'chain A and resid 601 and (name MN1 or name MN2 or name MN3 or name MN4)'
    isel = self.model.selection(sele_str).iselection()
    ph_sel = self.model.select(isel).get_hierarchy()
    atoms_sel = ph_sel.atoms()
    natoms = len(atoms_sel)
    xall = 0
    yall = 0
    zall = 0
    for atom_sel in atoms_sel:
      x,y,z=atom_sel.xyz
      xall += x
      yall += y
      zall += z
    return (xall/natoms, yall/natoms, zall/natoms)


  def run_single_atom_map_value(self):
    self.read_in_model_mtz()
    self.get_map_3d()
    x,y,z=self.get_single_atom_coordinates()
    map_value=self.get_map_value_at_point(x,y,z)
    print ('Map value at %.2f, %.2f, %.2f is %.2f'%(x,y,z,map_value))

  def boxgrid_about_point(self, central_point=None, box_length_angstrom = 5.0, grid_spacing_angstrom=0.20):
    gridpoints = flex.vec3_double()
    npoints = int(box_length_angstrom/grid_spacing_angstrom) 
    x0,y0,z0 = central_point
    xstart = x0-box_length_angstrom/2.0
    ystart = y0-box_length_angstrom/2.0
    zstart = z0-box_length_angstrom/2.0
    for ix in range(npoints):
      for iy in range(npoints):
        for iz in range(npoints):
          x = xstart + ix*grid_spacing_angstrom
          y = ystart + iy*grid_spacing_angstrom
          z = zstart + iz*grid_spacing_angstrom
          r = col((x,y,z))
          gridpoints.append(r)
    return gridpoints

  def get_map_flatness_about_point (self, central_point=None, box_length_angstrom=5.0):
    gridpoints = self.boxgrid_about_point(central_point=central_point, box_length_angstrom=box_length_angstrom)
    map_values = []
    with open('gridpoints.pdb', 'w') as fout:
      for gridpoint in gridpoints:
        x,y,z = gridpoint
        map_value = self.get_map_value_at_point(x,y,z) 
        map_values.append(map_value)
        fake_occu = (map_value+2.0)/25.0
        fake_adp = map_value + 2.0 # add an offset
        fout.write('HETATM53017  O   OOO Q1460     %6.3f %6.3f %6.3f   %1.2f %02.2f           O\n'%(x,y,z, fake_occu, fake_adp))
    # Get map flatness statistics here:
    neg_map_values = flex.double()
    pos_map_values = flex.double()
    map_threshold = 0.0 # threshold in sigma levels
    for map_value in map_values:
      if abs(map_value) > map_threshold:
        if map_value < 0:
          neg_map_values.append(map_value)
        if map_value >= 0:
          pos_map_values.append(map_value)
    neg_flatness = flex.mean(neg_map_values)
    pos_flatness = flex.mean(pos_map_values)
    from scipy.stats import skew
    skew_value = skew(map_values)
    print ('Average of Fo-Fc map values (signed) = %.2f, %.2f and skew = %.2f'%(neg_flatness, pos_flatness, skew_value))
    if False:
      import matplotlib.pyplot as plt
      plt.hist(map_values, bins=20, range=[-5, 5])
      plt.show()
    return skew_value, map_values

  def run_FBHR_rms_comparison(self, FBHR_pdb_file=None):
    # 
    print ('Running RMS comparison between HR and FBHR restraints')
    model_fns = [self.model_fn, FBHR_pdb_file]
    oec_sele_strs = ['chain A and resid 601'] #, 'chain a and resid 601']
    metal_oxy_list = ['MN1,O1', 'MN1,O3','MN2,O1', 'MN2,O2', 'MN2,O3', 'MN3,O2', 'MN3,O4', 'MN3,O5', 'MN4,O4', 'MN4,O5','MN4,W1', 'MN4,W2', 'CA1,O2', 'CA1,O1', 'CA1,O5', 'CA1,W3', 'CA1,W4']
    metal_metal_list = ['MN1,MN3', 'MN1,MN4', 'MN3,MN4','MN1,MN2', 'MN2,MN3', 'CA1,MN1', 'CA1,MN3', 'CA1,MN2']
    list_for_O6 = [] #'MN1,O6', 'CA1,O6', 'O5,O6']
    all_distances = []
    all_distances_dict = {}
    # Loop over all pdb files
    for model_fn in model_fns:
      pdb_inp = iotbx.pdb.input(file_name = model_fn)
      self.model = mmtbx.model.manager(model_input = pdb_inp)
      for oec_sele_str in oec_sele_strs:
        chainid = oec_sele_str[6]
        print ('Validating for selection: %s'%oec_sele_str)
        isel = self.model.selection(oec_sele_str).iselection()
        ph_sel = self.model.select(isel).get_hierarchy()
        atoms_sel = ph_sel.atoms()
        atom_names = atoms_sel.extract_name()
        atom_names = [atom_name.strip() for atom_name in atom_names]
        atom_coords = atoms_sel.extract_xyz()

        if '0F' in model_fn or '1F' in model_fn:
          final_pair_list = metal_metal_list+metal_oxy_list
        if '2F' in model_fn:
          final_pair_list = metal_metal_list+metal_oxy_list+list_for_O6

        for pair in final_pair_list:
          atom1, atom2 = pair.split(',')
          pair_nocomma = '%s:%s'%(atom1, atom2)
          distance = (col(atom_coords[atom_names.index(atom2)]) - col(atom_coords[atom_names.index(atom1)])).length()
          if '%s-%s'%(chainid, pair_nocomma) not in all_distances_dict.keys():
            all_distances_dict['%s-%s'%(chainid, pair_nocomma)] = []
          all_distances_dict['%s-%s'%(chainid, pair_nocomma)].append(distance)
        #all_distances.append(['%s-%s'%(chainid, pair),distance])
    #print (tabulate(all_distances))
    for pair in all_distances_dict.keys():
      str_to_write = ['%s'%pair]+all_distances_dict[pair]
      all_distances.append(str_to_write)

    delta_distances = flex.double([x[1]-x[2] for x in all_distances])
    print ('RMSD between dataset 1 and FBHR reference:%.3f'%(delta_distances.rms()))
    if False:
      with open('output.csv', 'w') as fout:
        for idx, entry in enumerate(all_distances):
          str_dist = [str('%.2f'%x) for x in all_distances[idx][1:]]
          fout.write('%s, '%entry[0]+', '.join(str_dist))
          fout.write('\n')
    print (tabulate(all_distances, headers=self.datasets, tablefmt='csv'))

  def get_residual_term_from_geo_file(self, has_slack=False, print_lines=False):
    bond_pdb_str = 'OEC A 601'  
    OEC_str = 'OEC'
    if '2F' in self.geo_file:
      bond_pdb_str = 'OEI A 601'  
      OEC_str = 'OEI'
    geo_residual_idx=5
    geo_delta_idx=2
    geo_sigma_idx=3 
    if has_slack:
      geo_residual_idx=6
      geo_delta_idx=3
      geo_sigma_idx=4 
    
    sum_residual = 0.0
    angle_sum_residual = 0.0
    zscores= flex.double()
    angle_zscores= flex.double()
    npairs = 0
    
    all_lines = []
    with open(self.geo_file, 'r') as fin:
      for line in fin:
        if line !='\n':
          all_lines.append(line.strip())
    for ii, line in enumerate(all_lines):
      if 'bond pdb=' in line:
        if bond_pdb_str in line:
          # Ensure only Ca1, MN4 and Mn3 are counted
          if OEC_str in all_lines[ii+1]: #'CA1' in all_lines[ii+1] or 'MN4' in all_lines[ii+1] or 'MN3' in all_lines[ii+1]: 
            if 'MN' in line and 'MN' in all_lines[ii+1]: continue
            residual = float(all_lines[ii+3].split()[geo_residual_idx])
            delta = float(all_lines[ii+3].split()[geo_delta_idx])
            sigma = float(all_lines[ii+3].split()[geo_sigma_idx])+1e-15
            sum_residual += abs(residual)
            zscores.append(delta/sigma)
            if print_lines:
              print (line)
              print (all_lines[ii+1])
        elif bond_pdb_str in all_lines[ii+1]:
          if OEC_str in line: #'CA1' in line or 'MN4' in line or 'MN3' in line: 
            if 'MN' in line and 'MN' in all_lines[ii+1]: continue
            residual = float(all_lines[ii+3].split()[geo_residual_idx])
            delta = float(all_lines[ii+3].split()[geo_delta_idx])
            sigma = float(all_lines[ii+3].split()[geo_sigma_idx])+1e-15
            sum_residual += abs(residual)
            zscores.append(delta/sigma)
            if print_lines:
              print (line)
              print (all_lines[ii+1])
      if ('angle pdb=' in line) and (not has_slack):
        if bond_pdb_str in line:
          if OEC_str in all_lines[ii+1] and OEC_str in all_lines[ii+2]:
            residual = float(all_lines[ii+4].split()[geo_residual_idx])
            delta = float(all_lines[ii+4].split()[geo_delta_idx])
            sigma = float(all_lines[ii+4].split()[geo_sigma_idx])+1e-15
            angle_sum_residual += abs(residual)
            angle_zscores.append(delta/sigma)
            if print_lines:
              print (line)
              print (all_lines[ii+1])
              print (all_lines[ii+2])
          
    abs_zscores = flex.abs(zscores)
    abs_angle_zscores = flex.abs(angle_zscores)
    print ('Total # of bonds= %d'%(len(zscores)))
    print ('Total # of angles= %d'%(len(angle_zscores)))
    #print ('Total residual term for %s is %.3f'%(bond_pdb_str, sum_residual))
    print ('# of bond pairs with |Z-scores| > 2.0=%d'%(sum(abs_zscores > 2.0)))
    print ('# of angle triples with |Z-scores| > 2.0=%d'%(sum(abs_angle_zscores > 2.0)))
    print ('Bond RMSZ = %.2f'%zscores.rms())
    if len(angle_zscores) > 0:
      print ('Angle RMSZ = %.2f'%angle_zscores.rms())
    return sum_residual
    # sqrt[sum((Z-Zmean)**2)/n]

  def run(self, has_slack=False, FBHR_pdb_file=None):
    box_length_angstrom = 5.0
    self.read_in_model_mtz()
    self.get_map_3d()
    print ('----------------------Getting Map Statistics from Bulk Solvent Region----------------------')
    central_point = (-136,-52.88,335.9) # bulk solvent region, OK for 0F,1F,2F
    skew_value, map_values = self.get_map_flatness_about_point(central_point=central_point, box_length_angstrom=box_length_angstrom)
    print ('----------------------Getting Map Statistics from OEC region----------------------')
    central_point = self.get_central_point()
    print ('Here is the central point of the box = %s'%(str(central_point)))
    skew_value, map_values = self.get_map_flatness_about_point(central_point=central_point, box_length_angstrom=box_length_angstrom)
    print ('----------------------Getting Statistics from Geo File for Restraints----------------------')
    residual = self.get_residual_term_from_geo_file(has_slack=has_slack)
    print ('----------------------Comparison of Harmonic Restraints (HR) with Flat Bottom Harmonic Restraints (FBHR)----------------------')
    if FBHR_pdb_file is not None:
      self.run_FBHR_rms_comparison(FBHR_pdb_file = FBHR_pdb_file)
    return (residual, skew_value, map_values)

if (__name__ == "__main__"):

  import sys
  t0 = time.time()
  # First file needs to be the HR file, the 2nd file can be the FBHR file
  timepoint, restraint_set = sys.argv[1], sys.argv[2]
  #timepoint = '0F'
  #restraint_set = 'HR'
  # 0F
  if timepoint == '0F':
    if restraint_set == 'HR':
      base_fnames = ['0F/L10136_0F_noanom_OEC_step7_optW_lessangles_1.97_118_15',] # step 4a with H
    if restraint_set == 'FBHR':
      base_fnames = ['0F/L10136_0F_noanom_OEC_waters_anom_addH_FBHR_1.97_118_30']
    if restraint_set == 'CCP4_default':
      base_fnames = ['0F/L10136_0F_noanom_OEC_waters_anom_addH_refmacrestr_1.97_118_30'] # sigma 0.02
    base_FBHR = ['0F/L10136_0F_noanom_OEC_waters_anom_addH_FBHR_1.97_118_30']

  # 1F
  if timepoint == '1F':
    if restraint_set == 'HR':
      base_fnames = ['1F/L10198_1F_noanom_OEC_step6_optW_lessangles_1.98_118_45']
    if restraint_set == 'FBHR':
      base_fnames = ['1F/L10198_1F_noanom_OEC_waters_addH_FBHR_1.98_118_30']
    if restraint_set == 'CCP4_default':
      base_fnames = ['1F/L10198_1F_noanom_OEC_waters_addH_refmacrestr_1.98_118_30']
    base_FBHR = ['1F/L10198_1F_noanom_OEC_waters_addH_FBHR_1.98_118_30']


  # 2F

  if timepoint == '2F':
    if restraint_set == 'HR':
      base_fnames = ['2F/L10198_2F_noanom_OEboth_1.93_step13_wCa1O6_0221_45']
    if restraint_set == 'FBHR':
      base_fnames = ['2F/L10198_2F_noanom_OEboth_1.93_addH_FBHR_0221_30'] 
    if restraint_set == 'CCP4_default':
      base_fnames = ['2F/L10198_2F_noanom_OEboth_1.93_addH_refmacrestr_0221_30']
    base_FBHR = ['2F/L10198_2F_noanom_OEboth_1.93_addH_FBHR_0221_30'] 

  all_residuals = []
  all_skews = []
  all_map_values = {}
  for base_fname in base_fnames:
    print ('#########################################################################################')
    has_slack=False
    pdb = base_fname+'.pdb'
    mtz = base_fname+'.mtz'
    geo = base_fname+'_final.geo'
    mvc = map_value_calculator(mtz, pdb, geo)
    if 'slack' in base_fname or 'FBHR' in base_fname: # and '1comp' in base_fname:
      has_slack=True
    if has_slack:
      residual, skew_value, map_values = mvc.run(has_slack=has_slack,)
    else:
      residual, skew_value, map_values = mvc.run(has_slack=has_slack, FBHR_pdb_file = base_FBHR[0]+'.pdb')
    all_residuals.append(residual)
    all_skews.append(skew_value)
    all_map_values[base_fname] = map_values

  if True:
    import matplotlib.pyplot as plt
    xvalues = np.arange(-4, 4, 0.1)
    for base_fname in base_fnames:
      plt.figure()
      if base_fname in ['2F_restr0.05_wO6', '2F_oldOEI']:
        data = plt.hist(all_map_values[base_fname], bins=60, range=[-4,4], histtype='bar', color='skyblue',) #, linewidth=2, label=base_fname,density=False)  
      else: # base_fname in ['2F_2comp_wO6_slack']:
        data = plt.hist(all_map_values[base_fname], bins=60, range=[-4,4], histtype='bar', color='wheat',) #, linewidth=2, label=base_fname,density=False)  
      value_at_zero = data[0][list(data[1]).index(0.0)]
      ygauss = gauss(xvalues, value_at_zero, 0.0, 1.0, )
      plt.plot(xvalues, ygauss, label=' $N(\mu=0, \sigma=1.0)$', color='k', alpha=0.6)
    #plt.xlabel('Fo-Fc Map values (sigma-scaled)')
    plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.ylim([0, 1250])
    plt.tight_layout()
    #plt.savefig('%s_hist.pdf'%base_fname)
    plt.show()
  if False:
    import matplotlib.pyplot as plt
    plt.scatter(all_residuals, all_skews)
    plt.xlabel('Residuals')
    plt.ylabel('Skew values of Fo-Fc map')
    plt.show()

  print("Total time: %-8.4f"%(time.time()-t0))
