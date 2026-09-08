import os
from iotbx.data_manager import DataManager
import mmtbx.model
import iotbx.pdb
from scitbx.array_family import flex
from iotbx import reflection_file_reader

class structural_similarity(object):
  def __init__(self, reference_pdb,  endrapid_prefix_path, sele_str_to_compare):
    self.reference_pdb = reference_pdb
    self.endrapid_prefix_path = endrapid_prefix_path
    self.sele_str = sele_str_to_compare
    self.mtzlabel='mFo-DFc_omit,PHImFo-DFc_omit'


  def get_reference_coordinates(self):
    f = os.path.join(self.reference_pdb)
    dm = DataManager()
    dm.set_overwrite(True)
    m=dm.get_model(filename=f)
    self.unit_cell = m.crystal_symmetry().unit_cell()
    m2=m.select(m.selection(self.sele_str))
    pdb_hierarchy = m2.get_hierarchy()
    names = pdb_hierarchy.atoms().extract_name()
    self.names = [x.strip() for x in names]
    self.reference_xyz = pdb_hierarchy.atoms().extract_xyz()
    
  def compare_with_endrapid(self):
    all_similarities = []
    all_omit_map_values = []
    nmax = 56
    for i in range(1, nmax, 1):
      timepoint='2F'
      f = os.path.join('%s_%d/2F_endrapid_0624_12.pdb'%(self.endrapid_prefix_path, i,))
      f_omit_O6 = os.path.join('%s_%d/omit_O6/2F_endrapid_0624_12_polder_map_coeffs.mtz'%(self.endrapid_prefix_path, i,))

      pdb_inp = iotbx.pdb.input(file_name = f)
      m = mmtbx.model.manager(model_input = pdb_inp)


      m2 = m.select(m.selection(self.sele_str))
      xyz = m2.get_hierarchy().atoms().extract_xyz()
      similarity = self.compare_coordinates(xyz)
      all_similarities.append(similarity)

      all_omit_map_values.append(self.get_map_value_at_oxygen_site(f_omit_O6, xyz[5], ))
    from IPython import embed; embed(); exit()

  def compare_coordinates(self, xyz):
    return (xyz - self.reference_xyz).norm()

  def get_map_value_at_oxygen_site(self, map_coeffs_fn, site_cart):
    miller_arrays = reflection_file_reader.any_reflection_file(file_name = map_coeffs_fn).as_miller_arrays()
    for m_array in miller_arrays:
      if m_array.info().label_string() == self.mtzlabel:
        ma = m_array
    fft_map = ma.fft_map(resolution_factor=0.25,)
    fft_map.apply_sigma_scaling()
    map_3d = fft_map.real_map_unpadded()
    site_frac = self.unit_cell.fractionalize(site_cart)
    map_value = map_3d.tricubic_interpolation(site_frac)
    return map_value



if __name__ == "__main__":

  pdb = '/pscratch/sd/a/asmit/end_rapid_highres_S3/pdb/L10198_2F_noanom_OEboth_1.93_MnOrestr_noO5O6_0221_45.pdb' 
  endrapid_prefix_path = '/pscratch/sd/a/asmit/end_rapid_highres_S3/calculations/2F_endrapidrefine' 
  sele_str = 'chain A and resname OEI'
  oec_ss = structural_similarity(pdb, endrapid_prefix_path, sele_str)
  oec_ss.get_reference_coordinates()
  oec_ss.compare_with_endrapid()
