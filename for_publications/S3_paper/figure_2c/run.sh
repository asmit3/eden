# Need to have both .pdb and .mtz (structure and map coeffs file from phenix.refine) files in the same folder with same prefix e.g. 0F.pdb and 0F.mtz
libtbx.python get_map_values.py 0F ANOM  # 7keV
libtbx.python get_map_values.py 2F ANOM  # 7keV
libtbx.python get_map_values.py 0F_L98 ANOM # 9.5keV 
libtbx.python get_map_values.py 2F_L98 ANOM # 9.5keV
# After writing out all map values, plot it as follows
libtbx.python plot_with_nodistnormalization.py MN1MN4 
