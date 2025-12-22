# I will use 0F as an example
# Create a 0F folder
# 
# HR restraints: Copy over the .pdb, .mtz (map coeffs) and _final.geo files from a phenix.refine run to 0F folder --> named as L10136_0F_noanom_OEC_step7_optW_lessangles_1.97_118_15 in the python script provided. Change as needed
# FBHR restraints: Copy over .pdb, .mtz, _final.geo files to 0F folder --> named as L10136_0F_noanom_OEC_waters_anom_addH_FBHR_1.97_118_30. Change as needed
# CCP4 monlib restraints: Copy over .pdb, .mtz, _final.geo files to 0F folder --> named as L10136_0F_noanom_OEC_waters_anom_addH_refmacrestr_1.97_118_30. Change as needed

# Now run the following script for 0F
libtbx.python get_restraints_report.py 0F HR
libtbx.python get_restraints_report.py 0F FBHR
libtbx.python get_restraints_report.py 0F CCP4_default
