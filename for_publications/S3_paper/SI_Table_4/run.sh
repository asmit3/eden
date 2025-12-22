# First you have to create a folder structure for each dataset 
# Copy over files to respective folders. Eg. for 0F (9.5keV) I have L10198_0F/0F.pdb and L10198_0F/0F.mtz (See how .mtz file was created below)
# Instructions on how to generate 0F.mtz file
# a. Copy over ANOM map coeffs from anomalous 0F dataset into a new .mtz file
# b. Copy over 2Fo-Fc map coeffs from conventional 0F dataset into the same .mtz file
# c. This is the .mtz data that is now used for the table (and also in figure 2a, 2b)
# d. The .pdb file is from refinement using conventional data



# Then just run the script after changing the following
# Before running the get_anom_map_centroid.py script, open the file and go to the bottom (__main__ function)
# You may change the following -  
# 1. Which experiment to run --> 9.5keV or 7keV. Uncomment as needed. You might have to change folder names based on your naming
# 2. Which intensities to use --> fobs (2Fo-Fc) or anom (dANO)
libtbx.python get_anom_map_centroid.py 
