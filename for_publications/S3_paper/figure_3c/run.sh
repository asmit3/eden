mkdir sigma_scale_csv
mkdir absolute_scale_csv
# Get sigmal scale csv first
libtbx.python get_map_values.py 2F fobs 0 0 0
libtbx.python get_map_values.py 2F fobs 1 0 0
mv *.csv sigma_scale_csv
# Get absolute scale csv
libtbx.python get_map_values.py 2F fobs 0 1 1
libtbx.python get_map_values.py 2F fobs 1 1 1
mv *.csv absolute_scale_csv
# Run the scripts to generate fourier images. These should generate .pkl files
libtbx.python run_O5O6.py
libtbx.python run_Ca1O6.py
libtbx.python run_Mn1O6.py
# Now generate plots; say for Mn1-O6
libtbx.python plot_simulation_lines_for_sigmascale.py data_sm_Mn1O6.pkl MN1O6


