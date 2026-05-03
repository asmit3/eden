from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
f = 'split_0000.expt'


# In MMCIF land, beam direction is -Z and Y is -g


#Looking from beam direction ie source to sample direction

#  X<-------
#         |
#         |
#         Y
# Z is out of the paper


expt = ExperimentListFactory.from_json_file(f, check_format=False)[0]
# Get the U-matrix i.e crystal orientation
U=expt.crystal.get_U()
# Get the orthogonalization matrix B
B=expt.crystal.get_B()







from IPython import embed; embed(); exit()

