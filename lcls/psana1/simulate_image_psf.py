from argparse import ArgumentParser
from iotbx import file_reader
from cctbx.eltbx import henke
import mmtbx.command_line.fmodel
import mmtbx.utils


parser = ArgumentParser()
parser.add_argument("--psf", action='store_true')
parser.add_argument("--rescale", action="store_true")
args = parser.parse_args()

def get_fcalc_from_pdb(
        pdb_file,
        wavelength=None,
        dmin=1,
        dmax=None,
        k_sol=0.435, b_sol=46):
    """
    produce a structure factor from PDB coords, see mmtbx/command_line/fmodel.py for formulation
    k_sol, b_sol form the solvent component of the Fcalc: Fprotein + k_sol*exp(-b_sol*s^2/4) (I think)
    """

    pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
    pdb_in.assert_file_type("pdb")
    xray_structure = pdb_in.file_object.xray_structure_simple()
    xray_structure.show_summary()
    for sc in xray_structure.scatterers():
        if wavelength is not None:
            expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
            sc.fp = expected_henke.fp()
            sc.fdp = expected_henke.fdp()
    phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
    params2 = phil2.extract()
    params2.high_resolution = dmin
    params2.low_resolution = dmax
    params2.fmodel.k_sol = k_sol
    params2.fmodel.b_sol = b_sol
    params2.structure_factors_accuracy.algorithm = 'fft'
    f_model = mmtbx.utils.fmodel_from_xray_structure(
        xray_structure=xray_structure,
        f_obs=None,
        add_sigmas=False,
        params=params2).f_model
    f_model = f_model.generate_bijvoet_mates()

    return f_model.amplitudes().set_observation_type_xray_amplitude()


from dxtbx.model.crystal import Crystal
from copy import deepcopy

from dxtbx.model import Panel
from cctbx import uctbx
from scitbx.matrix import sqr, rec, col
import numpy as np
from scipy.spatial.transform import Rotation
from IPython import embed

from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal

#from simtbx.diffBragg.sim_data import SimData
from simtbx.nanoBragg.sim_data import SimData
from simtbx.diffBragg import utils
from simtbx.diffBragg.utils import fcalc_from_pdb
from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
import time

t1=time.time()
ucell = (106.9, 106.9, 303.9, 90, 90, 90)
symbol = "P41212"
# generate a random raotation
rotation = Rotation.random(num=1, random_state=100)[0]
Q = rec(rotation.as_quat(), n=(4, 1))
rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

# generate a small perturbation rotation
np.random.seed(1)
perturb_rot_axis = np.random.random(3)
perturb_rot_axis /= np.linalg.norm(perturb_rot_axis)
perturb_rot_ang = 0.15  # degree random perturbtation

# make the ground truth crystal:
a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
C = Crystal(a_real, b_real, c_real, symbol)
C.rotate_around_origin(rot_axis, rot_ang)

# Setup the simulation and create a realistic image
# with background and noise
# <><><><><><><><><><><><><><><><><><><><><><><><><>
nbcryst = NBcrystal()
nbcryst.dxtbx_crystal = C   # simulate ground truth
nbcryst.thick_mm = 0.1
Ncells_gt = 12, 12, 12
nbcryst.Ncells_abc = Ncells_gt  # ground truth Ncells
#nbcryst.miller_array = fcalc_from_pdb(wavelength=1.7,dmin=2.0, dmax=30.0)
nbcryst.miller_array = get_fcalc_from_pdb('MMO.pdb', wavelength=1.73,dmin=2.0, dmax=30.0)
print("Ground truth ncells = %f" % (nbcryst.Ncells_abc[0]))

SIM = SimData()
SIM.detector = SimData.simple_detector(150, 0.133, (2560, 2560))

# TODO get the detector model
node = SIM.detector[0]
node_d = node.to_dict()
Origin = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
distance = Origin[2]
print ("Ground truth originZ=%f" % (SIM.detector[0].get_origin()[2]))

def normal(x,mu,sigma):
    return ( 2.*np.pi*sigma**2. )**-.5 * np.exp( -.5 * (x-mu)**2. / sigma**2. )

from cctbx import factor_ev_angstrom
waves=np.linspace(7080,7160, 10)
flux=1e12*normal(waves, 7122, 12) # 12 eV sigma, mean is 7122eV
#import matplotlib.pyplot as plt
#plt.scatter(waves, flux)
#plt.show()
waves=factor_ev_angstrom/waves
spectrum_arr = [(x,y) for x,y in zip(waves, flux)]
SIM.beam.spectrum=[(1.73, 1e12),]
#SIM.beam.spectrum=spectrum_arr

#from IPython import embed; embed(); exit()

SIM.crystal = nbcryst
SIM.instantiate_diffBragg(oversample=0, auto_set_spotscale=True)
SIM.D.nopolar = False
SIM.D.default_F = 0
SIM.D.progress_meter = False
SIM.water_path_mm = 0.005
SIM.air_path_mm = 0.1
SIM.add_air = False
SIM.add_Water = False
SIM.include_noise = False
SIM.D.verbose=2

SIM.D.add_diffBragg_spots()
spots = SIM.D.raw_pixels.as_numpy_array()
SIM.D.readout_noise_adu = 3
SIM._add_background()
SIM._add_noise()


from libtbx.easy_pickle import dump
img_pre_psf = SIM.D.raw_pixels.as_numpy_array()
dump('no_psf.pickle', img_pre_psf)

if True:
    img_pre_psf = SIM.D.raw_pixels.as_numpy_array()
    v = SIM.D.verbose
    SIM.D.verbose = 8
    SIM.D.detector_psf_kernel_radius_pixels = 0
    fwhm = 76./133
    radius = 3
    SIM.D.apply_psf(shapetype.Fiber, fwhm, radius)
    SIM.D.verbose = v

print("Using oversample %d" % SIM.D.oversample)

# This is the ground truth image:
img = SIM.D.raw_pixels.as_numpy_array()
dump('with_psf.pickle', img)
t2=time.time()
print ('Total time taken in seconds = %d'%(t2-t1))
import matplotlib.pyplot as plt
plt.imshow(img,vmax=200)
plt.show()
#from IPython import embed; embed(); exit()
