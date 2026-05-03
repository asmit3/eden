tmpdir='tmpiRA67k_tst_iota_ransac_193double_lattice'
json_glob=['/reg/neh/home/asmit/ExaFEL/dials/modules/xfel_regression/iota_test_data_and_scripts/LS49_sim/'+tmpdir+'/idx-step5_MPIbatch_000000.img_refined.expt']
json_glob=['/reg/neh/home/asmit/ExaFEL/dials/modules/xfel_regression/iota_test_data_and_scripts/LS49_sim/'+tmpdir+'/idx-composite_00000_refined.expt']+['/reg/neh/home/asmit/ExaFEL/dials/modules/xfel_regression/iota_test_data_and_scripts/LS49_sim/'+tmpdir+'/idx-composite_00000_refined.expt']
#json_glob='/reg/neh/home/asmit/ExaFEL/dials/modules/xfel_regression/iota_test_data_and_scripts/LS49_sim/01_stills_process/ransac_reference_answer/idx-step5_MPIbatch_000000.img_refined_experiments.json'
image_glob=['/reg/neh/home/asmit/ExaFEL/dials/modules/xfel_regression/iota_test_data_and_scripts/LS49_sim/01_stills_process/step5_MPIbatch_000000.img.gz']+['/reg/neh/home/asmit/ExaFEL/dials/modules/xfel_regression/iota_test_data_and_scripts/LS49_sim/01_stills_process/step5_MPIbatch_000001.img.gz']
lattice_dir='/reg/neh/home/asmit/ExaFEL/dials/modules/xfel_regression/iota_test_data_and_scripts/LS49_sim/01_stills_process'

from tst_iota_ransac import *

tag='idx-composite_'
compare_with_ground_truth(json_glob, image_glob, lattice_dir, tag) 
'''
tol_1=0.1
tol_2=0.1
from dxtbx.model.experiment_list import ExperimentListFactory
exp = ExperimentListFactory.from_json_file(json_glob, check_format=False)
exp0 = ExperimentListFactory.from_json_file(os.path.join(lattice_dir,'ransac_reference_answer', 'idx-step5_MPIbatch_000000.img_refined_experiments.json'), check_format=False)

crys = exp.crystals()
crys0 = exp0.crystals()
assert len(crys) == len(crys0) == 1, 'Number of indexed crystals should be 1 for single-lattice case'
from IPython import embed; embed(); exit()
assert Crystal(crys[0]).is_similar_to(crys0[0], angle_tolerance=tol_1, uc_rel_length_tolerance=tol_2), 'crystal model is dissimilar to reference model provided'
'''
