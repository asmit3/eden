#! /bin/csh


setenv state $1
setenv iter $2

source ~/.cshrc


# Change work_folder and work_folder2
setenv work_folder /net/cci/asmit/Y2K_projects/MCR_end_rapid_Dec2021 
setenv work_folder2 /net/cci-filer3/home/asmit/Y2K_projects/MCR_end_rapid_Dec2021


setenv phil_shake ${work_folder}/phils/${state}_iters_shake/${state}_shake_${iter}.phil
setenv phil_continue ${work_folder}/phils/${state}_iters_continue/${state}_continue_${iter}.phil
setenv data ${work_folder}/data/${state}_F-obs-filtered_Fmodel_and_Rfree.mtz

cd $SCRATCH 

# Change this to where you will be running end/rapid on the CCI machines
cd asmit/end_rapid/MCR_end_rapid_Dec2021/kicked_data_v1


setenv refinedir ${state}_endrapidrefine_${iter}
rm -r ${refinedir}
mkdir $refinedir
cd $refinedir

source /net/cci/iris/source_ccp4_env.csh
echo "executing /net/cci/iris/projects/experiments/LCLS/LQ39/refinement/20180302_kicked_data_refinement/kick_data_bydiff.com $data seed=$iter"
# Please make sure you specify both F and Fc labels here explicitly !!!!
${work_folder2}/kick_data_bydiff.com $data seed=$iter F=F-obs-filtered FC=F-model
source /net/cci/xp/phenix/phenix-1.20-4459/build/setpaths.csh 
echo "executing /net/cci/iris/projects/experiments/LCLS/LQ39/refinement/20180411_kicked_data_and_model_refinement/stepwise_phenix_refine.sh $phil_shake 12 >> log.out"
${work_folder2}/stepwise_phenix_refine.sh $phil_shake 1 >> log.out
echo "executing /net/cci/iris/projects/experiments/LCLS/LQ39/refinement/20180411_kicked_data_and_model_refinement/stepwise_phenix_refine.sh $phil_continue 12 >> log.out"
${work_folder2}/stepwise_phenix_refine.sh $phil_continue 12 >> log.out
