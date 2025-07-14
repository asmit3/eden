#! /bin/csh -q
#$ -cwd
#$ -o kickboth.out -j y -N MCR_XFEL


# Change MCR_XFEL to whatever ${tag} you have
setenv flashstate MCR_XFEL 

# Change the output path to wherever your output folder is  
./kick_and_refine.csh $flashstate $SGE_TASK_ID >& /net/cci/asmit/Y2K_projects/MCR_end_rapid_Dec2021/output/${flashstate}/kick.$SGE_TASK_ID.out

exit
