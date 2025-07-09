END/RAPID tutorial instructions:
Run the end/rapid calculations on an example case of MCR. 
Contact: Asmit Bhowmick [ email: abhowmick@lbl.gov ]

We will be doing all this work on the CCI machine Lewis because we need to use the cluster to run end/rapid. Log in 
<username>@lewis.lbl.gov where <username> is the same as used to log into dials or viper.

I am assuming you will be working in the home folder (i.e cd ~). Create your working folder there (I will call it "MCR_end_rapid_tutorial")

Now we begin:

1. First copy over the folder "tutorial_data" to your working folder. I am calling my working folder "MCR_end_rapid_tutorial".This tutorial_data folder contains all the data files required to run end/rapid on the MCR example here.
  A brief sampling of the files here -  
  (a) .eff, .mtz and .pdb files from the final phenix.refine run. 
  (b) The fasta sequence file for the protein

2. Now in the "MCR_end_rapid_tutorial" folder, create the following subfolders - (a) cifs (b) data (c) output (d) pdb (e) phils (f) results

3. Copy over the cif restraints files (f43.cif, com.cif, dya_ligands.cif, PEG.cif, tp7.cif) to the "cifs" folder. Once done copying, you should change the e.s.d value of the restraints to a higher number. In the tutorial data you should change the values of e.s.d values in the f43.cif, com.cif and tp7.cif files since these are related to our region of interest. In these files, you can increase all the bond e.s.d values to 0.200 and all the angle e.s.d values to 20.0 degrees. For your convenience, I have also provided the final modified files in case you want to consult it in the folder "final_cifs_to_be_used".
 By increasing the e.s.d values, we allow end/rapid to explore more space during refinement which can give us a better estimate of the error in our models. There is no definite answer on how much these e.s.d values should be increased. I have found it to be a tradeoff of resolution and how complex the region of interest you are studing is. It is instructive sometimes to run a couple of quick refinements with these looser restraints to make sure the system does not fall apart. Also it helps to talk to the person who did the original refinement to get a sense of what these loosened restraints should be. 

4. Now copy over the .mtz file (MCR_REV_refine.mtz) to the "data" folder. The file needs to be renamed in a very specific format namely ${tag}_F-obs-filtered_Fmodel_and_Rfree.mtz where ${tag} is the tag we will be using for the protein. For this tutorial I will use ${tag} as MCR_XFEL. Hence the .mtz file will be named as MCR_XFEL_F-obs-filtered_Fmodel_and_Rfree.mtz when copied over to the "data" folder.

5. Copy over the final refined .pdb file (MCR_REV_refine.pdb) to the "pdb" folder

6. I am assuming you have phenix sourced when doing this work. Make note that the version of phenix you are using is the EXACT version that was used to do the final refinement. For the tutorial data, this is phenix v.1.20.4459 and this can be sourced on the CCI computer as "source /net/cci/xp/phenix/phenix-1.20-4459/phenix_env.sh". 
  Once you have sourced phenix, run the following in the tutorial_data folder - "phenix.refine --diff_params "MCR_REV_refine_5.eff" >> MCR_XFEL.phil". This step writes out a smaller .phil file which only contains the non-default parameters that were used in the refinement. Also take note that I have named the phil file ${tag}.phil and this consistency in naming should be retained.

7. Once you have done step 6, copy over the MCR_XFEL.phil file over to the "phil" folder. Also copy over script "create_all_endrapid_phils.sh" from the "scripts" folder inside "tutorial_data"

8. Change directory to the "phils" folder. Create a copy of the MCR_XFEL.phil file and call it "MCR_XFEL_shake_0.phil". Again note the naming is ${tag}_shake_0.phil and this naming convention is required.
  Now open the MCR_XFEL_shake_0.phil file and we have to make a few changes here. The sections below in triple quotes show which parts of the file ned to be changed and describe the details. The EXACT lines that need to be changed are marked by X
  
  IMPORTANT ** make sure there is no space on either side of the = sign in the changes you are going to make below.
  For example if you have 
  file_name = my.mtz
  change that to 
  file_name=my.mtz
  ** 


 (a) Change the file_name line here (1X) to the full path of where your refined pdb file is. This should be in the "pdb" folder where you copied the .pdb file in step 5. 

  """
  input {
    pdb {
1X    file_name = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/Refine_4/MCR_REV_refine_4-coot-0.pdb"
    }
  """

  (b) Change the xray data file_name (1X) to just "kicked.mtz" (no paths). Change the labels (2X) to "F-obs-filtered,SIGF-obs-filtered". Change the r_free_flags file_name (3X) to the full path of the "MCR_XFEL_F-obs-filtered_Fmodel_and_Rfree.mtz" created in step 4.  

  """
  xray_data {
1X  file_name = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/Refine_4/MCR_REV_refine_4.mtz"
2X  labels = "I-obs,SIGI-obs"
    r_free_flags {
3X    file_name = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/Chris_1.19/Runs44_53_refine_data (3).mtz"
      label = "R-free-flags"
      test_flag_value = 1
  """

  (c) Change all the file_names below to the paths of the .cif files in the "cifs" folder (as done in step 3) 
 
  """
    monomers {
X     file_name = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/dya.ligands.cif"
X     file_name = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/COM.cif"
X     file_name = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/tp7_new.cif"
X     file_name = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/F43.cif"
    }
  """
  (d) Change the sequence file_name to the full path of where it is now. It should be in the "tutorial_data" folder

    sequence {
X     file_name = "/Users/mdasgupta/Desktop/MCR_OX/rcsb_pdb_3M1V.fasta"
    }

  (e) Change the prefix to "MCR_XFEL_refine_endrapid" and set serial to 0

  """
  output {
X   prefix = "MCR_REV_refine"
X   serial = 5
  """

  (f) Change the following blob within triple quotes to the following 
 

  --Original:
  """
  main {
X   number_of_macro_cycles = 5
X   nproc = 7
  }
  """

  --Modified:
  """
  main {
    number_of_macro_cycles=3
    nproc = 1
    random_seed=REPLACEME
    ordered_solvent=False
  }
  """


  (g) Set optimizing weights to be False. This is time consuming and should not be done for end/rapid

  """
  target_weights {
X   optimize_xyz_weight = True
X   optimize_adp_weight = True
  }
  """

  (h) Delete the following section marked. Be very careful of the curly braces you delete because exact brace closure is important
  Only one pair of curly braces should be taken out
  
  """
X  gui {
X    base_output_dir = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION"
X    tmp_dir = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/.phenix/tmp"
X    notify_email = "mdasgupta@lbl.gov"
X    skip_kinemage = True
X    phil_file = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/MCR_Ni_S1_COM.edits"
X    phil_file = "/Users/mdasgupta/Desktop/MCR_OX/MANUSCRIPT/REVISION/MCR_Ni_Q147.edits"
X  }
  """

9. Now that you are done changing the relevant parameters, Add the following lines within triple quotes at the end of your MCR_XFEL_shake_0.phil file. This is for shaking the starting model for end/rapid 

""" 
refinement.modify_start_model {
    modify {
      sites {
        atom_selection="not ((resname WAT or resname HOH or resname OOO))"
        shake=0.25
      }
    }
  }
"""

10. Copy the MCR_XFEL_shake_0.phil file to a new file MCR_XFEL_continue_0.phil file. Here you have to make the following change in
    the lines marked by X


  (a) Change the file_name that you had specified in step 8(a) to "MCR_XFEL_refine_endrapid_0.pdb" Note this is the output file that will be from the shake step and we are going to feed it for the refinement in the continue step.

  """
  input {
    pdb {
X     file_name=<whatever you had in step 8(a)>
    } 
  """

  (b) Change serial to 3. The script we have makes use of these numberings to do the shake and continue steps so conforming 
      to these conventions makes it easier.
  """
  output {
    prefix="MCR_XFEL_refine_endrapid"
X   serial=0
  """
  (c) Take out the modify_start_model section that was incorporated in step 9

  """
X  refinement.modify_start_model {
X      modify {
X        sites {
X          atom_selection="not ((resname WAT or resname HOH or resname OOO))"
X          shake=0.25
X        }
X      }
X    }
  """

11. Once you have the 2 phil files (i.e MCR_XFEL_shake_0.phil and MCR_XFEL_continue_0.phil) ready, make a new subfolder in the "phils" folder called "template". Copy over both the phil files to this template folder.

12. Run the script create_all_endrapid_phils.sh as "./create_all_endrapid_phils.sh MCR_XFEL". You might have to enable execute permissions with chmod +x prior to running this script. This script creates 100 phil files for the 100 trials of end/rapid. You should also notice two new folder with these 100 phil files.

13. Now change directory back to the base level "MCR_end_rapid_tutorial". Copy over the following scripts from "tutorial_data/scripts" to the current location - (a) kick_and_refine.csh (b) kick_data_bydiff.com (c) stepwise_phenix_refine.sh (d) submit_kick_MCR_XFEL.csh 

14. First open submit_kick_MCR_XFEL.csh. For new proteins you probably want to rename this file based on the ${tag} being used for the end/rapid. In the file I have commented where the changes need to be made. In the file note the line "setenv flashstate MCR_XFEL". This needs to be changed based on what the ${tag} is. In the next line, you also want to change the path to your "output" folder (created in step 2) 

  Original:
  /net/cci/asmit/Y2K_projects/MCR_end_rapid_Dec2021/output/${flashstate}/kick.$SGE_TASK_ID.out
  
  Modified:
  <Full path to your output directory>/${flashstate}/kick.$SGE_TASK_ID.out

You should also change the job name. See the following 
  Original:
  #$ -o kickboth.out -j y -N MCR_XFEL
  Modified:
  #$ -o kickboth.out -j y -N ${tag}

15. Now open the file "kick_and_refine.csh". I have commented where changes need to be made here in the file itself. First change is the full path of the work_folder and work_folder_2. This will be the full path of the "MCR_end_rapid_tutorial" folder. They are the same folder but the paths are slightly different because of access issues on different servers. Hence work_folder is /net/cci/... while work_folder_2 is /net/cci-filer3/... Just look at the full path you get for the folder (`pwd`) and type in the right path for the 2 variables here.
  Next change the folder where the end/rapid will be run. On the Lewis machine it will be run in the $SCRATCH folder (/net/lewis/scratch1/). Make sure you have created your own workspace in the $SCRATCH folder (for example I have created mine under /net/lewis/scratch1/asmit). That way users can work in their own spaces.
  Once you have created the scratch workspace, create the folder where you will run the end/rapid. For example you can do

  mkdir /net/lewis/scratch1/<username>/end_rapid_tutorial
  
  Once created, change the portion in the "kick_and_refine.csh" script to this path.


16. Once all of this is done, you should be set to run end/rapid. Now source the following to invoke the batch submission environment on Lewis. 
  "source /net/boa/srv/sge-6.1u2/default/common/settings.sh"
  Then run the following "qsub -q all.q@lewis -t 1-100 submit_kick_MCR_XFEL.csh"

  That should submit all your end/rapid jobs. You can check status by typing in qstat -u <username>. You can now log out of Lewis and check back in later.

  
17. The end/rapid jobs can take time (1-2 days) to finish. Once done, you need to gather your results which is the final .pdb file from the individual trials. Navigate to your home directory where you have the "MCR_end_rapid_tutorial". Change directory to the "results" folder. Copy the "gather_all.sh" from "tutorial_data/scripts". Open up the script and change the variable scratch_folder to be where you ran your end/rapid (the same path in step 15). Once you edit this, save the file and run it "./gather_all.sh". Make sure you have done chmod +x before running it. Thisw will copy the final .pdb files from each trial to this folder and create 100 such folders, one for each trial. 

18. You should now have all the output files from end/rapid and you can proceed with doing the statistical analysis you wish to do.

