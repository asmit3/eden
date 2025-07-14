#! /bin/bash

# Change this to be the folder where you ran your end/rapid
scratch_folder=scratch1/asmit/end_rapid/MCR_end_rapid_Dec2021/kicked_data_v1/
servers=(lewis)


for server in ${servers[@]}
  do
    for f in `ls /net/${server}/${scratch_folder}/`;
      do if [ -f /net/${server}/${scratch_folder}/"$f"/*12.pdb ];
        then if [ ! -d "$f" ];
          then mkdir "$f"
          fi
        if [ ! -f "$f"/*12.pdb ];
          then cp /net/${server}/${scratch_folder}/"$f"/*12.pdb "$f"
          fi
        fi
      done
  done
