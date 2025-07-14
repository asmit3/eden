#! /bin/bash

export phil_template=$1
export cycles=`grep "serial=" $phil_template | awk '{ print $1 }' | sed "s:serial=::"`
export cycles_total=$2
export pdb=`grep pdb\"$ $phil_template | grep file_name | awk '{ print $1 }' | sed "s:.*=\"::; s:.pdb\":.pdb:"`
export tag=`grep prefix $phil_template | awk '{ print $1 }' | sed "s:.*=\"::; s:\"::"`
export phil=${tag}_${cycles}.phil
cp $phil_template $phil
#export edits=$3

if [ $((cycles)) -ge 0 ]
  then while [ $((cycles_total-cycles)) -ge 0 ]
    do echo using noedits file
    echo starting refinement cycle $cycles of $cycles_total
    phenix.refine $phil --overwrite > ${tag}_${cycles}.out
    #phenix.refine $phil $edits --overwrite > ${tag}_${cycles}.out
    export new_cycles=$((cycles+3))
    export new_phil=${tag}_${new_cycles}.phil
    export new_pdb="${tag}_${cycles}.pdb"
    sed "s:serial=${cycles}:serial=${new_cycles}:; s:${pdb}:${new_pdb}:" $phil > $new_phil
    export cycles=$new_cycles
    export phil=$new_phil
    export pdb=$new_pdb
  done
fi
