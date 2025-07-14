state=$1

mkdir ${state}_iters_shake
mkdir ${state}_iters_continue

for i in `seq 1 100`
do
  cp template/${state}_shake_0.phil ${state}_iters_shake/${state}_shake_${i}.phil 
  sed -i "s/REPLACEME/${i}/g" ${state}_iters_shake/${state}_shake_${i}.phil 
  cp template/${state}_continue_0.phil ${state}_iters_continue/${state}_continue_${i}.phil 
  sed -i "s/REPLACEME/${i}/g" ${state}_iters_continue/${state}_continue_${i}.phil 
done
