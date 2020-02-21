#!/bin/bash

cd /home/nfs/jjvelthoen/

#cd /Users/jjvelthoen/Documents/simulations variable selection/sf_model_size


mkdir simulation_files
NRSIM=30

for ((i=1; i<=NRSIM; i++))
do
    cp script.R simulation_files/script_$i.R
    sed -i "7s/.*/sim_nr <- $i/" simulation_files/script_$i.R
    START=$(((i-1)*2 + 1))
    END=$((i*2))
    sed -n "$START,$END p" parameters.txt > temp.txt
    sed -i '10,11d' simulation_files/script_$i.R
    sed -i '9 r temp.txt' simulation_files/script_$i.R

    cp job.sh simulation_files/job_$i.sh
    sed -i "s/script.R/script_$i.R/g" simulation_files/job_$i.sh
done

rm temp.txt
