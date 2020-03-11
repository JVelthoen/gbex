#!/bin/bash

cd /home/nfs/jjvelthoen/

#cd /Users/jjvelthoen/Documents/simulations variable selection/sf_model_size


mkdir simulation_files
NRSIM=500

for ((i=1; i<=NRSIM; i++))
do
    cp script.R simulation_files/script_$i.R
    sed -i "4s/.*/sim_nr <- $i/" simulation_files/script_$i.R

    cp job.sh simulation_files/job_$i.sh
    sed -i "s/script.R/script_$i.R/g" simulation_files/job_$i.sh
done
