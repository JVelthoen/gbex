#!/bin/bash

#################
# SBATCH PARAMS #
#################
QOS='short'
STIME='01:30:00'
MEM=600

################
# TRAIN PARAMS #
################

################
#  TEST PARAMS #
################


# Loop through the above arrays
for ((i=1; i<=22; i++))
do
    JID=$(sbatch --qos=$QOS --time=$STIME --mem=$MEM --parsable simulation_files/job_$i.sh --output=sim_res_$i.rds)
    echo $JID
done

#JID=$(sbatch --qos=$QOS --time=$STIME --mem=$MEM --parsable simulation_files/job_609.sh --output=sim_res_609.rds)
# echo $JID
