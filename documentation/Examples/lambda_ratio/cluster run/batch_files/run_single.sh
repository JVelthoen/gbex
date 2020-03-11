#!/bin/bash

#################
# SBATCH PARAMS #
#################
QOS='short'
STIME='02:00:00'
MEM=1000

################
# TRAIN PARAMS #
################

################
#  TEST PARAMS #
################


# Loop through the above arrays
JID=$(sbatch --qos=$QOS --time=$STIME --mem=$MEM --parsable simulation_files/job_1.sh --output=sim_res_1.rds)
    echo $JID
