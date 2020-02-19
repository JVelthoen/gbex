#!/bin/bash

#################
# SBATCH PARAMS #
#################
QOS='short'
STIME='01:00:00'
MEM=500

################
# TRAIN PARAMS #
################

################
#  TEST PARAMS #
################


# Loop through the above arrays
JID=$(sbatch --qos=$QOS --time=$STIME --mem=$MEM --parsable simulation_files/job_1.sh --output=sim_res_1.rds)
    echo $JID
