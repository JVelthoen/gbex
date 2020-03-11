#!/bin/sh

# You can control the resources and scheduling with '#SBATCH' settings
# (see 'man sbatch' for more information on setting these parameters)

# The default partition is the 'general' partition
#SBATCH --partition=general

# The default number of parallel tasks per job is 1
#SBATCH --ntasks=1

# The default number of CPUs per task is 1, however CPUs are always allocated per 2, so for a single task you should use 2
#SBATCH --cpus-per-task=4

# Set mail type to 'END' to receive a mail when the job finishes (with usage statistics)
#SBATCH --mail-type=FAIL


#module use /opt/insy/modulefiles
# Uncomment these lines when your job requires this software
#module use /opt/insy/modulefiles
#module load cuda/10.0 cudnn/10.0-7.3.0.29
#module load matlab/R2018b



# Complex or heavy commands should be started with 'srun' (see 'man srun' for more information)
# (This is just an example, srun is of course not necessary for this command.)
Rscript /home/nfs/jjvelthoen/simulation_files/script.R
