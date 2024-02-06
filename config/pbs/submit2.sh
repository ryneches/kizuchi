#!/bin/bash

##################
####  PBS preamble

#### #### ####  These are the most frequently changing options

####  Job name
#PBS -N snakejob

####  Request resources here
####    These are typically, number of processors, amount of memory,
####    an the amount of time a job requires.  May include processor
####    type, too.

#PBS -l select=1:ncpus=1:mem=18gb
#PBS -l walltime=96:00:00

####  Flux account and queue specification here
####    These will change if you work on multiple projects, or need
####    special hardware, like large memory nodes or GPUs or,
####    or if you use software that is restricted to campus use.

#### PBS -A ACCOUNT
#PBS -q APC

#### #### ####  These are the least frequently changing options

####  Your e-mail address and when you want e-mail

#PBS -M hpc@vort.org
#PBS -m bea

####  Join output and error; pass environment to job

#PBS -e logs/pbs/
#PBS -o logs/pbs/
#PBS -V

# Add a note here to say what software modules should be loaded.
# for this job to run successfully.
# It will be convenient if you give the actual load command(s), e.g.,
#
# module load intel/16.0.4
# module load sratoolkit/2.8.2-1

####  End PBS preamble
##################

####  PBS job only tasks

##  Print the nodename(s) to the output in case needed for diagnostics,
##  or if you need information about the hardware after the job ran.
if [ -e "$PBS_NODEFILE" ] ; then
    echo "Running on"
    uniq -c $PBS_NODEFILE
fi

# if running a job, tell which machine its running on

##  Change to the directory from which you submit the job, if running
##  from within a job
if [ -d "$PBS_O_WORKDIR" ] ; then
    cd $PBS_O_WORKDIR
fi

####  Commands your job should run follow this line
##
##  Note:  In batch jobs, programs should always run in foreground.  Do
##         not use an & at the end of a command. Bad things will happen.

#### Create cluster log directory

mkdir -p logs
mkdir -p logs/pbs

# Initiating snakemake and running workflow in cluster mode
snakemake                                \
    --snakefile workflow/Snakefile2      \
    --profile config/pbs                 \
    --group-components phylogenetics=512 \
    --rerun-incomplete                   \

# Printing out job summary
qstat -f $PBS_JOBID

##  If you copied any files to /tmp, make sure you delete them here!
