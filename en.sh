#!/bin/bash
# parallel Environment JOB
#
#$ -cwd
#$ -j y
#$ -S /bin/bash

. ~/.bashrc
unset SGE_ROOT
# unset SGE_ROOT is a fix for
# error: executing task of job X failed: execution daemon on host Y didn't accept task

# NSLOTS and TMP/machines is generated by SGE Paralell Environment
ulimit -c unlimited
ulimit -s unlimited
gfortran -O3 -ffree-form h.f
./a.out
epstopdf energy.eps
