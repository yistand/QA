#!/bin/bash -l
#PBS -l nodes=1:ppn=1,walltime=100:00:00
#PBS -r n
#PBS -V
#PBS -q hep
#PBS -j oe


cd $PBS_O_WORKDIR
echo $PWD
echo "Job Start at `date`"

#module load Apps/ROOT/5.34.34
echo source /home/hep/caines/ly247/.bashrc
source /home/hep/caines/ly247/.bashrc


echo "root -l -q -b BemcMatchRate.C++"
root -l -q -b BemcMatchRate.C++

echo "Job End at `date`"




