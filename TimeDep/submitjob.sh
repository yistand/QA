#!/bin/bash
#
#BSUB -W 24:00



echo $PWD
echo "Job Start at `date`"

echo source /home/fas/caines/ly247/.bash_profile
source /home/fas/caines/ly247/.bash_profile


root -l <<EOF
gSystem->Load("libHist.so");
gSystem->Load("libPhysics"); 
gSystem->Load("~/Software/PicoCode/eventStructuredAu/libTStarJetPico.so");
.L TimeDep.C+
TimeDep()
.q
EOF

echo "Job End at `date`"

exit 0


