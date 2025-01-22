#!/bin/bash

# Comment
output_filepath="/nfs/dust/ilc/user/silvaand/PreSelection_logs/2f/ROOTFiles/PreSelection2f.root"

source /cvmfs/sw.hsf.org/key4hep/setup.sh 
cd /afs/desy.de/user/s/silvaand/silvaand/workarea/MarlinReco/
k4_local_repo

rm -f $output_filepath

cd /afs/desy.de/user/s/silvaand/silvaand/workarea/ILDConfig/StandardConfig/production

Marlin PreSelectionSteering.xml --global.LCIOInputFiles=${1} --constant.TTreeFileName=$output_filepath
