#!/bin/bash

# Comment
bname=$(basename -s .slcio ${1})
output_filepath="/nfs/dust/ilc/user/silvaand/4f/4f_ZZnunu/${bname}_mini-DST.slcio"

source /cvmfs/sw.hsf.org/key4hep/setup.sh
cd /afs/desy.de/user/s/silvaand/silvaand/workarea/MarlinReco/
k4_local_repo

rm -f $output_filepath

cd /afs/desy.de/user/s/silvaand/silvaand/workarea/ILDConfig/StandardConfig/production

Marlin MarlinStdRecoMiniDST.xml --constant.lcgeo_DIR=$lcgeo_DIR --global.LCIOInputFiles=${1} --constant.OutputFile=$output_filepath
