#!/bin/bash

# Extract the base name of the input file (without extension) to use for the output file name.
bname=$(basename -s .slcio ${1})
output_filepath="/nfs/dust/ilc/user/silvaand/SM@250_miniDST/2f/${bname}_mini-DST.slcio"

# Load the Key4HEP software environment.
source /cvmfs/sw.hsf.org/key4hep/setup.sh 
# Change to the directory where your MarlinReco is set up.
cd /afs/desy.de/user/s/silvaand/silvaand/workarea/MarlinReco/
k4_local_repo # Initialize local Key4HEP repository.

# Remove any existing output file with the same name to avoid conflicts.
rm -f $output_filepath

# Change to the directory containing the MiniDST production configuration for ILD.
cd /afs/desy.de/user/s/silvaand/silvaand/workarea/ILDConfig/StandardConfig/production

# Run Marlin with the MiniDSTMaker.xml configuration file.
# - `--constant.lcgeo_DIR`: Sets the directory for geometry files.
# - `--global.LCIOInputFiles`: Specifies the input LCIO file.
# - `--constant.OutputFile`: Specifies the output file path.
Marlin MarlinStdRecoMiniDST.xml --constant.lcgeo_DIR=$lcgeo_DIR --global.LCIOInputFiles=${1} --constant.OutputFile=$output_filepath
