#!/bin/bash
# call it by "source SetEnvironment.sh"

#ly 2015.07.09 modified from KKauder's SetEnvironment.csh
#ly set environment variables for linking library


# CHANGE the following to suit your environment
declare -x BASEDIR=/home/fas/caines/ly247/Software

#### ROOT
#declare -x ROOTSYS=/home/hep/share/app/root	# already set up in ~/.barshrc


### TStarJetPicoDst structure
declare -x STARPICOPATH=${BASEDIR}/PicoCode/eventStructuredAu


### underlying jet code
declare -x PPJET=/home/fas/caines/ly247/code/ppjet




############## Done with indivivual settings.

###### Update paths
if [[ -z $LD_LIBRARY_PATH ]] 
 then
	declare -x LD_LIBRARY_PATH
fi
if [[ -z $DYLD_LIBRARY_PATH ]] 
 then
	declare -x DYLD_LIBRARY_PATH
fi

declare -x PATH=./bin:${ROOTSYS}/bin:${PPJET}/include:${PATH}
declare -x CPATH=./bin:${ROOTSYS}/bin:${STARPICOPATH}:${PPJET}/include:${CPATH}	#ly
declare -x LD_LIBRARY_PATH=${ROOTSYS}/lib:${STARPICOPATH}:${PPJET}/src:${LD_LIBRARY_PATH}
declare -x DYLD_LIBRARY_PATH=${ROOTSYS}/lib:${STARPICOPATH}:${PPJET}/src:${DYLD_LIBRARY_PATH}

# ly if [[ $?TERM == 0 || $?prompt == 0 ]] exit 0

echo ''
echo 'Setup ktJet, STARJetPico '
echo '====================================='
echo ''
echo "<I>---------------Info--------------------<I>"
echo "Setting up the following environments: "
echo "ROOT: " $ROOTSYS
echo "STARPICOPATH: " $STARPICOPATH
echo "<I>---------------Info--------------------<I>"
echo ""
