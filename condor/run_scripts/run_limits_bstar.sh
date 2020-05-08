#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/10XwithNano.tgz ./
export SCRAM_ARCH=slc6_amd64_gcc700
scramv1 project CMSSW CMSSW_10_2_13
tar -xzf 10XwithNano.tgz
rm 10XwithNano.tgz

echo TEMPTAR
mkdir tardir; cp TEMPTAR tardir/; cd tardir
tar -xzvf TEMPTAR
cp -r * ../CMSSW_10_2_13/src/2DAlphabet/
cd ../CMSSW_10_2_13/src/2DAlphabet/
eval `scramv1 runtime -sh`
#scramv1 b clean; scramv1 b

# Grab the rootfiles
mkdir rootfiles
cd rootfiles
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/bstar_presel_rootfiles.tgz ./
tar -xzf bstar_presel_rootfiles.tgz
rm bstar_presel_rootfiles.tgz
cd ../

ls

echo run_Limit.py $*
python run_Limit.py $* 
sh ./shell_finisher.sh
