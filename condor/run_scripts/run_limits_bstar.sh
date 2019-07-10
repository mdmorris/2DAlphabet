#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/2Dalpha.tgz ./
export SCRAM_ARCH=slc6_amd64_gcc700
scramv1 project CMSSW CMSSW_10_2_13
tar -xzf 2Dalpha.tgz
rm 2Dalpha.tgz

mkdir tardir; cp tarball.tgz tardir/; cd tardir
tar -xzf tarball.tgz
cp -r * ../CMSSW_10_2_13/src/2DAlphabet/
cd ../CMSSW_10_2_13/src/2DAlphabet/
eval `scramv1 runtime -sh`
scramv1 b clean; scramv1 b

# Grab the rootfiles
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/bstar_presel_rootfiles.tgz ./
mkdir rootfiles
tar -xzf bstar_presel_rootfiles.tgz
rm bstar_presel_rootfiles.tgz

echo run_Limit.py $*
python run_Limit.py $* 
sh ./shell_finisher.sh
