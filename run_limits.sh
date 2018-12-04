source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/2DAlphaEnv.tgz ./
export SCRAM_ARCH=??
scramv1 project CMSSW CMSSW_8_1_0
tar -xzvf 2DAlphaEnv.tgz
rm 2DAlphaEnv.tgz

mkdir tardir; cp tarball.tgz tardir/; cd tardir
tar xzvf tarball.tgz
cp -r * ../CMSSW_8_1_0/src/2DAlphabet/
cd ../CMSSW_8_1_0/src/2DAlphabet/
eval `scramv1 runtime -sh`

echo Do_full_limits.py $*
python Do_full_limits.py $* #-s $1 -r $2 -d $3 -n $4 -j $5
cp ?? ../../../

