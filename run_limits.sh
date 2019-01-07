source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/lcorcodi/2Dalpha.tgz ./
export SCRAM_ARCH=slc6_amd64_gcc530
scramv1 project CMSSW CMSSW_8_1_0
tar xzf 2Dalpha.tgz
rm 2Dalpha.tgz

mkdir tardir; cp tarball.tgz tardir/; cd tardir
tar xzf tarball.tgz
cp -r * ../CMSSW_8_1_0/src/2DAlphabet/
cd ../CMSSW_8_1_0/src/2DAlphabet/
eval `scramv1 runtime -sh`

echo Sideband_wrapper.py $*
python Sideband_wrapper.py $* #-s $1 -r $2 -d $3 -n $4 -j $5
cp *.Asymptotic.* ../../../

