cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/Observation.h /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/interface/
cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/CombineHarvester.h /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/interface/
cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/Process.h /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/interface/

cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/Process.cc /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/src/
cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/Observation.cc /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/src/
cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/CombineHarvester_Evaluate.cc /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/src/
cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/CombineHarvester.cc /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/src
cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/PostFitShapes2D.cpp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/bin/

cp /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/setup/BuildFile.xml /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/bin/

cd /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_8_1_0/src/CombineHarvester/CombineTools/
scram b clean; scram b -j 10 USER_CXXFLAGS="-g"
cd /uscms_data/d3/lcorcodi/BStar13TeV/CMSSW_7_4_7_patch2/src/2DAlphabet/

