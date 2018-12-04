# For Combine (fitting)
cp chebyshevBasis.cxx $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/src/
cp chebyshevBasis.h $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/interface/
cp RooParametricHist2D.cxx $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/src/
cp RooParametricHist2D.h $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/interface/
cp ShapeTools.py $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/python/
cp BuildFile.xml $CMSSW_BASE/src/CombineHarvester/CombineTools/bin/

# For Combine Harvester (plotting)
cp Observation.h $CMSSW_BASE/src/CombineHarvester/CombineTools/interface/
cp CombineHarvester.h $CMSSW_BASE/src/CombineHarvester/CombineTools/interface/
cp Process.h $CMSSW_BASE/src/CombineHarvester/CombineTools/interface/
cp Process.cc $CMSSW_BASE/src/CombineHarvester/CombineTools/src/
cp Observation.cc $CMSSW_BASE/src/CombineHarvester/CombineTools/src/
cp CombineHarvester_Evaluate.cc $CMSSW_BASE/src/CombineHarvester/CombineTools/src/
cp CombineHarvester.cc $CMSSW_BASE/src/CombineHarvester/CombineTools/src
cp PostFitShapes2D.cpp $CMSSW_BASE/src/CombineHarvester/CombineTools/bin/


# Compile Combine
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/src/

if [ `grep "RooParametricHist2D" classes.h | wc -l` == 0 ] 
then
	sed -i '/RooParametricHist.h/ a\#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist2D.h"' classes.h
fi

if [ `grep "RooParametricHist2D" classes_def.xml | wc -l` == 0 ] 
then
	sed -i '/RooParametricHist/ a\        <class name="RooParametricHist2D" />' classes_def.xml
fi

if [ `grep "chebyshevBasis" classes.h | wc -l` == 0 ] 
then
	sed -i '/RooParametricHist2D.h/ a\#include "HiggsAnalysis/CombinedLimit/interface/chebyshevBasis.h"' classes.h
fi

if [ `grep "chebyshevBasis" classes_def.xml | wc -l` == 0 ] 
then
	sed -i '/RooParametricHist2D/ a\        <class name="chebyshevBasis" />' classes_def.xml
fi


# Compile Combine Harvester
cd $CMSSW_BASE/src/
scram b clean; scram b -j 10 USER_CXXFLAGS="-g"

# Come back
cd $CMSSW_BASE/src/2DAlphabet
