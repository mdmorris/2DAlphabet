cp RooParametricHist2D.cxx $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/src/
cp RooParametricHist2D.h $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/interface/
cp classes_def.xml $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/src/
cp classes.h $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/src/

cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b
cd $CMSSW_BASE/2DAlphabet