import ROOT
from ROOT import *

f = TFile.Open('../dataSidebandFewerbins_tt/fitDiagnostics.root')

rfr_old = f.Get('fit_s')
ral_old = rfr_old.floatParsFinal()
ral_new = RooArgList()
iterator = ral_old.createIterator()
var = iterator.Next()

while var:
    if var.GetName() in ['polyCoeff_x0y0','polyCoeff_x1y0','polyCoeff_x2y0','polyCoeff_x0y1','polyCoeff_x1y1','polyCoeff_x2y1','polyCoeff_x0y2','polyCoeff_x1y2','polyCoeff_x2y2']:
        new_param = var.Clone()
        new_param.SetName(var.GetName())
        new_param.setError(var.getError()*100)
        ral_new.add(new_param)

    else:
        ral_new.add(var)

    var = iterator.Next()

rfr_new = RooFitResult('fit_s_new')
rfr_new._RooFitResult__setInitParList(rfr_old.floatParsInit())
rfr_new._RooFitResult__setFinalParList(ral_new)

f_new = TFile.Open('test_RFR.root')
rfr_new.Write()
f_new.Close()
