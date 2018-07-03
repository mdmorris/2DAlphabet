#####################################################################################
# This is a standalone script that will put together a per-bin post-fit/pre-fit		#
# ratio of the ttbar MC shapes. It is NOT generic and so should only really be 		#
# used as a template. It will also quickly be outdated once a new generic script	#
# has been written that re-centers/fixes the nuisance parameters in the signal		#
# using the values found in the ttbar sideband fit.									#
#####################################################################################

import ROOT
from ROOT import *

import math
from math import sqrt
import array
from array import array

def copyHistWithNewXbounds(thisHist,copyName,newBinWidthX,xNewBinsLow,xNewBinsHigh):
    # Make a copy with the same Y bins but new X bins
    nBinsY = thisHist.GetNbinsY()
    yBinsLow = thisHist.GetYaxis().GetXmin()
    yBinsHigh = thisHist.GetYaxis().GetXmax()
    nNewBinsX = int((xNewBinsHigh-xNewBinsLow)/float(newBinWidthX))
    histCopy = TH2F(copyName,copyName,nNewBinsX,xNewBinsLow,xNewBinsHigh,nBinsY,yBinsLow,yBinsHigh)
    histCopy.Sumw2()
    
    histCopy.GetXaxis().SetName(thisHist.GetXaxis().GetName())
    histCopy.GetYaxis().SetName(thisHist.GetYaxis().GetName())


    # Loop through the old bins
    for binY in range(1,nBinsY+1):
        # print 'Bin y: ' + str(binY)
        for newBinX in range(1,nNewBinsX+1):
            newBinContent = 0
            newBinErrorSq = 0
            newBinXlow = histCopy.GetXaxis().GetBinLowEdge(newBinX)
            newBinXhigh = histCopy.GetXaxis().GetBinUpEdge(newBinX)

            # print '\t New bin x: ' + str(newBinX) + ', ' + str(newBinXlow) + ', ' + str(newBinXhigh)
            for oldBinX in range(1,thisHist.GetNbinsX()+1):
                if thisHist.GetXaxis().GetBinLowEdge(oldBinX) >= newBinXlow and thisHist.GetXaxis().GetBinUpEdge(oldBinX) <= newBinXhigh:
                    # print '\t \t Old bin x: ' + str(oldBinX) + ', ' + str(thisHist.GetXaxis().GetBinLowEdge(oldBinX)) + ', ' + str(thisHist.GetXaxis().GetBinUpEdge(oldBinX))
                    # print '\t \t Adding content ' + str(thisHist.GetBinContent(oldBinX,binY))
                    newBinContent += thisHist.GetBinContent(oldBinX,binY)
                    newBinErrorSq += thisHist.GetBinError(oldBinX,binY)**2

            # print '\t Setting content ' + str(newBinContent) + '+/-' + str(sqrt(newBinErrorSq))
            histCopy.SetBinContent(newBinX,binY,newBinContent)
            histCopy.SetBinError(newBinX,binY,sqrt(newBinErrorSq))

    return histCopy


def copyHistWithNewYbounds(thisHist,copyName,newBinWidthY,yNewBinsLow,yNewBinsHigh):
    # Make a copy with the same Y bins but new X bins
    nBinsX = thisHist.GetNbinsX()
    xBinsLow = thisHist.GetXaxis().GetXmin()
    xBinsHigh = thisHist.GetXaxis().GetXmax()
    nNewBinsY = int((yNewBinsHigh-yNewBinsLow)/float(newBinWidthY))
    histCopy = TH2F(copyName,copyName,nBinsX,xBinsLow,xBinsHigh,nNewBinsY,yNewBinsLow,yNewBinsHigh)
    histCopy.Sumw2()
    
    histCopy.GetXaxis().SetName(thisHist.GetXaxis().GetName())
    histCopy.GetYaxis().SetName(thisHist.GetYaxis().GetName())


    # Loop through the old bins
    for binX in range(1,nBinsX+1):
        # print 'Bin y: ' + str(binY)
        for newBinY in range(1,nNewBinsY+1):
            newBinContent = 0
            newBinErrorSq = 0
            newBinYlow = histCopy.GetYaxis().GetBinLowEdge(newBinY)
            newBinYhigh = histCopy.GetYaxis().GetBinUpEdge(newBinY)

            # print '\t New bin x: ' + str(newBinX) + ', ' + str(newBinXlow) + ', ' + str(newBinXhigh)
            for oldBinY in range(1,thisHist.GetNbinsY()+1):
                if thisHist.GetYaxis().GetBinLowEdge(oldBinY) >= newBinYlow and thisHist.GetYaxis().GetBinUpEdge(oldBinY) <= newBinYhigh:
                    # print '\t \t Old bin x: ' + str(oldBinX) + ', ' + str(thisHist.GetXaxis().GetBinLowEdge(oldBinX)) + ', ' + str(thisHist.GetXaxis().GetBinUpEdge(oldBinX))
                    # print '\t \t Adding content ' + str(thisHist.GetBinContent(oldBinX,binY))
                    newBinContent += thisHist.GetBinContent(binX,oldBinY)
                    newBinErrorSq += thisHist.GetBinError(binX,oldBinY)**2

            # print '\t Setting content ' + str(newBinContent) + '+/-' + str(sqrt(newBinErrorSq))
            histCopy.SetBinContent(binX,newBinY,newBinContent)
            histCopy.SetBinError(binX,newBinY,sqrt(newBinErrorSq))

    return histCopy


def rebinY(thisHist,name,new_y_bins_array):
    xnbins = thisHist.GetNbinsX()
    xmin = thisHist.GetXaxis().GetXmin()
    xmax = thisHist.GetXaxis().GetXmax()

    rebinned = TH2F(name,name,xnbins,xmin,xmax,len(new_y_bins_array)-1,new_y_bins_array)
    rebinned.Sumw2()

    for xbin in range(1,xnbins+1):
        newBinContent = 0
        newBinErrorSq = 0
        rebinHistYBin = 1
        for ybin in range(1,thisHist.GetNbinsY()+1):
            # If upper edge of old Rpf ybin is < upper edge of rebinHistYBin then add the Rpf bin to the count
            if thisHist.GetYaxis().GetBinUpEdge(ybin) < rebinned.GetYaxis().GetBinUpEdge(rebinHistYBin):
                newBinContent += thisHist.GetBinContent(xbin,ybin)
                newBinErrorSq += thisHist.GetBinError(xbin,ybin)**2
            # If ==, add to newBinContent, assign newBinContent to current rebinHistYBin, move to the next rebinHistYBin, and restart newBinContent at 0
            elif thisHist.GetYaxis().GetBinUpEdge(ybin) == rebinned.GetYaxis().GetBinUpEdge(rebinHistYBin):
                newBinContent += thisHist.GetBinContent(xbin,ybin)
                newBinErrorSq += thisHist.GetBinError(xbin,ybin)**2
                rebinned.SetBinContent(xbin, rebinHistYBin, newBinContent)
                rebinned.SetBinError(xbin, rebinHistYBin, sqrt(newBinErrorSq))# NEED TO SET BIN ERRORS
                rebinHistYBin += 1
                newBinContent = 0
                newBinErrorSq = 0
            else:
                print 'ERROR when doing psuedo-2D fit approximation. Slices do not line up on y bin edges'
                quit()

    return rebinned


if __name__ == "__main__":
    file_prefit = TFile.Open("/uscms_data/d3/lcorcodi/TTbar13TeV/CMSSW_7_4_1/src/TTbar13TeV/rootfiles/TT2Dalphabetweightedttbar_Trigger_nominal_none_PSET_default.root")

    new_y_bins = array('d',[800.0, 1100.0, 1300.0, 1600.0, 2000.0, 2800.0, 4000.0])

    pass_prefit_oldbins = file_prefit.Get('MttvMtPass')
    fail_prefit_oldbins = file_prefit.Get('MttvMtFail')
    # pass_prefit_x = copyHistWithNewXbounds(pass_prefit_oldbins,'pass_prefit_x',20,50,330)
    # fail_prefit_x = copyHistWithNewXbounds(fail_prefit_oldbins,'fail_prefit_x',20,50,330)
    pass_prefit = copyHistWithNewYbounds(pass_prefit_oldbins,'pass_prefit',100,800,4000)
    fail_prefit = copyHistWithNewYbounds(fail_prefit_oldbins,'fail_prefit',100,800,4000)

    file_postfit = TFile.Open('tt_xpol2_ypol2_widesig_nosig/MaxLikelihoodFitResult.root')
    w_postfit = file_postfit.Get('MaxLikelihoodFitResult')

    x_var = w_postfit.var('jetmass')
    y_var = w_postfit.var('resmass')

    x_nbins = 60
    x_low = 50
    x_high = 350

    y_nbins = 32
    y_low = 800
    y_high = 4000

    pass_postfit_unit = w_postfit.pdf('shapeBkg_pass_ttbar_morph')
    fail_postfit_unit = w_postfit.pdf('shapeBkg_fail_ttbar_morph')
    pass_postfit_norm = w_postfit.function('n_exp_final_binpass_proc_ttbar')
    fail_postfit_norm = w_postfit.function('n_exp_final_binfail_proc_ttbar')

    pass_postfit = pass_postfit_unit.createHistogram('ttbar_pass',x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)))
    fail_postfit = fail_postfit_unit.createHistogram('ttbar_fail',x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)))

    pass_postfit.Scale(pass_postfit_norm.getValV())
    fail_postfit.Scale(fail_postfit_norm.getValV())

    # for hist in [pass_prefit,pass_postfit]:
    #     print hist.GetNbinsX()
    #     print hist.GetNbinsY()
    #     print hist.GetYaxis().GetXmin()
    #     print hist.GetYaxis().GetXmax()

    pass_prefit_py = pass_prefit.ProjectionY()
    fail_prefit_py = fail_prefit.ProjectionY()
    pass_postfit_py = pass_postfit.ProjectionY()
    fail_postfit_py = fail_postfit.ProjectionY()

    pass_postfit_err = pass_postfit_py.Clone(pass_postfit_py.GetName()+'_err')
    fail_postfit_err = fail_postfit_py.Clone(fail_postfit_py.GetName()+'_err')



    # 1D res mass projections
    pass_prefit_py.SetLineColor(kBlack)
    pass_prefit_py.SetMarkerStyle(8)
    fail_prefit_py.SetLineColor(kBlack)
    fail_prefit_py.SetMarkerStyle(8)

    pass_postfit_py.SetFillColor(kRed)
    fail_postfit_py.SetFillColor(kRed)


    pass_postfit_err.SetFillColor(kBlack)
    pass_postfit_err.SetFillStyle(3354)
    fail_postfit_err.SetFillColor(kBlack)
    fail_postfit_err.SetFillStyle(3354)

    passCan = TCanvas('pass','pass',700,600)
    passCan.cd()
    pass_prefit_py.Draw('pe')
    pass_postfit_py.Draw('same hist')
    pass_postfit_err.Draw('same e2')

    failCan = TCanvas('fail','fail',700,600)
    failCan.cd()
    fail_prefit_py.Draw('pe')
    fail_postfit_py.Draw('same hist')
    fail_postfit_err.Draw('same e2')


    # Ratio 2D
    passRatio2D = pass_postfit.Clone('passRatio2D')
    failRatio2D = fail_postfit.Clone('failRatio2D')

    passRatio2D.Divide(pass_prefit)
    failRatio2D.Divide(fail_prefit)


    # Ratio 1D res
    r1Dres = TCanvas('r1Dres','r1Dres',1000,600)
    r1Dres.Divide(2,1)

    passRatio1Dres = passRatio2D.ProjectionY()
    failRatio1Dres = failRatio2D.ProjectionY()

    r1Dres.cd(1)
    passRatio1Dres.GetYaxis().SetRangeUser(0,15)
    passRatio1Dres.Draw('pe')

    r1Dres.cd(2)
    failRatio1Dres.GetYaxis().SetRangeUser(0,15)
    failRatio1Dres.Draw('pe')

    # Ratio 1D jet
    r1Djet = TCanvas('r1Djet','r1Djet',1000,600)
    r1Djet.Divide(2,1)

    passRatio1Djet = passRatio2D.ProjectionX()
    failRatio1Djet = failRatio2D.ProjectionX()

    r1Djet.cd(1)
    passRatio1Djet.GetYaxis().SetRangeUser(0,25)
    passRatio1Djet.Draw('pe')

    r1Djet.cd(2)
    failRatio1Djet.GetYaxis().SetRangeUser(0,25)
    failRatio1Djet.Draw('pe')

    passRatio2D.SaveAs('ttSF/passRatio2D.root','root')
    failRatio2D.SaveAs('ttSF/failRatio2D.root','root')

    raw_input('waiting')
    # # Fit the 1D ratios
    # xformula = TF1('xformula','[0]',x_low,x_high)
    # passRatio1Djet.Fit('xformula')
    # passXjet = passRatio1Djet.GetFunction('xformula')
    # failRatio1Djet.Fit('xformula')
    # failXjet = failRatio1Djet.GetFunction('xformula')
    
    # print passXjet.GetParameter(0)
    # print failXjet.GetParameter(0)


    # yformula = TF1('yformula','[0]',y_low,y_high)
    # passRatio1Dres.Fit('yformula')
    # passYres = passRatio1Dres.GetFunction('yformula')
    # failRatio1Dres.Fit('yformula')
    # failYres = failRatio1Dres.GetFunction('yformula')
    
    # print passYres.GetParameter(0)
    # print failYres.GetParameter(0)

    # r1Djet.cd(1)
    # passRatio1Djet.Draw('pe')

    # r1Djet.cd(2)
    # failRatio1Djet.Draw('pe')


    # # Fit the 2D ratios
    # formula2Dpass = TF2('formula2Dpass','[0]',x_low,x_high,y_low,y_high)
    # formula2Dfail = TF2('formula2Dfail','[0]',x_low,x_high,y_low,y_high)


    # passRatio2D.Fit('formula2Dpass')
    # pass2Dresult = passRatio2D.GetFunction('formula2Dpass')
    # failRatio2D.Fit('formula2Dfail')
    # fail2Dresult = failRatio2D.GetFunction('formula2Dfail')

    # for p in range(1):
    #     print pass2Dresult.GetParameter(p)

    # for p in range(1):
    #     print fail2Dresult.GetParameter(p)

    # r2D = TCanvas('r2D','r2D',1000,600)

    # r2D.Divide(2,1)
    # r2D.cd(1)
    # passRatio2D.Draw('lego')

    # r2D.cd(2)
    # failRatio2D.Draw('lego')

    # raw_input('waiting')