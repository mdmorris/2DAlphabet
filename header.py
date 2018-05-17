import ROOT
from ROOT import *
import os
import math
from math import sqrt


# Function stolen from https://stackoverflow.com/questions/9590382/forcing-python-json-module-to-work-with-ascii
def ascii_encode_dict(data):    
    ascii_encode = lambda x: x.encode('ascii') if isinstance(x, unicode) else x 
    return dict(map(ascii_encode, pair) for pair in data.items())

def getRRVs(inputConfig, blinded):
    xname = inputConfig['BINNING']['X']['NAME']
    xlowname = xname + '_Low'
    xhighname = xname + '_High'
    yname = inputConfig['BINNING']['Y']['NAME']

    xlow = inputConfig['BINNING']['X']['LOW']
    xhigh = inputConfig['BINNING']['X']['HIGH']

    ylow = inputConfig['BINNING']['Y']['LOW']
    yhigh = inputConfig['BINNING']['Y']['HIGH']

    xRRV = RooRealVar(xname,xname,xlow,xhigh)
    yRRV = RooRealVar(yname,yname,ylow,yhigh)

    if blinded:
        xSigStart = inputConfig['BINNING']['X']['SIGSTART']
        xSigEnd = inputConfig['BINNING']['X']['SIGEND']
        xLowRRV = RooRealVar(xlowname,xlowname,xlow,xSigStart)
        xHighRRV = RooRealVar(xhighname,xhighname,xSigEnd,xhigh)

        return xRRV,xLowRRV,xHighRRV,yRRV
    
    else:
        return xRRV,yRRV

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

def makeBlindedHist(nomHist,lowHist,highHist):
    # Grab stuff to make it easier to read
    xlow = nomHist.GetXaxis().GetXmin()
    xhigh = nomHist.GetXaxis().GetXmax()
    xnbins = nomHist.GetNbinsX()
    ylow = nomHist.GetYaxis().GetXmin()
    yhigh = nomHist.GetYaxis().GetXmax()
    ynbins = nomHist.GetNbinsY()
    blindName = nomHist.GetName()

    # Need to change nominal hist name or we'll get a memory leak
    nomHist.SetName(blindName+'_unblinded')

    blindedHist = TH2F(blindName,blindName,xnbins,xlow,xhigh,ynbins,ylow,yhigh)

    for binY in range(1,ynbins+1):
        # First fill with the lowHist bins since this is easy
        for binX in range(1,lowHist.GetNbinsX()+1):
            blindedHist.SetBinContent(binX,binY,lowHist.GetBinContent(binX,binY))
       
        # Now figure out how many bins we have to jump to get to highHist
        bins2jump = xnbins - highHist.GetNbinsX()
        for binX in range(1,highHist.GetNbinsX()+1):
            blindedHist.SetBinContent(binX+bins2jump,binY,highHist.GetBinContent(binX,binY))

    return blindedHist

def makeRDH(myTH2,RAL_vars):
    name = myTH2.GetName()
    thisRDH = RooDataHist(name,name,RAL_vars,myTH2)
    return thisRDH


def makeRHP(myRDH,RAL_vars):
    name = myRDH.GetName()
    thisRAS = RooArgSet(RAL_vars)
    thisRHP = RooHistPdf(name,name,thisRAS,myRDH)
    return thisRHP


def colliMate(myString,width=18):
    sub_strings = myString.split(' ')
    new_string = ''
    for i,sub_string in enumerate(sub_strings):
        string_length = len(sub_string)
        n_spaces = width - string_length
        if i != len(sub_strings)-1:
            if n_spaces <= 0:
                n_spaces = 2
            new_string += sub_string + ' '*n_spaces
        else:
            new_string += sub_string
    return new_string

def dictStructureCopy(inDict):
    newDict = {}
    for k1,v1 in inDict.items():
        if type(v1) == dict:
            newDict[k1] = dictStructureCopy(v1)
        else:
            newDict[k1] = 0
    return newDict

def dictCopy(inDict):
    newDict = {}
    for k1,v1 in inDict.items():
        if type(v1) == dict:
            newDict[k1] = dictStructureCopy(v1)
        else:
            newDict[k1] = v1
    return newDict

def printWorkspace(myfile,myworkspace):
    myf = TFile.Open(myfile)
    myw = myf.Get(myworkspace)
    myw.Print()


# def make4x4Can(name,h1,h2,h3,h4):
#     histlist = [h1,h2,h3,h4]
#     my4x4Can = TCanvas(name,name,1200,1000)
#     my4x4Can.Divide(2,2)

#     for index, hist in enumerate(histlist):
#         my4x4Can.cd(index+1)
#         if hist.ClassName().find('TH2') != -1:
#             hist.Draw('lego')
#         else:
#             hist.Draw('hist e')

#     my4x4Can.Print('plots/'+name+'.pdf','pdf')

# def make1x1Can(name,h1):
#     my1x1Can = TCanvas(name,name,800,700)
#     my1x1Can.cd()
#     h1.Draw('lego')
#     my1x1Can.Print('plots/'+name+'.pdf','pdf')

def makeCan(name, tag, histlist, bkglist=[],logy=False):  
    # histlist is just the generic list but if bkglist is specified (non-empty)
    # then this function will stack the backgrounds and compare against histlist as if 
    # it is data. The imporant bit is that bkglist is a list of lists. The first index
    # of bkglist corresponds to the index in histlist (the corresponding data). 
    # For example you could have:
    #   histlist = [data1, data2]
    #   bkglist = [[bkg1_1,bkg2_1],[bkg1_2,bkg2_2]]

    if len(histlist) == 1:
        width = 800
        height = 700
        padx = 1
        pady = 1
    elif len(histlist) == 2:
        width = 1200
        height = 700
        padx = 2
        pady = 1
    elif len(histlist) == 3:
        width = 1600
        height = 700
        padx = 3
        pady = 1
    elif len(histlist) == 4:
        width = 1200
        height = 1000
        padx = 2
        pady = 2
    elif len(histlist) == 6 or len(histlist) == 5:
        width = 1600
        height = 1000
        padx = 3
        pady = 2
    else:
        print 'histlist of size ' + len(histlist) + ' not currently supported'
        return 0

    myCan = TCanvas(name,name,width,height)
    myCan.Divide(padx,pady)

    # Just some colors that I think work well together and a bunch of empty lists for storage if needed
    colors = [kRed,kMagenta,kGreen,kCyan]
    stacks = []
    tot_hists = []
    legends = []
    mains = []
    subs = []
    pulls = []
    logString = ''

    # For each hist/data distribution
    for hist_index, hist in enumerate(histlist):
        # Grab the pad we want to draw in
        myCan.cd(hist_index+1)
        if len(histlist) > 1:
            thisPad = myCan.GetPrimitive(name+'_'+str(hist_index+1))
            thisPad.cd()
        
        

        # If this is a TH2, just draw the lego
        if hist.ClassName().find('TH2') != -1:
            if logy == True:
                gPad.SetLogy()
                logString = '_semilog'
            hist.Draw('lego')
            if len(bkglist) > 0:
                print 'ERROR: It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.'
        
        # Otherwise it's a TH1 hopefully
        else:
            hist.SetLineColor(kBlack)
            hist.SetMarkerStyle(8)
            
            # If there are no backgrounds, only plot the data (semilog if desired)
            if len(bkglist) == 0:
                hist.Draw('p e')
            
            # Otherwise...
            else:
                # Create some subpads, a legend, a stack, and a total bkg hist that we'll use for the error bars
                mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.3, 1, 1))
                subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 1, 0.3))
                legends.append(TLegend(0.7,0.8,0.95,0.93))
                stacks.append(THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
                tot_hists.append(TH1F(hist.GetName()+'_tot',hist.GetName()+'_tot',hist.GetNbinsX(),hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax()))

                # Set margins and make these two pads primitives of the division, thisPad
                mains[hist_index].SetBottomMargin(0.0)

                mains[hist_index].SetLeftMargin(0.16)
                mains[hist_index].SetRightMargin(0.05)
                mains[hist_index].SetTopMargin(0.1)

                subs[hist_index].SetLeftMargin(0.16)
                subs[hist_index].SetRightMargin(0.05)
                subs[hist_index].SetTopMargin(0)
                subs[hist_index].SetBottomMargin(0.3)
                mains[hist_index].Draw()
                subs[hist_index].Draw()

                # Build the stack
                for bkg_index,bkg in enumerate(bkglist[hist_index]):     # Won't loop if bkglist is empty
                    bkg.Sumw2()
                    tot_hists[hist_index].Add(bkg)
                    bkg.SetLineColor(kBlack)

                    if bkg.GetName().find('qcd') != -1:
                        bkg.SetFillColor(kYellow)

                    else:
                        bkg.SetFillColor(colors[bkg_index])

                    stacks[hist_index].Add(bkg)

                    legends[hist_index].AddEntry(bkg,bkg.GetName()[:bkg.GetName().find('__')],'f')
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()
                if logy == True:
                    mains[hist_index].SetLogy()
                    logString = '_semilog'

                # Set y max of all hists to be the same to accomodate the tallest
                histList = [stacks[hist_index],tot_hists[hist_index],hist]

                yMax = histList[0].GetMaximum()
                maxHist = histList[0]
                for h in range(1,len(histList)):
                    if histList[h].GetMaximum() > yMax:
                        yMax = histList[h].GetMaximum()
                        maxHist = histList[h]
                for h in histList:
                    h.SetMaximum(yMax*1.1)

                # Now draw the main pad
                stacks[hist_index].Draw('hist')

                tot_hists[hist_index].SetFillColor(kBlack)
                tot_hists[hist_index].SetFillStyle(3354)

                tot_hists[hist_index].Draw('e2 same')
                legends[hist_index].Draw()

                hist.Draw('p e same')

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull
                pulls.append(Make_Pull_plot(hist,tot_hists[hist_index]))
                pulls[hist_index].SetFillColor(kBlue)
                pulls[hist_index].SetTitle(";"+hist.GetXaxis().GetTitle()+";(Data-Bkg)/#sigma")
                pulls[hist_index].SetStats(0)

                LS = .13

                pulls[hist_index].GetYaxis().SetRangeUser(-2.9,2.9)
                pulls[hist_index].GetYaxis().SetTitleOffset(0.4)
                pulls[hist_index].GetXaxis().SetTitleOffset(0.9)
                             
                pulls[hist_index].GetYaxis().SetLabelSize(LS)
                pulls[hist_index].GetYaxis().SetTitleSize(LS)
                pulls[hist_index].GetYaxis().SetNdivisions(306)
                pulls[hist_index].GetXaxis().SetLabelSize(LS)
                pulls[hist_index].GetXaxis().SetTitleSize(LS)
                pulls[hist_index].Draw('hist')

    myCan.Print(tag+'/plots/'+name+logString+'.pdf','pdf')

# Not yet integrated/used
def Make_Pull_plot( DATA,BKG):
    BKGUP, BKGDOWN = Make_up_down(BKG)
    pull = DATA.Clone(DATA.GetName()+"_pull")
    pull.Add(BKG,-1)
    sigma = 0.0
    FScont = 0.0
    BKGcont = 0.0
    for ibin in range(1,pull.GetNbinsX()+1):
        FScont = DATA.GetBinContent(ibin)
        BKGcont = BKG.GetBinContent(ibin)
        if FScont>=BKGcont:
            FSerr = DATA.GetBinErrorLow(ibin)
            BKGerr = abs(BKGUP.GetBinContent(ibin)-BKG.GetBinContent(ibin))
        if FScont<BKGcont:
            FSerr = DATA.GetBinErrorUp(ibin)
            BKGerr = abs(BKGDOWN.GetBinContent(ibin)-BKG.GetBinContent(ibin))
        sigma = sqrt(FSerr*FSerr + BKGerr*BKGerr)
        if FScont == 0.0:
            pull.SetBinContent(ibin, 0.0 )  
        else:
            if sigma != 0 :
                pullcont = (pull.GetBinContent(ibin))/sigma
                pull.SetBinContent(ibin, pullcont)
            else :
                pull.SetBinContent(ibin, 0.0 )
    return pull


def Make_up_down(hist):
    hist_up = hist.Clone(hist.GetName()+'_up')
    hist_down = hist.Clone(hist.GetName()+'_down')

    for xbin in range(1,hist.GetNbinsX()+1):
        errup = hist.GetBinErrorUp(xbin)
        errdown = hist.GetBinErrorLow(xbin)
        nom = hist.GetBinContent(xbin)

        hist_up.SetBinContent(xbin,nom+errup)
        hist_down.SetBinContent(xbin,nom-errdown)

    return hist_up,hist_down

