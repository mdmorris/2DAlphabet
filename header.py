import ROOT
from ROOT import *
import os
import math
from math import sqrt


# Function stolen from https://stackoverflow.com/questions/9590382/forcing-python-json-module-to-work-with-ascii
def ascii_encode_dict(data):    
    ascii_encode = lambda x: x.encode('ascii') if isinstance(x, unicode) else x 
    return dict(map(ascii_encode, pair) for pair in data.items())

def getRRVs(inputConfig):#, blinded):
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

    # if blinded:
    #     xSigStart = inputConfig['BINNING']['X']['SIGSTART']
    #     xSigEnd = inputConfig['BINNING']['X']['SIGEND']
    #     xLowRRV = RooRealVar(xlowname,xlowname,xlow,xSigStart)
    #     xHighRRV = RooRealVar(xhighname,xhighname,xSigEnd,xhigh)

    #     return xRRV,xLowRRV,xHighRRV,yRRV
    
    # else:
    return xRRV,yRRV

def copyHistWithNewXbounds(thisHist,copyName,newBinWidthX,xNewBinsLow,xNewBinsHigh):
    # Make a copy with the same Y bins but new X bins
    nBinsY = thisHist.GetNbinsY()
    yBinsLow = thisHist.GetYaxis().GetXmin()
    yBinsHigh = thisHist.GetYaxis().GetXmax()
    nNewBinsX = int((xNewBinsHigh-xNewBinsLow)/float(newBinWidthX))
    # Use copyName with _temp to avoid overwriting if thisHist has the same name
    # We can do this at the end but not before we're finished with thisHist
    histCopy = TH2F(copyName+'_temp',copyName+'_temp',nNewBinsX,xNewBinsLow,xNewBinsHigh,nBinsY,yBinsLow,yBinsHigh)
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

    # Will now set the copyName which will overwrite thisHist if it has the same name
    histCopy.SetName(copyName)
    histCopy.SetTitle(copyName)

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

def rebinY(thisHist,name,tag,new_y_bins_array):
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

    makeCan(name+'_rebin_compare',tag,[rebinned,thisHist])
    return rebinned

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



def makeCan(name, tag, histlist, bkglist=[],signals=[],colors=[],titles=[],logy=False,rootfile=False,xtitle='',ytitle='',dataOff=False):  
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
        print 'histlist of size ' + str(len(histlist)) + ' not currently supported'
        return 0

    myCan = TCanvas(name,name,width,height)
    myCan.Divide(padx,pady)

    # Just some colors that I think work well together and a bunch of empty lists for storage if needed
    default_colors = [kRed,kMagenta,kGreen,kCyan,kBlue]
    if len(colors) == 0:   
        colors = default_colors
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
            gPad.SetLeftMargin(0.2)
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)
            hist.GetXaxis().SetTitleOffset(1.5)
            hist.GetYaxis().SetTitleOffset(2.3)
            hist.GetZaxis().SetTitleOffset(1.8)
            if len(titles) > 0:
                hist.SetTitle(titles[hist_index])

            hist.Draw('lego')
            if len(bkglist) > 0:
                print 'ERROR: It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.'
        
        # Otherwise it's a TH1 hopefully
        else:
            alpha = 1
            if dataOff:
                alpha = 0
            hist.SetLineColorAlpha(kBlack,alpha)
            hist.SetMarkerColorAlpha(kBlack,alpha)
            hist.SetMarkerStyle(8)
            
            # If there are no backgrounds, only plot the data (semilog if desired)
            if len(bkglist) == 0:
                hist.GetXaxis().SetTitle(xtitle)
                hist.GetYaxis().SetTitle(ytitle)
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.Draw('p e')
            
            # Otherwise...
            else:
                # Create some subpads, a legend, a stack, and a total bkg hist that we'll use for the error bars
                if not dataOff:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.3, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 1, 0.3))

                else:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.1, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 0, 0))

                legends.append(TLegend(0.65,0.6,0.95,0.93))
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
                    # bkg.Sumw2()
                    tot_hists[hist_index].Add(bkg)
                    bkg.SetLineColor(kBlack)
                    if logy:
                        bkg.SetMinimum(1e-3)

                    if bkg.GetName().find('qcd') != -1:
                        bkg.SetFillColor(kYellow)

                    else:
                        if colors[bkg_index] != None:
                            bkg.SetFillColor(colors[bkg_index])
                        else:
                            bkg.SetFillColor(default_colors[bkg_index])

                    stacks[hist_index].Add(bkg)

                    legends[hist_index].AddEntry(bkg,bkg.GetName().split('_')[0],'f')
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()

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
                    if logy == True:
                        h.SetMaximum(yMax*10)

                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.SetTitleOffset(1.5,"xy")
                hist.GetYaxis().SetTitle('Events')
                if logy == True:
                    hist.SetMinimum(1e-3)
                hist.Draw('pe')

                stacks[hist_index].Draw('same hist')

                # Do the signals
                if len(signals) > 0: 
                    signals[hist_index].SetLineColor(kBlue)
                    signals[hist_index].SetLineWidth(2)
                    if logy == True:
                        signals[hist_index].SetMinimum(1e-3)
                    legends[hist_index].AddEntry(signals[hist_index],signals[hist_index].GetName().split('_')[0],'L')
                    signals[hist_index].Draw('hist same')

                tot_hists[hist_index].SetFillColor(kBlack)
                tot_hists[hist_index].SetFillStyle(3354)

                tot_hists[hist_index].Draw('e2 same')
                legends[hist_index].Draw()

                if not dataOff:
                    legends[hist_index].AddEntry(hist,'data','p')
                    hist.Draw('p e same')

                gPad.RedrawAxis()

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

                pulls[hist_index].GetXaxis().SetTitle(xtitle)
                pulls[hist_index].GetYaxis().SetTitle("(Data-Bkg)/#sigma")
                pulls[hist_index].Draw('hist')

                if logy == True:
                    mains[hist_index].SetLogy()

    if rootfile:
        myCan.Print(tag+'plots/'+name+'.root','root')
    else:
        myCan.Print(tag+'plots/'+name+'.png','png')


def FindCommonString(string_list):
    to_match = ''   # initialize the string we're looking for/building
    for s in string_list[0]:    # for each character in the first string
        passed = True
        for istring in range(1,len(string_list)):   # compare to_match+s against strings in string_list
            string = string_list[istring]
            if to_match not in string:                  # if in the string, add more
                passed = False
            
        if passed == True:
            to_match+=s

    if to_match[-2] == '_':
        return to_match[:-2] 
    else:
        return to_match[:-1]                # if not, return to_match minus final character

    return to_match[:-2]
        

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


def checkFitForm(xFitForm,yFitForm):
    if '@0' in xFitForm:
        print 'ERROR: @0 in XFORM. This is reserved for the variable x. Please either replace it with x or start your parameter naming with @1. Quitting...'
        quit()
    if '@0' in yFitForm:
        print 'ERROR: @0 in YFORM. This is reserved for the variable y. Please either replace it with y or start your parameter naming with @1. Quitting...'
        quit()
    if 'x' in yFitForm and 'exp' not in yFitForm:
        print 'ERROR: x in YFORM. Did you mean to write "y"? Quitting...'
        quit()
    if 'y' in xFitForm:
        print 'ERROR: y in XFORM. Did you mean to write "x"? Quitting...'
        quit()


def RFVform2TF1(RFVform,shift=0):
    TF1form = ''
    lookingForDigit = False
    for index,char in enumerate(RFVform):
        # print str(index) + ' ' + char
        if char == '@':
            TF1form+='['
            lookingForDigit = True
        elif lookingForDigit:
            if char.isdigit():
                if RFVform[index+1].isdigit() == False:     # if this is the last digit
                    TF1form+=str(int(char)+shift)
            else:
                TF1form+=']'+char
                lookingForDigit = False
        else:
            TF1form+=char

    return TF1form


# Taken from Kevin's limit_plot_shape.py
def make_smooth_graph(h2,h3):
    h2 = TGraph(h2)
    h3 = TGraph(h3)
    npoints = h3.GetN()
    h3.Set(2*npoints+2)
    for b in range(npoints+2):
        x1, y1 = (ROOT.Double(), ROOT.Double())
        if b == 0:
            h3.GetPoint(npoints-1, x1, y1)
        elif b == 1:
            h2.GetPoint(npoints-b, x1, y1)
        else:
            h2.GetPoint(npoints-b+1, x1, y1)
        h3.SetPoint(npoints+b, x1, y1)
    return h3

# Built to wait for condor jobs to finish and then check that they didn't fail
# The script that calls this function will quit if there are any job failures
# listOfJobs input should be whatever comes before '.listOfJobs' for the set of jobs you submitted
def WaitForJobs( listOfJobs ):
    # Runs grep to count the number of jobs - output will have non-digit characters b/c of wc
    preNumberOfJobs = subprocess.check_output('grep "python" '+listOfJobs+' | wc -l', shell=True)
    commentedNumberOfJobs = subprocess.check_output('grep "# python" '+listOfJobs+'.listOfJobs | wc -l', shell=True)

    # Get rid of non-digits and convert to an int
    preNumberOfJobs = int(filter(lambda x: x.isdigit(), preNumberOfJobs))
    commentedNumberOfJobs = int(filter(lambda x: x.isdigit(), commentedNumberOfJobs))
    numberOfJobs = preNumberOfJobs - commentedNumberOfJobs

    finishedJobs = 0
    # Rudementary progress bar
    while finishedJobs < numberOfJobs:
        # Count how many output files there are to see how many jobs finished
        # the `2> null.txt` writes the stderr to null.txt instead of printing it which means
        # you don't have to look at `ls: output_*.log: No such file or directory`
        finishedJobs = subprocess.check_output('ls output_*.log 2> null.txt | wc -l', shell=True)
        finishedJobs = int(filter(lambda x: x.isdigit(), finishedJobs))
        sys.stdout.write('\rProcessing ' + str(listOfJobs) + ' - ')
        # Print the count out as a 'progress bar' that refreshes (via \r)
        sys.stdout.write("%i / %i of jobs finished..." % (finishedJobs,numberOfJobs))
        # Clear the buffer
        sys.stdout.flush()
        # Sleep for one second
        time.sleep(1)


    print 'Jobs completed. Checking for errors...'
    numberOfTracebacks = subprocess.check_output('grep -i "Traceback" output*.log | wc -l', shell=True)
    numberOfSyntax = subprocess.check_output('grep -i "Syntax" output*.log | wc -l', shell=True)

    numberOfTracebacks = int(filter(lambda x: x.isdigit(), numberOfTracebacks))
    numberOfSyntax = int(filter(lambda x: x.isdigit(), numberOfSyntax))

    # Check there are no syntax or traceback errors
    # Future idea - check output file sizes
    if numberOfTracebacks > 0:
        print str(numberOfTracebacks) + ' job(s) failed with traceback error'
        quit()
    elif numberOfSyntax > 0:
        print str(numberOfSyntax) + ' job(s) failed with syntax error'
        quit()
    else:
        print 'No errors!'

def Inter(g1,g2):
    xaxisrange = g1.GetXaxis().GetXmax()-g1.GetXaxis().GetXmin()
    xaxismin = g1.GetXaxis().GetXmin()
    inters = []
    for x in range(0,10000):
        xpoint = xaxismin + (float(x)/1000.0)*xaxisrange
        xpoint1 = xaxismin + (float(x+1)/1000.0)*xaxisrange
        Pr1 = g1.Eval(xpoint)
        Pr2 = g2.Eval(xpoint)
        Po1 = g1.Eval(xpoint1)
        Po2 = g2.Eval(xpoint1)
        if (Pr1-Pr2)*(Po1-Po2)<0:
            inters.append(0.5*(xpoint+xpoint1))
        
    return inters

def applyFitMorph(process, region ,has_shape_uncert, inputConfig, dists, w, suffix=''):
    # This function will apply the findings of the blinded fit to the signal region 
    # This is necessary because Combine knows NOTHING about this region

    # dists - dictionary with distributions (nominal, syst_up, syst_down,...) for this process and region
    # w - workspace output from combine fit

    # Hack
    allVars = []

    # Setup the axes
    xvar,yvar = getRRVs(inputConfig)
    allVars.extend([xvar,yvar])
    x_low = inputConfig['BINNING']['X']['LOW']
    x_high = inputConfig['BINNING']['X']['HIGH']
    x_nbins = inputConfig['BINNING']['X']['NBINS']
    x_binWidth = float(x_high-x_low)/float(x_nbins)
    sigstart = inputConfig['BINNING']['X']['SIGSTART']
    sigend = inputConfig['BINNING']['X']['SIGEND']
    y_low = inputConfig['BINNING']['Y']['LOW']
    y_high = inputConfig['BINNING']['Y']['HIGH']
    y_nbins = inputConfig['BINNING']['Y']['NBINS']

    normlist = RooArgList() # Stores ALL norms for RooProduct at the end
    allVars.append(normlist)


    # Morph only if we had shape based systematics
    if has_shape_uncert:
        # Storage for the pre-fit signal region for nominal and each systematic shape
        signal_hists = {}
        for dist in dists.keys():
            if dist.find('_unblinded') != -1:
                newkey = dist[:dist.find('_unblinded')]
                if dist[:dist.find('_unblinded')] != 'nominal':
                    hist_name = process+'_'+region+suffix+'_'+newkey+'_presignal'
                else:
                    hist_name = process+'_'+region+suffix+'_presignal'
                signal_hists[newkey] = copyHistWithNewXbounds(dists[dist],hist_name,x_binWidth,sigstart,sigend)

        # Build the pieces for the FVIHP2D2 along with the AsymPow (normalization) pieces
        # Get TList with pdfs and RooArgList with RooRealVars for FVIHP2D2 construction
        coefList = RooArgList()                                     # RAL is ordered [a,b,c,...]
        pdfList = TList()                                           # TList is ordered [nom, a_hi, a_lo, b_hi, b_lo,...]
        pdfList.Add(signal_hists['nominal'])
        asympows = {}
        for n in inputConfig['PROCESS'][process]['SYSTEMATICS']:
            if inputConfig['SYSTEMATIC'][n]['CODE'] > 1: # 2 or 3
                thisVar = w.var(n)
                if not thisVar:
                    continue
                allVars.append(thisVar)
                # FVIHP2D2 pieces
                coefList.add(thisVar)
                allVars.append(thisVar)
                thisTH2up = signal_hists[n+'Up']
                thisTH2down = signal_hists[n+'Down']
                pdfList.Add(thisTH2up)
                pdfList.Add(thisTH2down)
                # AsymPow pieces
                # for su in inputConfig['PROCESS'][process]['SYSTEMATICS']:
                kappaLowVal = thisTH2down.Integral()/signal_hists['nominal'].Integral()
                kappaHighVal = thisTH2up.Integral()/signal_hists['nominal'].Integral()
                if abs(kappaHighVal-1) > 1e-3 and abs(kappaLowVal-1) > 1e-3:    # Combine won't count the systematic unless this condition is satisfied so we shouldn't either
                    kappaLow = RooConstVar(process+'_'+region+suffix+'_'+n+'_kappaLow',process+'_'+region+suffix+'_'+n+'_kappaLow',kappaLowVal)
                    kappaHigh = RooConstVar(process+'_'+region+suffix+'_'+n+'_kappaHigh',process+'_'+region+suffix+'_'+n+'_kappaHigh',kappaHighVal)
                    allVars.extend([kappaLow,kappaHigh])
                    asympows[n] = AsymPow('systeff_'+region+suffix+'_'+process+'_'+n, 'systeff_'+region+suffix+'_'+process+'_'+n, kappaLow, kappaHigh, thisVar)
                    allVars.extend([kappaLow,kappaHigh])
                    normlist.add(asympows[n])

        # Conditional should be FALSE (I checked this - I think it normalizes each x slice to 1 independently)
        # Last two arguments are the same values used by Combine
        shape_FVIHP = FastVerticalInterpHistPdf2D2(process+'_'+region+suffix+'_FVIHP', process+'_'+region+suffix+'_FVIHP', xvar, yvar, False, pdfList, coefList, 1, 0)
        shape_TH2 = shape_FVIHP.createHistogram(process+'_'+region+suffix+'_signal',xvar,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(yvar,RooFit.Binning(y_nbins,y_low,y_high)))


    # Otherwise just grab the nominal histogram
    else:
        shape_TH2 = copyHistWithNewXbounds(dists['nominal_unblinded'],process+'_'+region+suffix+'_presignal',x_binWidth,sigstart,sigend)


    # Next the lnN normalization contribution can just be grabbed from w
    try:
        thislnNNorm = w.function('n_exp_bin'+region+suffix+'_proc_'+process)
        allVars.append(thislnNNorm)
        normlist.add(thislnNNorm)
    except:
        print 'MESSAGE: n_exp_bin'+region+suffix+'_proc_'+process + ' does not exist. Assuming no lnN uncertainties for this.'

    full_norm_name = process+'_'+region+suffix+'_fullNorm'
    full_norm = RooProduct(full_norm_name, full_norm_name, normlist)
    allVars.append(full_norm)

    shape_TH2.Scale(abs(full_norm.getValV()))
    return shape_TH2

# Right after I wrote this I realized it's obsolete... It's cool parentheses parsing though so I'm keeping it
def separateXandYfromFitString(fitForm):
    # Look for all opening and closing parentheses
    openIndexes = []
    closedIndexes = []
    for index,char in enumerate(fitForm):
        if char == '(': # start looking for ")" after this "("
            openIndexes.append(index)
        if char == ')':
            closedIndexes.append(index)


    # Now pair them by looking at the first in closedIndexes and finding the closest opening one (to the left)
    # Remove the pair from the index lists and repeat
    innerContent = []
    for iclose in closedIndexes:
        diff = len(fitForm)     # max length to start because we want to minimize
        for iopen in openIndexes:
            if iclose > iopen:
                this_diff = iclose - iopen
                if this_diff < diff:
                    diff = this_diff
                    candidateOpen = iopen
        openIndexes.remove(candidateOpen)
        innerContent.append(fitForm[iclose-diff+1:iclose])


    outerContent = []
    for c in innerContent:
        keep_c = True
        for d in innerContent:
            if '('+c+')' in d and c != d:
                keep_c = False
                break
        if keep_c:
            outerContent.append(c)

    if len(outerContent) != 2:
        print 'ERROR: Form of the fit did not factorize correctly. Please make sure it is in (x part)(y part) form. Quitting...'
        quit()
    else:
        for c in outerContent:
            if 'x' in c and 'y' not in c:
                xPart = c
            elif 'x' not in c and 'y' in c:
                yPart = c
            else:
                print 'ERROR: Form of the fit did not factorize correctly. Please make sure it is in (x part)(y part) form. Quitting...'
                quit()

    return xPart,yPart
