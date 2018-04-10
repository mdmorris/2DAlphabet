import ROOT
from ROOT import *
import os


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
    
    histCopy.GetXaxis().SetName(thisHist.GetXaxis().GetName())
    histCopy.GetYaxis().SetName(thisHist.GetYaxis().GetName())


    # Loop through the old bins
    for binY in range(1,nBinsY+1):
        # print 'Bin y: ' + str(binY)
        for newBinX in range(1,nNewBinsX+1):
            newBinContent = 0
            newBinXlow = histCopy.GetXaxis().GetBinLowEdge(newBinX)
            newBinXhigh = histCopy.GetXaxis().GetBinUpEdge(newBinX)

            # print '\t New bin x: ' + str(newBinX) + ', ' + str(newBinXlow) + ', ' + str(newBinXhigh)
            for oldBinX in range(1,thisHist.GetNbinsX()+1):
                if thisHist.GetXaxis().GetBinLowEdge(oldBinX) >= newBinXlow and thisHist.GetXaxis().GetBinUpEdge(oldBinX) <= newBinXhigh:
                    # print '\t \t Old bin x: ' + str(oldBinX) + ', ' + str(thisHist.GetXaxis().GetBinLowEdge(oldBinX)) + ', ' + str(thisHist.GetXaxis().GetBinUpEdge(oldBinX))
                    # print '\t \t Adding content ' + str(thisHist.GetBinContent(oldBinX,binY))
                    newBinContent += thisHist.GetBinContent(oldBinX,binY)

            # print '\t Setting content ' + str(newBinContent)
            histCopy.SetBinContent(newBinX,binY,newBinContent)

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


def make4x4Can(name,h1,h2,h3,h4):
    my4x4Can = TCanvas(name,name,1200,1000)
    my4x4Can.Divide(2,2)

    my4x4Can.cd(1)
    h1.Draw('lego')
    my4x4Can.cd(2)
    h2.Draw('lego')
    my4x4Can.cd(3)
    h3.Draw('lego')
    my4x4Can.cd(4)
    h4.Draw('lego')

    my4x4Can.Print('plots/'+name+'.pdf','pdf')

def make1x1Can(name,h1):
    my1x1Can = TCanvas(name,name,800,700)
    my1x1Can.cd()
    h1.Draw('lego')
    my1x1Can.Print('plots/'+name+'.pdf','pdf')




