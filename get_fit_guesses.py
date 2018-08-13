#####################################################################################################
# get_fit_guesses.py - written by Lucas Corcodilos, 5/10/18                                         #
# --------------------------------------------------------                                          #
# This script takes 2D distributions as input, derives a binning scheme for the y-axis that tries   #
# to have similar statistics per bin, fits the x-axis in each y-bin, and then fits the x-axis fit   #
# coefficients as a function of the y-axis. This is also known as psuedo-2D Alphabet since it does  #
# the 2D fit in slices.                                                                             #
#####################################################################################################

import ROOT
from ROOT import *
import array
from array import array
import math
from math import sqrt
import pprint
pp = pprint.PrettyPrinter(indent = 2)


import header
from header import *


def main(inputConfig,blinded,tag,nslices=6,sigma=25):

    # Grab everything
    file = TFile.Open(inputConfig['PROCESS']['data_obs']['FILE'])
    
    thisPass = file.Get(inputConfig['PROCESS']['data_obs']['HISTPASS']).Clone('qcd_pass')
    thisFail = file.Get(inputConfig['PROCESS']['data_obs']['HISTFAIL']).Clone('qcd_fail')
    thisTotal = thisPass.Clone('thisTotal')     # Reconstruct pre-pass/fail amount for statistics reasons later
    thisTotal.Add(thisFail)

    xmax = inputConfig['BINNING']['X']['HIGH']
    xmin = inputConfig['BINNING']['X']['LOW']
    xnbins = inputConfig['BINNING']['X']['NBINS']
    xBinWidth = float(xmax - xmin)/float(xnbins)
    try:
        x_title = inputConfig['BINNING']['X']['TITLE']
    except:
        x_title = ''

    ymin = inputConfig['BINNING']['Y']['LOW']
    ymax = inputConfig['BINNING']['Y']['HIGH']
    try:
        y_title = inputConfig['BINNING']['Y']['TITLE']
    except:
        y_title = ''

    # Subtract away non-qcd processes
    for nonqcd in [process for process in inputConfig['PROCESS'].keys() if process != 'HELP' and process != 'data_obs']:
        nonqcd_file = TFile.Open(inputConfig['PROCESS'][nonqcd]['FILE'])
        nonqcd_pass = nonqcd_file.Get(inputConfig['PROCESS'][nonqcd]['HISTPASS'])
        nonqcd_fail = nonqcd_file.Get(inputConfig['PROCESS'][nonqcd]['HISTFAIL'])
        thisPass.Add(nonqcd_pass,-1)
        thisFail.Add(nonqcd_fail,-1)


    # Zero any negative bins
    for this in [thisPass,thisFail]:
        for xbin in range(1,thisPass.GetNbinsX()+1):
            for ybin in range(1,thisPass.GetNbinsY()+1):
                if this.GetBinContent(xbin,ybin) < 0:
                    this.SetBinContent(xbin,ybin,0)
                if this.GetBinContent(xbin,ybin) < 0:
                    this.SetBinContent(xbin,ybin,0)


    # Need to figure out how to slice up the y-axis based on the total statistics
    # (pass+fail) for each 2D bin

    # Get an average number of events per slice
    total_events = thisTotal.Integral()
    events_per_slice = total_events/float(nslices)

    #############################
    #   Derive the new y bins   #
    #############################

    # Loop over ybins, get the number of events in that ybin, and check if less than events_per_slice
    ysum = 0            # ysum for after we've confirmed that the next bin does not lead to a number greater than events_per_slice
    ysum_temp = 0       # ysum allowed to overflow
    new_y_bins = [float(inputConfig['BINNING']['Y']['LOW'])]
    for ybin in range(1,thisTotal.GetNbinsY()+1):                 # For each ybin
        # print 'ybin = ' + str(ybin) + ', ' + str(ysum_temp)

        # If on last ybin, add the final edge to the list
        # - - - - - - - - - - - - - - - - - - - - - - - - 
        # There's an interesting note to make here about the number of slices that one can do. There are three important facts here:
        # 1) this method only reduces the number of bins
        # 2) every bin edge in the new binning was a bin edge in the previous binning
        # 3) new bins have fewer than the target events_per_slice UNLESS the new bin is only made up of one bin in which case it has more than events_per_slice
        #
        # This means two scenarios happen
        # Example 1) if old bins 1-3 have fewer events than events_per_slice but old bins 1-4 have more than events_per_slice,
        #            then new bin 1 will contain the contents of old bins 1-3 (no 4). 
        # Example 2) But if the content of bin 4 is GREATER THAN events_per_slice,
        #            then new bin 2 will contain the contents of ONLY old bin 4. 
        #
        # Here is the interesting part. If the Example 2 scenario happens enough times (as is possible with a peaked distribution), then 
        # there won't be enough events to "go around" for the tail of the distribution. So you might ask for 9 slices but because 
        # slices 2-4 have many more events in them than 3*events_per_slice, then the next slices will have to count more events in
        # the tail and you'll "use" all of the events in the tail before you get to the 9th slice. In other words, this method naturally caps itself at a certain number of slices.
        if ybin == thisTotal.GetNbinsY() or len(new_y_bins) == nslices:
            new_y_bins.append(inputConfig['BINNING']['Y']['HIGH'])  # Add the final edge
            break                                                   # This final bin will most likely have more than events_per_slice in it

        # Otherwise, if we're still building the slice list...
        else:
            for xbin in range(1,thisTotal.GetNbinsX()):             # For each xbin
                ysum_temp += thisTotal.GetBinContent(xbin,ybin)     # Add the 2D bin content to the temp sum for the ybin

            # If less, set ysum and go onto the next bin
            if ysum_temp < events_per_slice:
                ysum = ysum_temp
                continue
            # Otherwise, cut off the slice at the previous ybin and restart the sum with this ybin
            else:
                # print 'Done with this bin total ' + str(ysum)
                ysum_temp = 0
                for xbin in range(1,thisTotal.GetNbinsX()):
                    ysum_temp += thisTotal.GetBinContent(xbin,ybin)
                new_y_bins.append(thisTotal.GetYaxis().GetBinLowEdge(ybin))
                # print new_y_bins

    print 'Will bin y-axis for fit guesses using bins ',
    print new_y_bins

    # new_y_bins_array = array('d',new_y_bins)
    new_y_bins_array = array('d',[800.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3200.0, 4000.0])    

    ##############################
    #   Rebin the distributions  #
    ##############################

    # Blind if necessary and calculate Rpf in 2D
    if blinded:
        sigstart = inputConfig['BINNING']['X']['SIGSTART']
        sigend = inputConfig['BINNING']['X']['SIGEND']

        # X rebin and blind
        rebinnedXPass = copyHistWithNewXbounds(thisPass,'rebinnedXPass',xBinWidth,xmin,xmax)
        rebinnedXFail = copyHistWithNewXbounds(thisFail,'rebinnedXFail',xBinWidth,xmin,xmax)

        lowPass = copyHistWithNewXbounds(rebinnedXPass,'lowPass',xBinWidth,xmin,sigstart)
        lowFail = copyHistWithNewXbounds(rebinnedXFail,'lowFail',xBinWidth,xmin,sigstart)

        highPass = copyHistWithNewXbounds(rebinnedXPass,'highPass',xBinWidth,sigend,xmax)
        highFail = copyHistWithNewXbounds(rebinnedXFail,'highFail',xBinWidth,sigend,xmax)

        blindedPass = makeBlindedHist(rebinnedXPass,lowPass,highPass)
        blindedFail = makeBlindedHist(rebinnedXFail,lowFail,highFail)


        # Rp/f before y rebin
        thisTotal = blindedPass.Clone('thisTotal')  # Need to remake this if we've blinded
        thisTotal.Add(blindedFail)

        Rpf = blindedPass.Clone('Rpf')
        Rpf.Divide(blindedFail)

        # Rp/f after y rebin
        rebinnedPass = rebinY(blindedPass,'rebinnedPass',tag,new_y_bins_array)
        rebinnedFail = rebinY(blindedFail,'rebinnedFail',tag,new_y_bins_array)

        rebinnedRpf = rebinnedPass.Clone('rebinnedRpf')
        rebinnedRpf.Divide(rebinnedFail)


    # Otherwise just get the Rpf
    else:
        # Rp/f before y rebin
        rebinnedXPass = copyHistWithNewXbounds(thisPass,'rebinnedXPass',xBinWidth,xmin,xmax)
        rebinnedXFail = copyHistWithNewXbounds(thisFail,'rebinnedXFail',xBinWidth,xmin,xmax)

        Rpf = rebinnedXPass.Clone('Rpf')
        Rpf.Divide(rebinnedXFail)

        # Rp/f after y rebin
        rebinnedPass = rebinY(rebinnedXPass,'rebinnedPass',tag,new_y_bins_array)
        rebinnedFail = rebinY(rebinnedXFail,'rebinnedFail',tag,new_y_bins_array)

        rebinnedRpf = rebinnedPass.Clone('rebinnedRpf')
        rebinnedRpf.Divide(rebinnedFail)

    # Plot comparisons out
    # makeCan('rebinned_dists',tag,[rebinnedXPass,rebinnedXFail,rebinnedPass,rebinnedFail])
    makeCan('Rpfs',tag,[rebinnedRpf],xtitle=x_title,ytitle=y_title)

    ###############################################
    # Determine fit function from the inputConfig #
    ###############################################
    if 'XFORM' in inputConfig['FIT'].keys() and 'YFORM' in inputConfig['FIT'].keys():
        # Do some quick checks to make sure these are formatted correctly
        checkFitForm(inputConfig['FIT']['XFORM'],inputConfig['FIT']['YFORM'])
        # Determine number of params in each direction
        nxparams = max([int(param[1:]) for param in inputConfig['FIT'].keys() if param.find('X') != -1 and param != 'XFORM'])
        nyparams = max([int(param[1:]) for param in inputConfig['FIT'].keys() if param.find('Y') != -1 and param != 'YFORM'])
        # Get the strings
        xFuncString = RFVform2TF1(inputConfig['FIT']['XFORM'],-1)
        yFuncString = RFVform2TF1(inputConfig['FIT']['YFORM'],-1)
        yFuncString = yFuncString.replace('y','x')

        # Make the fixed form of the formula
        funcString = ''
        paramIndex = 0
        for xparam in range(nxparams):
            for yparam in range(nyparams):
                funcString += '['+str(paramIndex)+']*x**'+str(xparam)+'*y**'+str(yparam)+'+'
                paramIndex += 1
        funcString = funcString[:-1]


    elif 'FORM' in inputConfig['FIT'].keys():
        funcString = inputConfig['FIT']['FORM']
        # Reconstruct x
        xFuncString = ''
        nxparams = 0
        for xparam in inputConfig['FIT'].keys():
            if xparam.find('X') != -1:                                           # For each X*Y0
                if xparam.find('Y0') != -1:
                    powerIndex = xparam[xparam.find('X')+1]
                    xFuncString += '['+str(powerIndex)+']*x**'+str(powerIndex)+'+'
                    nxparams += 1
        xFuncString = xFuncString[:-1]

        # Reconstruct y
        yFuncString = ''
        nyparams = 0
        for yparam in inputConfig['FIT'].keys():
            if yparam.find('X0Y') != -1:                                           # For each X0Y*
                powerIndex = yparam[yparam.find('Y')+1]
                yFuncString += '['+str(powerIndex)+']*x**'+str(powerIndex)+'+'
                nyparams += 1
        # Will have a trailing '+' that needs to be removed
        yFuncString = yFuncString[:-1]


    ##################################
    # Now do the fit in the y slices #
    ##################################
    fitResults = {}

    # Book TGraphs to store the fit results as a function of y bins
    for xparam in range(nxparams):
        fitResults['xparam_'+str(xparam)+'_vs_y'] = TH1F('xparam_'+str(xparam)+'_vs_y','xparam_'+str(xparam)+'_vs_y',len(new_y_bins_array)-1,new_y_bins_array)

    # Project each y-axis bin and fit along x - save out coefficients to booked tgraph
    for ybin in range(1,rebinnedRpf.GetNbinsY()+1):
        fitResults['fitSlice_'+str(ybin)] = TF1('fitSlice_'+str(ybin),xFuncString,xmin,xmax)
        projX = rebinnedRpf.ProjectionX('rebinnedRpf_sliceX_'+str(ybin),ybin,ybin,'e o')
        projX.Draw('p e')
        # binCenter = rebinnedRpf.GetYaxis().GetBinCenter(ybin)
        # binLowEdge = rebinnedRpf.GetYaxis().GetBinLowEdge(ybin)
        # binHighEdge = rebinnedRpf.GetYaxis().GetBinUpEdge(ybin)
        projX.Fit(fitResults['fitSlice_'+str(ybin)],'E')
        projX.Draw('p e')
        for ix in range(nxparams):
            fitResults['xparam_'+str(ix)+'_vs_y'].SetBinContent(ybin,fitResults['fitSlice_'+str(ybin)].GetParameter(ix))
            fitResults['xparam_'+str(ix)+'_vs_y'].SetBinError(ybin,fitResults['fitSlice_'+str(ybin)].GetParError(ix))
 

    ########################################################
    # And now fit these parameters as a function of y axis #
    ########################################################

    # Build fit for each parameter distribution along y-axis
    drawList = []
    for xparam in range(nxparams):
        fitResults['fitParam_'+str(xparam)] = TF1('yFunc_'+str(xparam),yFuncString,ymin,ymax)
        # Do the fit
        fitResults['xparam_'+str(xparam)+'_vs_y'].Fit(fitResults['fitParam_'+str(xparam)],"E")
        drawList.append(fitResults['xparam_'+str(xparam)+'_vs_y'])
        # Get and store parameters found
        for iy in range(fitResults['fitParam_'+str(xparam)].GetNpar()):
            fitResults['X'+str(xparam)+'Y'+str(iy)] = fitResults['fitParam_'+str(xparam)].GetParameter(iy)
            fitResults['X'+str(xparam)+'Y'+str(iy)+'err'] = fitResults['fitParam_'+str(xparam)].GetParError(iy)

    makeCan('xparam_v_y',tag,drawList,xtitle=y_title)

    # Remove old fit values and store new ones in inputConfig
    print 'Resetting fit parameters in input config'
    inputConfig['FIT'] = {'FORM':funcString}

    psuedo2D_Rpf = TF2('psuedo2D_Rpf',funcString,xmin,xmax,ymin,ymax)
    paramIndex = 0
    for ix in range(nxparams):
        for iy in range(nyparams):
            param = 'X'+str(ix)+'Y'+str(iy)
            psuedo2D_Rpf.SetParameter(paramIndex,fitResults[param])
            psuedo2D_Rpf.SetParError(paramIndex,fitResults[param+'err'])
            inputConfig['FIT'][param] = {'NOMINAL':None,'LOW':None,'HIGH':None}
            if fitResults[param] < 0:
                inputConfig['FIT'][param]['NOMINAL'] = fitResults[param]
                inputConfig['FIT'][param]['HIGH'] = fitResults[param]+sigma*fitResults[param+'err']#min(fitResults[param]+sigma*fitResults[param+'err'],0)
                inputConfig['FIT'][param]['LOW'] = fitResults[param]-sigma*fitResults[param+'err']
            elif fitResults[param] > 0:
                inputConfig['FIT'][param]['NOMINAL'] = fitResults[param]
                inputConfig['FIT'][param]['HIGH'] = fitResults[param]+sigma*fitResults[param+'err']
                inputConfig['FIT'][param]['LOW'] = fitResults[param]-sigma*fitResults[param+'err']#max(fitResults[param]-sigma*fitResults[param+'err'],0)

            inputConfig['FIT'][param]['ERROR'] = fitResults[param+'err']

            paramIndex+=1

    # Finally draw the surface
    makeCan('psuedo2D_Rpf',tag,[psuedo2D_Rpf],xtitle=x_title,ytitle=y_title)

    return inputConfig
