#####################################################################################################
# get_fit_guesses.py - written by Lucas Corcodilos, 5/10/18                                         #
# --------------------------------------------------------                                          #
# This script takes 2D distributions as input, derives a binning scheme for the y-axis that tries   #
# to have similar statistics per bin, fits the x-axis in each y-bin, and then fits the x-axis fit   #
# coefficients as a function of the y-axis. This is also known as pseudo-2D Alphabet since it does  #
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


def main(inputConfig,blinded,tag,nslices=0,sigma=5):

    # Grab everything
    infile = TFile.Open(inputConfig['PROCESS']['data_obs']['FILE'])
    
    dataPass = infile.Get(inputConfig['PROCESS']['data_obs']['HISTPASS'])
    dataPass = dataPass.Clone('qcd_pass')
    dataFail = infile.Get(inputConfig['PROCESS']['data_obs']['HISTFAIL'])
    dataFail = dataFail.Clone('qcd_fail')

    # Get BINNING info
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
    ynbins = inputConfig['BINNING']['Y']['NBINS']
    yBinWidth = float(ymax - ymin)/float(ynbins)
    try:
        y_title = inputConfig['BINNING']['Y']['TITLE']
    except:
        y_title = ''

    # Subtract away non-qcd processes
    for nonqcd in [process for process in inputConfig['PROCESS'].keys() if process != 'HELP']:
        if inputConfig['PROCESS'][nonqcd]['CODE'] > 1: # remove data and signal from considerations
            print 'Subtracting ' + nonqcd
            nonqcd_file = TFile.Open(inputConfig['PROCESS'][nonqcd]['FILE'])
            nonqcd_pass = nonqcd_file.Get(inputConfig['PROCESS'][nonqcd]['HISTPASS'])
            nonqcd_fail = nonqcd_file.Get(inputConfig['PROCESS'][nonqcd]['HISTFAIL'])
            dataPass.Add(nonqcd_pass,-1)
            dataFail.Add(nonqcd_fail,-1)


    # Zero any negative bins
    for reghist in [dataPass,dataFail]:
        for xbin in range(1,dataPass.GetNbinsX()+1):
            for ybin in range(1,dataPass.GetNbinsY()+1):
                if reghist.GetBinContent(xbin,ybin) < 0:
                    reghist.SetBinContent(xbin,ybin,0)
                if reghist.GetBinContent(xbin,ybin) < 0:
                    reghist.SetBinContent(xbin,ybin,0)



    # Rebin x and y immediately - now we have binning comparable to full fit
    rebinnedXPass = copyHistWithNewXbounds(dataPass,'rebinnedXPass',xBinWidth,xmin,xmax)
    rebinnedXFail = copyHistWithNewXbounds(dataFail,'rebinnedXFail',xBinWidth,xmin,xmax)

    rebinnedPass = copyHistWithNewYbounds(rebinnedXPass,'rebinnedPass',yBinWidth,ymin,ymax)
    rebinnedFail = copyHistWithNewYbounds(rebinnedXFail,'rebinnedFail',yBinWidth,ymin,ymax)

    # Make a total to see what the statistics are like
    dataTotal = rebinnedPass.Clone('dataTotal')     
    dataTotal.Add(rebinnedFail)

    # Need to figure out how to slice up the y-axis based on the total statistics (pass+fail) for each 2D bin
    
    #############################
    #   Derive the new y bins   #
    #############################
    if nslices > 0:
        # Get an average number of events per slice
        total_events = dataTotal.Integral()
        events_per_slice = total_events/float(nslices)

        # Loop over ybins, get the number of events in that ybin, and check if less than events_per_slice
        ysum = 0            # ysum for after we've confirmed that the next bin does not lead to a number greater than events_per_slice (for printing only)
        ysum_temp = 0       # ysum allowed to overflow (for counting)
        new_y_bins = []
        for ybin in range(1,dataTotal.GetNbinsY()+1):                 # For each ybin
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
            
            if ybin == dataTotal.GetNbinsY() or len(new_y_bins) == nslices: # If final bin or max number of slices reached
                new_y_bins.append(inputConfig['BINNING']['Y']['HIGH'])  # Add the final edge
                break                                                   # This final bin will most likely have more than events_per_slice in it

            # Otherwise, if we're still building the slice list...
            else:
                for xbin in range(1,dataTotal.GetNbinsX()):             # For each xbin
                    ysum_temp += dataTotal.GetBinContent(xbin,ybin)     # Add the 2D bin content to the temp sum for the ybin

                # If less, set ysum and go onto the next bin
                if ysum_temp < events_per_slice:
                    ysum = ysum_temp
                    
                # Otherwise, cut off the slice at the previous ybin and restart the sum with this ybin
                else:
                    # print 'Done with this bin total ' + str(ysum)
                    ysum_temp = 0
                    for xbin in range(1,dataTotal.GetNbinsX()):
                        ysum_temp += dataTotal.GetBinContent(xbin,ybin)
                    new_y_bins.append(dataTotal.GetYaxis().GetBinLowEdge(ybin))
                    # print new_y_bin
        

    else:
        new_y_bins = [yBinWidth*i+ymin for i in range(ynbins+1)]

    # new_y_bins = [1000,1200,1400,1600,1800,2000,4000]

    print 'Will bin y-axis for fit guesses using bins ',
    print new_y_bins

    new_y_bins_array = array('d',new_y_bins)

    ######################################################
    #   Rebin the distributions according to new y bins  #
    ######################################################

    # Blind if necessary
    if blinded:
        sigstart = inputConfig['BINNING']['X']['SIGSTART']
        sigend = inputConfig['BINNING']['X']['SIGEND']

        lowPass = copyHistWithNewXbounds(rebinnedPass,'lowPass',xBinWidth,xmin,sigstart)
        lowFail = copyHistWithNewXbounds(rebinnedFail,'lowFail',xBinWidth,xmin,sigstart)
        
        highPass = copyHistWithNewXbounds(rebinnedPass,'highPass',xBinWidth,sigend,xmax)
        highFail = copyHistWithNewXbounds(rebinnedFail,'highFail',xBinWidth,sigend,xmax)

        blindedPass = makeBlindedHist(rebinnedPass,lowPass,highPass)
        blindedFail = makeBlindedHist(rebinnedFail,lowFail,highFail)

        # Rp/f after y rebin
        if nslices > 0:
            finalPass = rebinY(blindedPass,'finalPass',tag,new_y_bins_array)
            finalFail = rebinY(blindedFail,'finalFail',tag,new_y_bins_array)

        else:
            finalPass = blindedPass
            finalFail = blindedFail

        finalPass.SetName('finalPass')
        finalPass.SetTitle('finalPass')
        finalFail.SetName('finalFail')
        finalFail.SetTitle('finalFail')

    # Otherwise just get the Rpf
    else:
        if nslices > 0:
            finalPass = rebinY(rebinnedPass,'finalPass',tag,new_y_bins_array)
            finalFail = rebinY(rebinnedFail,'finalFail',tag,new_y_bins_array)

        else:
            finalPass = rebinnedPass
            finalFail = rebinnedFail

        finalPass.SetName('finalPass')
        finalPass.SetTitle('finalPass')
        finalFail.SetName('finalFail')
        finalFail.SetTitle('finalFail')


    RpfToRemap = finalPass.Clone('RpfToRemap')
    RpfToRemap.Divide(finalFail)

    Rpf = remapToUnity(RpfToRemap)

    # Plot comparisons out
    # makeCan('rebinned_dists',tag,[rebinnedXPass,rebinnedXFail,rebinnedPass,rebinnedFail])
    makeCan('prefit_rpf_lego',tag,[RpfToRemap],xtitle=x_title,ytitle=y_title)
    makeCan('prefit_rpfunit_lego',tag,[Rpf],xtitle=x_title,ytitle=y_title)
    out_rpf = TFile(tag+'/plots/rpf_comparisons.root','RECREATE')
    RpfToRemap.Write()
    Rpf.Write()
    out_rpf.Close()

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
            if 'X' in xparam:                                           # For each X*Y0
                if 'Y0' in xparam:
                    powerIndex = xparam[xparam.find('X')+1]
                    xFuncString += '['+str(powerIndex)+']*x**'+str(powerIndex)+'+'
                    nxparams += 1
        xFuncString = xFuncString[:-1]

        # Reconstruct y
        yFuncString = ''
        nyparams = 0
        for yparam in inputConfig['FIT'].keys():
            if 'X0Y' in yparam:                                           # For each X0Y*
                powerIndex = yparam[yparam.find('Y')+1]
                yFuncString += '['+str(powerIndex)+']*x**'+str(powerIndex)+'+'
                nyparams += 1
        yFuncString = yFuncString[:-1]


    elif 'GFORM' in inputConfig['FIT'].keys():
        # Need to take a 2D function and freeze it in one dimension for each y slice
        funcString = RFVform2TF1(inputConfig['FIT']['GFORM'],0)
        yFuncString = '[0]' # since all y dependence will be absorbed by the funcString, there should be no y dependence on each of the parameters and thus we should fit each with a constant
        nxparams = max([int(param) for param in inputConfig['FIT'].keys() if param != 'GFORM']) +1

    else:
        print 'Fit form not supported in get_fit_guesses.py. Quitting...'
        quit()

    ##################################
    # Now do the fit in the y slices #
    ##################################
    fitResults = {}

    # Book TGraphs to store the fit results as a function of y bins
    unitYbins = array('d',[Rpf.GetYaxis().GetBinLowEdge(b) for b in range(1,Rpf.GetNbinsY()+1)]+[1])
    print unitYbins
    for xparam in range(nxparams):
        fitResults['xparam_'+str(xparam)+'_vs_y'] = TH1F('xparam_'+str(xparam)+'_vs_y','xparam_'+str(xparam)+'_vs_y',Rpf.GetNbinsY(),unitYbins)

    # Project each y-axis bin and fit along x - save out coefficients to booked tgraph
    projXs = []
    pseudo2D_results = TFile.Open(tag+'/pseudo2D_results.root','RECREATE')
    pseudo2D_results.cd()
    for ybin in range(1,Rpf.GetNbinsY()+1):
        # If doing GFORM, define xFuncString here with y bin center plugged in
        if 'GFORM' in inputConfig['FIT'].keys():
            thisYBinCenter = Rpf.GetYaxis().GetBinCenter(ybin)
            xFuncString = funcString.replace('y',str(thisYBinCenter))

        fitResults['fitSlice_'+str(ybin)] = TF1('fitSlice_'+str(ybin),xFuncString,0,1)

        if 'GFORM' in inputConfig['FIT'].keys():
            for p in range(nxparams):
                fitResults['fitSlice_'+str(ybin)].SetParameter(p,inputConfig['FIT'][str(p)]['NOMINAL'])
                fitResults['fitSlice_'+str(ybin)].SetParLimits(p,inputConfig['FIT'][str(p)]['LOW'],inputConfig['FIT'][str(p)]['HIGH'])

        projX = Rpf.ProjectionX('rebinnedRpf_sliceX_'+str(ybin),ybin,ybin,'e o')
        # projX.Draw('p e')
        projX.Write()
        projX.Fit(fitResults['fitSlice_'+str(ybin)],'EM')
        
        # projX.Draw('p e')
        projX.SetMaximum(1.1)
        projX.SetMinimum(0.0)
        projX.SetTitle('fitSlice_'+str(ybin))

        projXs.append(projX)

        for ix in range(nxparams):
            fitResults['xparam_'+str(ix)+'_vs_y'].SetBinContent(ybin,fitResults['fitSlice_'+str(ybin)].GetParameter(ix))
            fitResults['xparam_'+str(ix)+'_vs_y'].SetBinError(ybin,fitResults['fitSlice_'+str(ybin)].GetParError(ix))
 
    if len(projXs) <= 6:
        makeCan('fitSlices_0-'+str(len(projXs)),tag,projXs,xtitle=x_title)

    else:
        chunkedProjX = [projXs[i:i + 6] for i in xrange(0, len(projXs), 6)]
        for i,chunk in enumerate(chunkedProjX):
            makeCan('fitSlices_'+str(i*6)+'-'+str(6*(i+1)),tag,chunk,xtitle=x_title)

    pseudo2D_results.Close()

    ########################################################
    # And now fit these parameters as a function of y axis #
    ########################################################

    # Build fit for each parameter distribution along y-axis
    drawList = []
    for xparam in range(nxparams):
        fitResults['fitParam_'+str(xparam)] = TF1('yFunc_'+str(xparam),yFuncString,0,1)
        # Do the fit
        fitResults['xparam_'+str(xparam)+'_vs_y'].Fit(fitResults['fitParam_'+str(xparam)],"EM")
        drawList.append(fitResults['xparam_'+str(xparam)+'_vs_y'])
        # Get and store parameters found
        for iy in range(fitResults['fitParam_'+str(xparam)].GetNpar()):
            fitResults['X'+str(xparam)+'Y'+str(iy)] = fitResults['fitParam_'+str(xparam)].GetParameter(iy)
            fitResults['X'+str(xparam)+'Y'+str(iy)+'err'] = fitResults['fitParam_'+str(xparam)].GetParError(iy)

    if len(drawList) <= 6:
        makeCan('xparam_v_y',tag,drawList,xtitle=y_title)
    else:
        chunkedDrawList = [drawList[i:i + 6] for i in xrange(0, len(drawList), 6)]
        for i,chunk in enumerate(chunkedDrawList):
            makeCan('xparam_v_y_'+str(i),tag,chunk,xtitle=y_title)


    # Remove old fit values and store new ones in inputConfig
    print 'Resetting fit parameters in input config'

    pseudo2D_Rpf = TF2('pseudo2D_Rpf',funcString,0,1,0,1)
    paramIndex = 0

    if 'GFORM' in inputConfig['FIT'].keys():
        inputConfig['FIT'] = {'GFORM':funcString}
        for p in range(nxparams):
            param = 'X'+str(p)+'Y0'
            pseudo2D_Rpf.SetParameter(paramIndex,fitResults[param])
            pseudo2D_Rpf.SetParError(paramIndex,fitResults[param+'err'])
            inputConfig['FIT'][str(p)] = {'NOMINAL':None,'LOW':None,'HIGH':None}
            inputConfig['FIT'][str(p)]['NOMINAL'] = fitResults[param]
            inputConfig['FIT'][str(p)]['HIGH'] = fitResults[param]+sigma*fitResults[param+'err']
            inputConfig['FIT'][str(p)]['LOW'] = fitResults[param]-sigma*fitResults[param+'err']
            inputConfig['FIT'][str(p)]['ERROR'] = fitResults[param+'err']

            paramIndex+=1

    else:
        inputConfig['FIT'] = {'FORM':funcString}
        for ix in range(nxparams):
            for iy in range(nyparams):
                param = 'X'+str(ix)+'Y'+str(iy)
                pseudo2D_Rpf.SetParameter(paramIndex,fitResults[param])
                pseudo2D_Rpf.SetParError(paramIndex,fitResults[param+'err'])
                inputConfig['FIT'][param] = {'NOMINAL':None,'LOW':None,'HIGH':None}
                # if fitResults[param] < 0:
                inputConfig['FIT'][param]['NOMINAL'] = fitResults[param]
                # if fitResults[param]+sigma*fitResults[param+'err'] < 1.0:
                #     inputConfig['FIT'][param]['HIGH'] = 1.0
                # else:
                inputConfig['FIT'][param]['HIGH'] = fitResults[param]+sigma*fitResults[param+'err']

                # if fitResults[param]-sigma*fitResults[param+'err'] > -1.0:
                #     inputConfig['FIT'][param]['LOW'] = -1.0
                # else:
                inputConfig['FIT'][param]['LOW'] = fitResults[param]-sigma*fitResults[param+'err']
                
                inputConfig['FIT'][param]['ERROR'] = fitResults[param+'err']

                paramIndex+=1

    # Finally draw the surface
    makeCan('pseudo2D_Rpf',tag,[pseudo2D_Rpf],xtitle=x_title,ytitle=y_title)
    # quit()
    return inputConfig
