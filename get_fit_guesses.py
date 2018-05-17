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


def main(inputConfig,blinded,tag,nslices=8,sigma=1):

    # Grab everything
    file = TFile.Open(inputConfig['PROCESS']['data_obs']['FILE'])
    
    thisPass = file.Get(inputConfig['PROCESS']['data_obs']['HISTPASS'])
    thisFail = file.Get(inputConfig['PROCESS']['data_obs']['HISTFAIL'])
    thisTotal = thisPass.Clone('thisTotal')     # Reconstruct pre-pass/fail amount for statistics reasons later
    thisTotal.Add(thisFail)

    xmax = inputConfig['BINNING']['X']['HIGH']
    xmin = inputConfig['BINNING']['X']['LOW']
    xnbins = inputConfig['BINNING']['X']['NBINS']
    xBinWidth = float(xmax - xmin)/float(xnbins)

    ymin = inputConfig['BINNING']['Y']['LOW']
    ymax = inputConfig['BINNING']['Y']['HIGH']

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


    


    # Blind if necessary and calculate Rpf in 2D
    if blinded:
        sigstart = inputConfig['BINNING']['X']['SIGSTART']
        sigend = inputConfig['BINNING']['X']['SIGEND']

        rebinnedPass = copyHistWithNewXbounds(thisPass,'rebinnedPass',xBinWidth,xmin,xmax)
        rebinnedFail = copyHistWithNewXbounds(thisFail,'rebinnedFail',xBinWidth,xmin,xmax)

        lowPass = copyHistWithNewXbounds(rebinnedPass,'lowPass',xBinWidth,xmin,sigstart)
        lowFail = copyHistWithNewXbounds(rebinnedFail,'lowFail',xBinWidth,xmin,sigstart)

        highPass = copyHistWithNewXbounds(rebinnedPass,'highPass',xBinWidth,sigend,xmax)
        highFail = copyHistWithNewXbounds(rebinnedFail,'highFail',xBinWidth,sigend,xmax)

        blindedPass = makeBlindedHist(rebinnedPass,lowPass,highPass)
        blindedFail = makeBlindedHist(rebinnedFail,lowFail,highFail)

        thisTotal = blindedPass.Clone('thisTotal')  # Need to remake this if we've blinded
        thisTotal.Add(blindedFail)

        Rpf = blindedPass.Clone('Rpf')
        Rpf.Divide(blindedFail)

    # Otherwise just get the Rpf
    else:
        Rpf = thisPass.Clone('Rpf')
        Rpf.Divide(thisFail)


    # Now we need to figure out how to slice up the y-axis based on the total statistics
    # (pass+fail) for each 2D bin

    # Get an average number of events per slice
    total_events = thisTotal.Integral()
    events_per_slice = total_events/float(nslices)

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

    new_y_bins_array = array('d',new_y_bins)

    #################
    # Rebin the Rpf #
    #################

    rebinnedRpf = TH2F('rebinnedRpf','rebinnedRpf',Rpf.GetNbinsX(),xmin,xmax,len(new_y_bins_array)-1,new_y_bins_array)
    rebinnedRpf.Sumw2()


    for xbin in range(1,Rpf.GetNbinsX()+1):
        newBinContent = 0
        newBinErrorSq = 0
        rebinHistYBin = 1
        for ybin in range(1,Rpf.GetNbinsY()+1):
            # If upper edge of old Rpf ybin is < upper edge of rebinHistYBin then add the Rpf bin to the count
            if Rpf.GetYaxis().GetBinUpEdge(ybin) < rebinnedRpf.GetYaxis().GetBinUpEdge(rebinHistYBin):
                newBinContent += Rpf.GetBinContent(xbin,ybin)
                newBinErrorSq += Rpf.GetBinError(xbin,ybin)**2
            # If ==, add to newBinContent, assign newBinContent to current rebinHistYBin, move to the next rebinHistYBin, and restart newBinContent at 0
            elif Rpf.GetYaxis().GetBinUpEdge(ybin) == rebinnedRpf.GetYaxis().GetBinUpEdge(rebinHistYBin):
                newBinContent += Rpf.GetBinContent(xbin,ybin)
                newBinErrorSq += Rpf.GetBinError(xbin,ybin)**2
                rebinnedRpf.SetBinContent(xbin, rebinHistYBin, newBinContent)
                rebinnedRpf.SetBinError(xbin, rebinHistYBin, sqrt(newBinErrorSq))# NEED TO SET BIN ERRORS
                rebinHistYBin += 1
                newBinContent = 0
                newBinErrorSq = 0
            else:
                print 'ERROR when doing psuedo-2D fit approximation. Slices do not line up on y bin edges'
                quit()

    makeCan('rebinCompare',tag,[rebinnedRpf,Rpf])

    # Determine fit function (a polynomial) from the inputConfig
    xFuncString = ''
    nxcoeffs = 0
    for xcoeff in inputConfig['FIT'].keys():
        if xcoeff.find('X') != -1:                                           # For each X*Y0
            if xcoeff.find('Y0') != -1:
                powerIndex = xcoeff[xcoeff.find('X')+1]
                xFuncString += '['+str(powerIndex)+']*x**'+str(powerIndex)+'+'
                nxcoeffs += 1
    # Will have a trailing '+' that needs to be removed
    xFuncString = xFuncString[:-1]

    print 'Will use ' + xFuncString + ' for the x axis fit'

    fitOutput = TFile('fitOutput.root',"RECREATE")
    fitOutput.cd()



    ##################################
    # Now do the fit in the y slices #
    ##################################
    xFunc = TF1('xFunc',xFuncString,xmin,xmax)
    fitResults = {}

    # Book TGraphs to store the fit results as a function of y bins
    for xcoeff in range(nxcoeffs):
        fitResults['xcoeff_'+str(xcoeff)+'_vs_y'] = TH1F('xcoeff_'+str(xcoeff)+'_vs_y','xcoeff_'+str(xcoeff)+'_vs_y',len(new_y_bins_array)-1,new_y_bins_array)

    # Project each y-axis bin and fit along x - save out coefficients to booked tgraph
    for ybin in range(1,rebinnedRpf.GetNbinsY()+1):
        fitResults['fitSlice_'+str(ybin)] = TF1('fitSlice_'+str(ybin),xFuncString,xmin,xmax)
        projX = rebinnedRpf.ProjectionX('rebinnedRpf_sliceX_'+str(ybin),ybin,ybin,'e o')
        projX.Draw('p e')
        binCenter = rebinnedRpf.GetYaxis().GetBinCenter(ybin)
        binLowEdge = rebinnedRpf.GetYaxis().GetBinLowEdge(ybin)
        binHighEdge = rebinnedRpf.GetYaxis().GetBinUpEdge(ybin)
        projX.Fit(xFunc,'E')
        projX.Draw('p e')
        for ix in range(nxcoeffs):
            print 'y' + str(ybin)+ ' xcoeff'+str(ix)+': ' + str(xFunc.GetParameter(ix))
            fitResults['xcoeff_'+str(ix)+'_vs_y'].SetBinContent(ybin,xFunc.GetParameter(ix))
            fitResults['xcoeff_'+str(ix)+'_vs_y'].SetBinError(ybin,xFunc.GetParError(ix))
 
    

            # OLD - USES FitSlicesX
                                            # fitResultsArray = TObjArray()
                                            # rebinnedRpf.FitSlicesX(xFunc,0,-1,0,'',fitResultsArray)

                                            # # Grab fit values from TObjArray
                                            # iterator = fitResultsArray.MakeIterator()
                                            # param = iterator.Next()
                                            # coeffCount = 0
                                            # while param:
                                            #     print param
                                            #     # Just grab the coeff (no ChiSquare, etc) which are the first ncoeff items in the TObjArray
                                            #     if coeffCount <= nxcoeffs-1:
                                            #         param.SetName('slice_X'+str(coeffCount))
                                            #         fitResults['slice_X'+str(coeffCount)] = param
                                            #         coeffCount += 1
                                            #         param.Write()
                                                
                                            #     param = iterator.Next()

    ########################################################
    # And now fit these parameters as a function of y axis #
    ########################################################
    # Build function string
    yFuncString = ''
    nycoeffs = 0
    for ycoeff in inputConfig['FIT'].keys():
        if ycoeff.find('X0Y') != -1:                                           # For each X0Y*
            powerIndex = ycoeff[ycoeff.find('Y')+1]
            yFuncString += '['+str(powerIndex)+']*x**'+str(powerIndex)+'+'
            nycoeffs += 1
    # Will have a trailing '+' that needs to be removed
    yFuncString = yFuncString[:-1]
    print 'Will use ' + yFuncString + ' for the y axis fit'

    # yFunc = TF1('yFunc',yFuncString,ymin,ymax)
    fitResults['xcoeff_0_vs_y'].SetMarkerStyle(8)
    fitResults['xcoeff_0_vs_y'].Draw('p e')

    # Build fit for each parameter distribution along y-axis
    drawList = []
    for xcoeff in range(nxcoeffs):
        fitResults['fit_'+str(xcoeff)] = TF1('yFunc_'+str(xcoeff),yFuncString,ymin,ymax)
        # Do the fit
        fitResults['xcoeff_'+str(xcoeff)+'_vs_y'].Fit(fitResults['fit_'+str(xcoeff)],"E")
        drawList.append(fitResults['xcoeff_'+str(xcoeff)+'_vs_y'])
        # Get and store parameters found
        for iy in range(fitResults['fit_'+str(xcoeff)].GetNpar()):
            fitResults['X'+str(xcoeff)+'Y'+str(iy)] = fitResults['fit_'+str(xcoeff)].GetParameter(iy)
            fitResults['X'+str(xcoeff)+'Y'+str(iy)+'err'] = fitResults['fit_'+str(xcoeff)].GetParError(iy)

    makeCan('xcoeff_v_y',tag,drawList)

    # Store new values in inputConfig
    print 'Resetting fit parameters in input config'
    for key in inputConfig['FIT'].keys():
        if key != 'HELP':
            print 'Changing fit guesses for ' + key
            print str(inputConfig['FIT'][key]['NOMINAL']) + ' -> ' + str(max(0,fitResults[key]))
            print str(inputConfig['FIT'][key]['HIGH']) + ' -> ' + str(fitResults[key]+sigma*fitResults[key+'err'])
            print str(inputConfig['FIT'][key]['LOW']) + ' -> ' + str(max(fitResults[key]-sigma*fitResults[key+'err'],0))
            inputConfig['FIT'][key]['NOMINAL'] = max(0,fitResults[key])
            inputConfig['FIT'][key]['HIGH'] = fitResults[key]+sigma*fitResults[key+'err']
            inputConfig['FIT'][key]['LOW'] = max(fitResults[key]-sigma*fitResults[key+'err'],0)

