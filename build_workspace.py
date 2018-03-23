#####################################################################################################
# build_workspace.py - written by Lucas Corcodilos, 3/13/18                                         #
# ---------------------------------------------------------                                         #
# This script builds the RooWorkspace with all of the RooFit objects that Combine needs to run the  #
# fit. Most of the TH2s are just converted to RooDataHists (RDH) but we also build the QCD          #
# background estimate using the custom-made RooParametricHist2D class.                              #
#####################################################################################################

import ROOT
from ROOT import *
import os
import pickle
from optparse import OptionParser
import header
from header import getRRVs,makeRDH,dictCopy,makeRHP
import pprint
pp = pprint.PrettyPrinter(indent = 2)


#########################
#       Start Here      #
#########################
def main(dictTH2s,inputConfig,blinded,tag):

    allVars = []    # This is a list of all RooFit objects made. It never gets used for anything but if the
                    # objects never get saved here, then the python memory management will throw them out
                    # because of conflicts with the RooFit memory management. It's a hack.

    ################################
    # Establish our axis variables #
    ################################
    x_var,x_var_low,x_var_high,y_var = getRRVs(inputConfig) # only have to do this once (make sure this TH2 has correct axis names!)
    var_list = RooArgList(x_var,y_var)
    low_var_list = RooArgList(x_var_low,y_var)
    high_var_list = RooArgList(x_var_high,y_var)

    allVars.extend([x_var,x_var_low,x_var_high,y_var,var_list,low_var_list,high_var_list])

    signal_region_start = inputConfig['BINNING']['X']['SIGSTART']
    signal_region_end   = inputConfig['BINNING']['X']['SIGEND']
    x_nbins = inputConfig['BINNING']['X']['NBINS']

    #########################
    #   Make RooDataHists   #
    #########################

    # It may have seemed crazy to keep this dictionary of TH2s around but it has two great things
    # 1 - structure, 2 - the TH2s we need to make into RDHs
    # However, we will do one thing for convenience - copy it and replace the TH2s in the copy with RDHs
    # if the process has CODE 0,1,2 and a PDF with a normalization if the CODE is 3

    catagories = ['pass','fail']
    if blinded:
        catagories.extend(['passLow','failLow','passHigh','failHigh'])

    Roo_dict = dictCopy(dictTH2s)
    rateParam_lines = []

    # For procees, cat, dict...
    for process in Roo_dict.keys():
        for cat in catagories:
            for dist in Roo_dict[process][cat].keys():

                if cat.find('Low') != -1:
                    this_var_list = low_var_list
                elif cat.find('High') != -1:
                    this_var_list = high_var_list
                else:
                    this_var_list = var_list


                if inputConfig["PROCESS"][process]["CODE"] != 3:
                    Roo_dict[process][cat][dist] = {}
                    Roo_dict[process][cat][dist]['RDH'] = makeRDH(dictTH2s[process][cat][dist],this_var_list)


                elif inputConfig["PROCESS"][process]["CODE"] == 3:         # Allowing normalization of this process to float
                    # Need to renormalize the TH2 to 1 so the norm can dictate the PDF normalization (1*norm=norm)
                    # Normalization will come from rateParam in card
                    dictTH2s[process][cat][dist+'_unitNorm'] =  dictTH2s[process][cat][dist].Clone()
                    dictTH2s[process][cat][dist+'_unitNorm'].Scale(1/dictTH2s[process][cat][dist+'_unitNorm'].Integral())

                    Roo_dict[process][cat][dist] = {}
                    Roo_dict[process][cat][dist]['RDH'] = makeRDH(dictTH2s[process][cat][dist+'_unitNorm'],this_var_list) # Has to be made and saved because the RHP references it
                    Roo_dict[process][cat][dist]['RHP'] = makeRHP(Roo_dict[process][cat][dist]['RDH'],this_var_list)

                    # Make line for the card
                    norm_name = Roo_dict[process][cat][dist]['RHP'].GetName()
                    norm_start = dictTH2s[process][cat][dist].Integral()
                    norm_low = 0.0
                    norm_high = 2.0 * norm_start

                    rateParam_lines.append(norm_name + ' rateParam * ' + process + ' ' + str(norm_start) + ' [0,' + str(norm_high)+ ']')


#############################################################################################
# Everything from here on is only dealing with the QCD estimate - everything else is done   #
#############################################################################################

    ####################################################
    # Get the fit information and store as RooRealVars #
    ####################################################
    # Polynomial Order
    polXO = 0
    polYO = 0
    for param_name in [key for key in inputConfig['FIT'].keys() if key != 'HELP']:
        # Assuming poly order is a single digit (pretty reasonable I think...)
        tempXorder = int(param_name[param_name.find('X')+1])
        tempYorder = int(param_name[param_name.find('Y')+1])
        if tempXorder > polXO:
            polXO = tempXorder
        if tempYorder > polYO:
            polYO = tempYorder

    PolyCoeffs = {}
    for yi in range(polYO+1):
        for xi in range(polXO+1):

            input_param_vals = inputConfig['FIT']['X'+str(xi)+'Y'+str(yi)]
            thisNom = input_param_vals['NOMINAL']
            thisLow = input_param_vals['LOW']
            thisHigh = input_param_vals['HIGH']
            name = 'polyCoeff_'+'x'+str(xi)+'y'+str(yi)

            PolyCoeffs['x'+str(xi)+'y'+str(yi)] = RooRealVar(name,name,thisNom,thisLow,thisHigh)
            allVars.append(PolyCoeffs['x'+str(xi)+'y'+str(yi)])


    ######################################
    # Build the RooParametricHist2D bins #
    ######################################
    
    # If we are doing a blinded search, need to split this up into two parts (low and high x-axis regions) and loop over
    if blinded:
        listToEnumerate = [dictTH2s['data_obs']['failLow']['nominal'],dictTH2s['data_obs']['failHigh']['nominal']]
    else:
        listToEnumerate = dictTH2s['data_obs']['fail']['nominal']

    Roo_dict['qcd'] = {}

    # Loop over Low and High categories if blinded
    for iband, TH2_data_fail in enumerate(listToEnumerate):
        if len(listToEnumerate) == 1:       # Not blinded
            sband = ''
            x_var_band = x_var
        else:                               # Blinded
            if iband == 0:   
                sband = 'Low'
                x_var_band = x_var_low
            elif iband == 1:
                sband = 'High'
                x_var_band = x_var_high

        binListFail = RooArgList()
        binListPass = RooArgList()

        # Get each bin
        for ybin in range(1,TH2_data_fail.GetNbinsY()+1):
            for xbin in range(1,TH2_data_fail.GetNbinsX()+1):
                # Now that we're in a specific bin, we need to process it
                
                # First make a name for the bin RRV
                name = 'Fail'+sband+'_bin_'+str(xbin)+'-'+str(ybin)

                # Initialize contents
                binContent = TH2_data_fail.GetBinContent(xbin,ybin)
                binErrUp = binContent + TH2_data_fail.GetBinErrorUp(xbin,ybin)
                binErrDown = binContent - TH2_data_fail.GetBinErrorLow(xbin,ybin)
                
                # Now subtract away the parts that we don't want - namely backgrounds of CODE 2,3
                for process in dictTH2s.keys():
                    thisTH2 = dictTH2s[process]['fail'+sband]['nominal']

                    # Check the code and change bin content and errors accordingly
                    if inputConfig['PROCESS'][process]['CODE'] == 0: # signal
                        continue
                    elif inputConfig['PROCESS'][process]['CODE'] == 1: # data
                        continue
                    elif inputConfig['PROCESS'][process]['CODE'] == 2: # unchanged MC
                        binContent  = binContent    - thisTH2.GetBinContent(xbin,ybin)
                        binErrUp    = binErrUp      - thisTH2.GetBinContent(xbin,ybin) + thisTH2.GetBinErrorUp(xbin,ybin)              # Just propagator errors normally
                        binErrDown  = binErrDown    - thisTH2.GetBinContent(xbin,ybin) - thisTH2.GetBinErrorLow(xbin,ybin)
                    elif inputConfig['PROCESS'][process]['CODE'] == 3: # MC to be renormalized
                        binContent  = binContent    - thisTH2.GetBinContent(xbin,ybin)
                        binErrUp    = binErrUp                                                  # If CODE 3 MC was normalized to 0              - this will dominated so ignoring
                        binErrDown  = binErrDown    - 2*thisTH2.GetBinContent(xbin,ybin)        # If CODE 3 MC was doubled                        bin errors from these CODE 3 MCs
                

                binRRV = RooRealVar(name, name, binContent, max(binErrDown,0), max(binErrUp,0))
                
                # Store the bin
                binListFail.add(binRRV)
                allVars.append(binRRV)

                # Then get bin center and assign it to a RooConstVar
                xCenter = TH2_data_fail.GetXaxis().GetBinCenter(xbin)
                yCenter = TH2_data_fail.GetYaxis().GetBinCenter(ybin)

                xConst = RooConstVar("ConstVar"+sband+"_x_"+str(xbin)+'_'+str(ybin),"ConstVar"+sband+"_x_"+str(xbin)+'_'+str(ybin),xCenter)
                yConst = RooConstVar("ConstVar"+sband+"_y_"+str(xbin)+'_'+str(ybin),"ConstVar"+sband+"_y_"+str(xbin)+'_'+str(ybin),yCenter)

                allVars.append(xConst)
                allVars.append(yConst)

                # And now make a polynomial for this bin
                xPolyList = RooArgList()
                for yCoeff in range(polYO+1):
                    xCoeffList = RooArgList()

                    # Get each x coefficient for this y
                    for xCoeff in range(polXO+1):                    
                        xCoeffList.add(PolyCoeffs['x'+str(xCoeff)+'y'+str(yCoeff)])

                    # Make the polynomial and save it to the list of x polynomials
                    thisXPolyVarLabel = "xPol"+sband+"_y_"+str(yCoeff)+"_Bin_"+str(int(xbin))+"_"+str(int(ybin))
                    xPolyVar = RooPolyVar(thisXPolyVarLabel,thisXPolyVarLabel,xConst,xCoeffList)
                    xPolyList.add(xPolyVar)
                    allVars.append(xPolyVar)

                # Now make a polynomial out of the x polynomials
                thisYPolyVarLabel = "FullPol"+sband+"_Bin_"+str(int(xbin))+"_"+str(int(ybin))
                thisFullPolyVar = RooPolyVar(thisYPolyVarLabel,thisYPolyVarLabel,yConst,xPolyList)

                allVars.append(thisFullPolyVar)


                # Finally make the pass distribution
                formulaArgList = RooArgList(binRRV,thisFullPolyVar)
                thisBinPass = RooFormulaVar('Pass'+sband+'_bin_'+str(xbin)+'-'+str(ybin),'Pass'+sband+'_bin_'+str(xbin)+'-'+str(ybin),"@0*@1",formulaArgList)
                binListPass.add(thisBinPass)
                allVars.append(thisBinPass)


        print "Making RPH2Ds"
        Roo_dict['qcd']['fail'+sband] = {}
        Roo_dict['qcd']['pass'+sband] = {}

        Roo_dict['qcd']['fail'+sband]['RPH2D'] = RooParametricHist2D('qcd_fail'+sband,'qcd_fail'+sband,x_var_band, y_var, binListFail, TH2_data_fail)
        Roo_dict['qcd']['fail'+sband]['norm']  = RooAddition('qcd_fail'+sband+'_norm','qcd_fail'+sband+'_norm',binListFail)
        Roo_dict['qcd']['pass'+sband]['RPH2D'] = RooParametricHist2D('qcd_pass'+sband,'qcd_pass'+sband,x_var_band, y_var, binListPass, TH2_data_fail)
        Roo_dict['qcd']['pass'+sband]['norm']  = RooAddition('qcd_pass'+sband+'_norm','qcd_pass'+sband+'_norm',binListPass)


    print "Making workspace..."
    # Make workspace to save in
    myWorkspace = RooWorkspace("w_2D")
    for process in Roo_dict.keys():
        for cat in [k for k in Roo_dict[process].keys() if k.find('file') == -1]:
            for dist in Roo_dict[process][cat].keys():
                rooObj = Roo_dict[process][cat][dist]
                try: 
                    print "Importing " + rooObj.GetName() + ' from ' + process + ', ' + cat + ', ' + dist
                    getattr(myWorkspace,'import')(rooObj,RooFit.RecycleConflictNodes())
                except:
                    for itemkey in rooObj.keys():
                        print "Importing " + rooObj[itemkey].GetName() + ' from ' + process + ', ' + cat + ', ' + dist + ', ' + itemkey
                        getattr(myWorkspace,'import')(rooObj[itemkey],RooFit.RecycleConflictNodes())
                


    # Now save out the RooDataHists
    myWorkspace.writeToFile('base_'+tag+'.root',True)  

    return myWorkspace, rateParam_lines



    