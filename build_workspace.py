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
    x_var,y_var = getRRVs(inputConfig,False)
    var_list = RooArgList(x_var,y_var)

    allVars.extend([x_var,y_var,var_list])


    #########################
    #   Make RooDataHists   #
    #########################

    # It may have seemed crazy to keep this dictionary of TH2s around but it has two great things
    # 1 - structure, 2 - the TH2s we need to make into RDHs
    # However, we will do one thing for convenience - copy it and replace the TH2s in the copy with RDHs
    # if the process has CODE 0,1,2 and a PDF with a normalization if the CODE is 3

    Roo_dict = dictCopy(dictTH2s)
    rateParam_lines = []

    # For procees, cat, dict...
    for process in Roo_dict.keys():
        # Make a normalization for CODE 3 process that floats between 0 and double
        if inputConfig["PROCESS"][process]["CODE"] == 3:
            Roo_dict[process]['NORM'] = RooRealVar(process+'_norm',process+'_norm',0.0,2.0)
        for cat in ['pass','fail']:
            for dist in Roo_dict[process][cat].keys():

                if inputConfig["PROCESS"][process]["CODE"] != 3:
                    Roo_dict[process][cat][dist] = {}
                    Roo_dict[process][cat][dist]['RDH'] = makeRDH(dictTH2s[process][cat][dist],var_list)

                elif inputConfig["PROCESS"][process]["CODE"] == 3:
                    print 'Will have floating process ' + process,         
                    Roo_dict[process][cat][dist] = {}
                    Roo_dict[process][cat][dist]['RDH'] = makeRDH(dictTH2s[process][cat][dist],var_list) 
                    Roo_dict[process][cat][dist]['RDH'].SetName(dictTH2s[process][cat][dist].GetName()+'_RDH') 
                    Roo_dict[process][cat][dist]['RHP'] = makeRHP(Roo_dict[process][cat][dist]['RDH'],var_list)
                    Roo_dict[process][cat][dist]['RHP'].SetName(dictTH2s[process][cat][dist].GetName().replace('_RDH',''))
                    Roo_dict[process][cat][dist]['RHP'].Print()

                    # Make normalization
                    norm_name = Roo_dict[process][cat][dist]['RHP'].GetName() + '_norm'
                    norm_start = RooConstVar(norm_name+'_start',norm_name+'_start',float(dictTH2s[process][cat][dist].Integral()))
                    allVars.append(norm_start)

                    norm = RooProduct(norm_name,norm_name,RooArgList(Roo_dict[process]['NORM'],norm_start))
                    Roo_dict[process][cat][dist]['NORM'] = norm


#############################################################################################
# Everything from here on is only dealing with the QCD estimate - everything else is done   #
#############################################################################################

    ####################################################
    # Get the fit information and store as RooRealVars #
    ####################################################
    xFitForm = inputConfig['FIT']['XFORM']    
    yFitForm = inputConfig['FIT']['YFORM']

    # Do some quick checks to make sure these are formatted correctly
    checkFitForm(xFitForm,yFitForm)


    nxparams = max([int(param[1:]) for param in inputConfig['FIT'].keys() if param.find('X') != -1 ])
    nyparams = max([int(param[1:]) for param in inputConfig['FIT'].keys() if param.find('Y') != -1 ])
    print 'Total number of x fit parameters is ' + str(nxparams)
    print 'Total number of y fit parameters is ' + str(nyparams)

    fitParamRRVs = {}

    # Make and store RRVs for each param
    for nparams in [nxparams, nyparams]:
        if nparams == nxparams:
            thisVar = 'X'
        else:
            thisVar = 'Y'
        for ip in range(nparams+1):
            input_param_vals = inputConfig['FIT'][str(ip)]
            thisNom = input_param_vals['NOMINAL']

            name = 'fitParam'+thisVar+'_'+str(ip)

            if 'LOW' in input_param_vals.keys() and 'HIGH' in input_param_vals.keys():
                thisLow = input_param_vals['LOW']
                thisHigh = input_param_vals['HIGH']
                fitParamRRVs['fitParam'+thisVar+'_'+str(ip)] = RooRealVar(name,name,thisNom,thisLow,thisHigh)
            else:
                input_break = raw_input('WARNING: Upper and lower bounds on global fit parameter ' +str(input_param_vals['GLOBALPARAM']) + ' not defined. Please type "N" if you would like to quit or enter if you would like to treat the parameter as constant')
                if input_break == 'N':
                    quit()
                else:
                    fitParamRRVs['fitParam'+thisVar+'_'+str(ip)] = RooConstVar(name,name,thisNom)

            allVars.append(fitParamRRVs['fitParam'+thisVar+'_'+str(ip)])


    ######################################
    # Build the RooParametricHist2D bins #
    ######################################

    # Start by copying the failing distribution of data
    TH2_data_fail = dictTH2s['data_obs']['fail']['nominal']

    Roo_dict['qcd'] = {}

    binListFail = RooArgList()
    binListPass = RooArgList()

    # Get each bin
    for ybin in range(1,TH2_data_fail.GetNbinsY()+1):
        for xbin in range(1,TH2_data_fail.GetNbinsX()+1):
            # Now that we're in a specific bin, we need to process it
            
            # Now that we know we aren't in the blinded region, make a name for the bin RRV
            name = 'Fail_bin_'+str(xbin)+'-'+str(ybin)

            # First let's check if we are blinded and if so if we're in a blinded region
            upperEdge_in_sig = False
            lowerEdge_in_sig = False
            if blinded:
                upperEdge_in_sig = TH2_data_fail.GetXaxis().GetBinUpEdge(xbin) > inputConfig['BINNING']['X']['SIGSTART'] and TH2_data_fail.GetXaxis().GetBinUpEdge(xbin) <= inputConfig['BINNING']['X']['SIGEND']
                lowerEdge_in_sig = TH2_data_fail.GetXaxis().GetBinLowEdge(xbin) >= inputConfig['BINNING']['X']['SIGSTART'] and TH2_data_fail.GetXaxis().GetBinLowEdge(xbin) < inputConfig['BINNING']['X']['SIGEND']
            

            # If we have zero bins at any point, the fail and pass bins must two properties 
            # 1) Have a negligible number of events but NOT zero (will break Combine)
            # 2) Have no dependence on the Rp/f (otherwise they will influence the Rp/f shape)
            if blinded and (upperEdge_in_sig or lowerEdge_in_sig):
                binRRV = RooRealVar(name, name, 0.001,0.0,0.01)
                binListFail.add(binRRV)
                allVars.append(binRRV)

                thisBinPass = RooRealVar('Pass_bin_'+str(xbin)+'-'+str(ybin), 'Pass_bin_'+str(xbin)+'-'+str(ybin), 0.001,0.0,0.01)
                binListPass.add(thisBinPass)
                allVars.append(thisBinPass)

            # Otherwise
            else:
                # Initialize contents
                binContent = TH2_data_fail.GetBinContent(xbin,ybin)
                binErrUp = binContent + TH2_data_fail.GetBinErrorUp(xbin,ybin)
                binErrDown = binContent - TH2_data_fail.GetBinErrorLow(xbin,ybin)
                
                # Now subtract away the parts that we don't want
                for process in dictTH2s.keys():
                    thisTH2 = dictTH2s[process]['fail']['nominal']

                    # Check the code and change bin content and errors accordingly
                    if inputConfig['PROCESS'][process]['CODE'] == 0: # signal
                        continue
                    elif inputConfig['PROCESS'][process]['CODE'] == 1: # data
                        continue
                    elif inputConfig['PROCESS'][process]['CODE'] == 2: # unchanged MC
                        binContent  = binContent    - thisTH2.GetBinContent(xbin,ybin)
                        binErrUp    = binErrUp      - thisTH2.GetBinContent(xbin,ybin) + thisTH2.GetBinErrorUp(xbin,ybin)              # Just propagator errors normally
                        binErrDown  = binErrDown    - thisTH2.GetBinContent(xbin,ybin) - thisTH2.GetBinErrorLow(xbin,ybin)
                    elif inputConfig['PROCESS'][process]['CODE'] == 3: # renorm MC
                        binErrUp    = binContent                                               # Err up is no MC subtraction
                        binErrDown  = binContent    - 2.0 * thisTH2.GetBinContent(xbin,ybin)   # Err down is double MC subtraction
                        binContent  = binContent    - thisTH2.GetBinContent(xbin,ybin)         # Nominal is nominal MC subtraction

                # If bin content is <= 0, treat this bin as a RooConstVar at 0
                if binContent <= 0:
                    binRRV = RooRealVar(name, name, 0.001,0.0,0.01)
                    binListFail.add(binRRV)
                    allVars.append(binRRV)

                    thisBinPass = RooRealVar('Pass_bin_'+str(xbin)+'-'+str(ybin), 'Pass_bin_'+str(xbin)+'-'+str(ybin), 0.001,0.0,0.01)
                    binListPass.add(thisBinPass)
                    allVars.append(thisBinPass)

                else:
                    binRRV = RooRealVar(name, name, binContent, max(binErrDown,0), max(binErrUp,0))

                    # Store the bin
                    binListFail.add(binRRV)
                    allVars.append(binRRV)

                    # Then get bin center and assign it to a RooConstVar
                    xCenter = TH2_data_fail.GetXaxis().GetBinCenter(xbin)
                    yCenter = TH2_data_fail.GetYaxis().GetBinCenter(ybin)

                    xConst = RooConstVar("ConstVar_x_"+str(xbin)+'_'+str(ybin),"ConstVar_x_"+str(xbin)+'_'+str(ybin),xCenter)
                    yConst = RooConstVar("ConstVar_y_"+str(xbin)+'_'+str(ybin),"ConstVar_y_"+str(xbin)+'_'+str(ybin),yCenter)

                    allVars.append(xConst)
                    allVars.append(yConst)

                    # And now make an Rpf function for this bin 

                    # X                 
                    xParamList = RooArgList(xConst)
                    for xParam in range(nxparams+1):                    
                        xParamList.add(fitParamRRVs['fitParamX_'+str(xParam)])
                    xFormulaName = 'xFunc_Bin_'+str(int(xbin))+"_"+str(int(ybin))
                    xFormulaVar = RooFormulaVar(xFormulaName,xFormulaName,xFitForm.replace('x','@0'),xParamList)
                    allVars.extend([xFormulaVar,xParamList])

                    # Y                 
                    yParamList = RooArgList(yConst)
                    for yParam in range(nyparams+1):                    
                        yParamList.add(fitParamRRVs['fitParamY_'+str(yParam)])
                    yFormulaName = 'yFunc_Bin_'+str(int(xbin))+"_"+str(int(ybin))
                    yFormulaVar = RooFormulaVar(yFormulaName,yFormulaName,yFitForm.replace('y','@0'),yParamList)
                    allVars.extend([yFormulaVar,yParamList])

                    # X*Y
                    fullFormulaName = 'FullFunc_Bin_'+str(int(xbin))+"_"+str(int(ybin))
                    fullFormulaVar = RooProduct(fullFormulaName,fullFormulaName,RooArgList(xFormulaVar,yFormulaVar))
                    allVars.append(fullFormulaVar)

                    # Finally make the pass distribution
                    formulaArgList = RooArgList(binRRV,fullFormulaVar)
                    thisBinPass = RooFormulaVar('Pass_bin_'+str(xbin)+'-'+str(ybin),'Pass_bin_'+str(xbin)+'-'+str(ybin),"@0*@1",formulaArgList)
                    binListPass.add(thisBinPass)
                    allVars.append(thisBinPass)


    print "Making RPH2Ds"
    Roo_dict['qcd']['fail'] = {}
    Roo_dict['qcd']['pass'] = {}

    Roo_dict['qcd']['fail']['RPH2D'] = RooParametricHist2D('qcd_fail','qcd_fail',x_var, y_var, binListFail, TH2_data_fail)
    Roo_dict['qcd']['fail']['norm']  = RooAddition('qcd_fail_norm','qcd_fail_norm',binListFail)
    Roo_dict['qcd']['pass']['RPH2D'] = RooParametricHist2D('qcd_pass','qcd_pass',x_var, y_var, binListPass, TH2_data_fail)
    Roo_dict['qcd']['pass']['norm']  = RooAddition('qcd_pass_norm','qcd_pass_norm',binListPass)


    print "Making workspace..."
    # Make workspace to save in
    myWorkspace = RooWorkspace("w_2D")
    for process in Roo_dict.keys():
        for cat in [k for k in Roo_dict[process].keys() if k.find('file') == -1]:
            if cat == 'NORM':
                # continue
                print "Importing " + Roo_dict[process][cat].GetName() + ' from ' + process + ', ' + cat + ', ' + dist
                getattr(myWorkspace,'import')(Roo_dict[process][cat],RooFit.RecycleConflictNodes(),RooFit.Silence())
            else:
                for dist in Roo_dict[process][cat].keys():
                    rooObj = Roo_dict[process][cat][dist]
                    try: 
                        print "Importing " + rooObj.GetName() + ' from ' + process + ', ' + cat + ', ' + dist
                        getattr(myWorkspace,'import')(rooObj,RooFit.RecycleConflictNodes(),RooFit.Silence())
                    except:
                        for itemkey in rooObj.keys():
                            print "Importing " + rooObj[itemkey].GetName() + ' from ' + process + ', ' + cat + ', ' + dist + ', ' + itemkey
                            getattr(myWorkspace,'import')(rooObj[itemkey],RooFit.RecycleConflictNodes(),RooFit.Silence())
                    


    # Now save out the RooDataHists
    myWorkspace.writeToFile('base_'+tag+'.root',True)  

    return myWorkspace#, rateParam_lines



    