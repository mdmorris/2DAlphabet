import ROOT
from ROOT import *
import header
from header import getRRVs, dictStructureCopy, makeCan, copyHistWithNewXbounds, Make_up_down
import pprint
pp = pprint.PrettyPrinter(indent = 2)

def main(inputConfig, organizedDict, blinded, tag):
    allVars = []

    #####################
    #   Get everything  #
    #####################

    # Vars
    x_var,y_var = getRRVs(inputConfig,False)
    var_list = RooArgList(x_var,y_var)

    allVars.extend([x_var,y_var,var_list])

    # Binning
    x_low = inputConfig['BINNING']['X']['LOW']
    x_high = inputConfig['BINNING']['X']['HIGH']
    x_nbins = inputConfig['BINNING']['X']['NBINS']
    sigstart = inputConfig['BINNING']['X']['SIGSTART']
    sigend = inputConfig['BINNING']['X']['SIGEND']
    y_low = inputConfig['BINNING']['Y']['LOW']
    y_high = inputConfig['BINNING']['Y']['HIGH']
    y_nbins = inputConfig['BINNING']['Y']['NBINS']
    x_binWidth = float(x_high-x_low)/float(x_nbins)


    # Open up our files and workspaces
    new_file = TFile.Open(tag+'/MaxLikelihoodFitResult.root')
    new_w = new_file.Get('MaxLikelihoodFitResult')

    old_file = TFile.Open(tag+'/base_'+tag+'.root')    # Need to get the data_obs pass and fail
    old_w = old_file.Get('w_2D')                # which is not stored in the output of Combine
                                                # (the sum of the two is stored there)

    # Build another dictionary that can sort and save our pdfs and RDHs
    # based upon info from the config
    new_roo_dict = dictStructureCopy(organizedDict)

    # Need to add qcd to this
    new_roo_dict['qcd'] = {'pass':0,'fail':0}

    # Grab the RDH for data and PDFs and normalizations for all else
    # For each process, catagory, and distribution...
    for proc in new_roo_dict.keys():
        for cat in new_roo_dict[proc].keys():
            if cat not in ['pass','fail']:                   # Remove any keys that aren't part of catagories
                del new_roo_dict[proc][cat]             # so we don't accidentally access them later
                continue

            # Empty out the old stuff
            new_roo_dict[proc][cat] = {}

            # Grab the correct Combine output

            # There's an important bug fix here. If not renormalizing or using shape based uncertainties
            # then the output of Combine is not a pdf and norm with these names so we need to test if a
            # shape based uncertainty exists for a given process (non-qcd bkg and signal) to determine 
            # what needs to be grabbed
            has_shape_uncert = False
            if proc != 'qcd':
                for syst in inputConfig['PROCESS'][proc]['SYSTEMATICS']:
                    if inputConfig['SYSTEMATIC'][syst]["CODE"] == 2 or inputConfig['SYSTEMATIC'][syst]["CODE"] == 3:
                        has_shape_uncert = True

            # Check for qcd first
            if proc == 'qcd':
                new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeBkg_'+proc+'_'+cat)
                new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_qcd')

            # Now the rest

            elif inputConfig['PROCESS'][proc]['CODE'] == 0:       # If signal
                if has_shape_uncert:
                    new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeSig_'+cat+'_'+proc+'_morph')
                    new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_'+proc)     # normalization
                else:
                    print "Attempting to grab shapeSig_"+proc+'_'+cat+'Pdf'
                    new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeSig_'+proc+'_'+cat+'Pdf')
                    new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_bin'+cat+'_proc_'+proc)
                    new_roo_dict[proc][cat]['PDF'].Print()

            elif inputConfig['PROCESS'][proc]['CODE'] == 1:     # If data
                new_roo_dict[proc][cat]['RDH'] = old_w.data('data_obs_'+cat)

            elif inputConfig['PROCESS'][proc]['CODE'] == 2:     # If unchanged MC bkg
                if has_shape_uncert:
                    new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeBkg_'+cat+'_'+proc+'_morph')
                    new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_'+proc)
                else:
                    print 'Attempting to grab shapeBkg_'+proc+'_'+cat+'Pdf'
                    new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeBkg_'+proc+'_'+cat+'Pdf')
                    new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_bin'+cat+'_proc_'+proc)
                    new_roo_dict[proc][cat]['PDF'].Print()

            elif inputConfig['PROCESS'][proc]['CODE'] == 3:     # If renormalized MC bkg
                if has_shape_uncert:
                    new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeBkg_'+cat+'_'+proc+'_morph')
                else:
                    new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeBkg_'+proc+'_'+cat)
                new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_'+proc)     # normalization


            else: 
                print 'Process ' + proc + ' has code ' + str(inputConfig['PROCESS'][proc]['CODE']) + ' in the input configuration which is not valid. Quitting...'
                quit()

    # Get fit parameters - there's a trick here: if the params float then they're in new_w 
    #                      but if they're RooConstVars then the values need to be grabbed from inputConfig

    # Build a dictionary to store coefficients based on the inputConfig
    PolyCoeffs = {}
    for coeffName in inputConfig['FIT'].keys():
        if coeffName != 'HELP':
            PolyCoeffs[coeffName.lower()] = 0

    # Floating parameters of the fit (store them in python list immediately)
    RAS_rpfParams = new_w.allVars().selectByName('polyCoeff_x*y*',True)
    iter_params = RAS_rpfParams.createIterator()
    RPV_par = iter_params.Next()
    while RPV_par:
        coeffName = RPV_par.GetName()[RPV_par.GetName().find('x'):] # returns "x#y#"
        PolyCoeffs[coeffName] = RPV_par
        # print coeffName + ': ',
        RPV_par.Print()
        allVars.append(RPV_par)
        RPV_par = iter_params.Next()

    # Get remaining constant parameters from old_w
    for coeffName in PolyCoeffs.keys():
        if PolyCoeffs[coeffName] == 0:
            PolyCoeffs[coeffName] = inputConfig['FIT'][coeffName.upper()]['NOMINAL']

    ####################################################################
    #   Do some rebuilding of the polynomial and make a TF2 and histo  #
    ####################################################################

    # Polynomial Order
    polXO = 0
    polYO = 0
    for param_name in PolyCoeffs.keys():
        # Assuming poly order is a single digit (pretty reasonable I think...)
        tempXorder = int(param_name[param_name.find('x')+1])
        tempYorder = int(param_name[param_name.find('y')+1])
        if tempXorder > polXO:
            polXO = tempXorder
        if tempYorder > polYO:
            polYO = tempYorder


    ############
    # Make TF2 #
    ############
    # Build a formula string for the TF2 - needs to be done in this order so that the params are in the right spots
    param_count = 0
    formula_string = ''
    for iy in range(polYO+1):
        for ix in range(polXO+1):
            formula_string += '['+str(param_count)+']*x**'+str(ix)+'*y**'+str(iy)+'+'
            param_count += 1

    # Remove trailing '+' and construct TF2
    formula_string = formula_string[:-1]
    TF2_rpf = TF2("TF2_rpf",formula_string,x_low,x_high,y_low,y_high)

    # Set params and store in text file for later viewing
    param_out = open(tag+'/rpf_params.txt','w')
    param_count = 0
    for iy in range(polYO+1):
        for ix in range(polXO+1):
            if isinstance(PolyCoeffs['x'+str(ix)+'y'+str(iy)], (int, float)):        # For constant parameter
                TF2_rpf.SetParameter(param_count,PolyCoeffs['x'+str(ix)+'y'+str(iy)])
                TF2_rpf.SetParError(param_count,0)
                param_out.write('x'+str(ix)+'y'+str(iy) + ': ' + str(PolyCoeffs['x'+str(ix)+'y'+str(iy)]) + ' +/- 0.0\n')
            else:     # For floating parameter
                TF2_rpf.SetParameter(param_count,PolyCoeffs['x'+str(ix)+'y'+str(iy)].getValV())
                TF2_rpf.SetParError(param_count,PolyCoeffs['x'+str(ix)+'y'+str(iy)].getError())
                param_out.write('x'+str(ix)+'y'+str(iy) + ': ' + str(PolyCoeffs['x'+str(ix)+'y'+str(iy)].getValV()) + ' +/- ' + str(PolyCoeffs['x'+str(ix)+'y'+str(iy)].getError())+'\n')
            param_count += 1              

    param_out.close()

    # Rp/f
    derived_rpf = TCanvas('derived_rpf','derived_rpf',800,700)
    derived_rpf.cd()
    TF2_rpf.SetTitle('Derived R_{P/F}')
    TF2_rpf.Draw('lego')
    derived_rpf.Print(tag+'/plots/derived_rpf.pdf','pdf')


    #################################################################################
    # Old RooPolyVar method that gave too large of errors or improperly scaled hist #
    #################################################################################

    # # Rebuild the RooPolyVar (can't just grab since we have one in each bin stored! Need something over whole space)
    # xPolyList = RooArgList()
    
    # for yCoeff in range(polYO+1):
    #     xCoeffList = RooArgList()

    #     # Get each x coefficient for this y
    #     for xCoeff in range(polXO+1):                    
    #         xCoeffList.add(PolyCoeffs['x'+str(xCoeff)+'y'+str(yCoeff)])


    #     # Make the RooPolyVar in x and save it to the list of x polynomials
    #     thisXPolyVarLabel = "xPol_y_"+str(yCoeff)
    #     xPolyVar = RooPolyVar(thisXPolyVarLabel,thisXPolyVarLabel,x_var,xCoeffList)
    #     xPolyList.add(xPolyVar)
    #     allVars.append(xPolyVar)

    # # Now make a RooPolyVar out of the x polynomials
    # RPV_rpf_func = RooPolyVar("FullPol","FullPol",y_var,xPolyList)
    # allVars.append(RPV_rpf_func)

    # And make a histogram from that
    # VERY IMPORTANT NOTE: You need to call RooFit.Scaling(False) here otherwise it will scale each bin by the xBinWidth*yBinWidth and you'll get huge values. 
    # TH2_rpf_func = RPV_rpf_func.createHistogram("TH2_rpf",x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)))#,RooFit.Scaling(False))



    ###############################################
    #   Make the blinded (closure) distributions  #
    ###############################################
    # If you've blinded the fit, then you'll get the x-axis sideband distributions
    # If you haven't blinded the fit, then you'll get the full x-axis distributions

    final_hists = {}

    # All of this info is inside the Combine output so just need to grab from dictionary and scale if needed
    for proc in new_roo_dict.keys():
        final_hists[proc] = {}
        # First make the histograms
        for cat in new_roo_dict[proc].keys():

            thisDist = new_roo_dict[proc][cat]

            if 'RDH' in thisDist.keys():
                # thisDist['TH2'] = thisDist['RDH'].createHistogram('data_obs_'+cat,x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)))
                thisDist['TH2'] = thisDist['RDH'].createHistogram(proc+'_'+cat,x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)))
                makeCan(proc+'_'+cat,tag,[thisDist['TH2']])

            # PDFs need to be scaled
            else:
                thisDist['TH2'] = thisDist['PDF'].createHistogram(proc +'_'+cat,x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)))
                makeCan(proc +'_'+cat,tag,[thisDist['TH2']])
                if abs(1.0-thisDist['TH2'].Integral()) > 0.001:
                    print 'ERROR: Double check PDF ' + thisDist['PDF'].GetName() + '. It integrated to ' + str(thisDist['TH2'].Integral()) + ' instead of 1'
                
                try:    # inputConfig['PROCESS']['qcd'] doesn't exist so need to try and except (for qcd)
                    if inputConfig['PROCESS'][proc]['CODE'] == 0:   # If signal, need to scale by negative of norm
                        thisDist['TH2'].Scale(thisDist['NORM'].getValV())
                    else:
                        thisDist['TH2'].Scale(thisDist['NORM'].getValV())
                except:
                    thisDist['TH2'].Scale(thisDist['NORM'].getValV())
                makeCan(proc +'_'+cat+'_scaled',tag,[thisDist['TH2']])

        for reg in ['pass','fail']:
            final_hists[proc][reg] = {}
            final_hists[proc][reg]['closure'] = new_roo_dict[proc][reg]['TH2']


    #############################################
    #   Make the signal region distributions    #
    #############################################

    # Take all of the input processes and derive the signal region only
    for proc in organizedDict.keys():
        for reg in ['pass','fail']:
            if blinded: 
                thisUnblindedHist = organizedDict[proc][reg]['nominal_unblinded']
            else:
                thisUnblindedHist = organizedDict[proc][reg]['nominal']
            thisUnblindedHist.Sumw2()
            thisSignalHist = copyHistWithNewXbounds(thisUnblindedHist,proc+'_'+reg+'_signalRegion',x_binWidth,sigstart,sigend)
            final_hists[proc][reg]['signal'] = thisSignalHist
            

    # Now estimate the QCD - need to rebuild from data
    qcd_estimate_fail_signal = final_hists['data_obs']['fail']['signal'].Clone('qcd_fail_signalRegion')

    for proc in organizedDict.keys():
        this_nonqcd = final_hists[proc]['fail']['signal']
        
        # Check the code and change bin content and errors accordingly
        if inputConfig['PROCESS'][proc]['CODE'] == 0: # signal
            continue
        elif inputConfig['PROCESS'][proc]['CODE'] == 1: # data
            continue
        elif inputConfig['PROCESS'][proc]['CODE'] == 2: # unchanged MC
            qcd_estimate_fail_signal.Add(this_nonqcd,-1)
        elif inputConfig['PROCESS'][proc]['CODE'] == 3: # renorm MC
            this_nonqcd.Scale(new_roo_dict[proc]['fail']['NORM'].getValV())
            qcd_estimate_fail_signal.Add(this_nonqcd,-1)

    # Get the Rp/f for the signal region bins
    # TH2_rpf_func_signal = copyHistWithNewXbounds(TH2_rpf_func,TH2_rpf_func.GetName()+'_signalRegion',x_binWidth,sigstart,sigend)

    # Multiply by the qcd fail bins to make the pass estimate
    qcd_estimate_pass_signal = qcd_estimate_fail_signal.Clone('qcd_pass_signalRegion')
    qcd_estimate_pass_signal.Multiply(TF2_rpf)
    qcd_estimate_pass_signal.Draw('p e')

    final_hists['qcd']['fail']['signal'] = qcd_estimate_fail_signal
    final_hists['qcd']['pass']['signal'] = qcd_estimate_pass_signal


    #########################
    #   Plot on canvases    #
    #########################

    for test in ['closure','signal']:
        # Simultaneously build qcd_estimate + non_qcd and data - non_qcd
        full_bkg_fail = final_hists['qcd']['fail'][test].Clone('full_bkg_fail_'+test)
        full_bkg_pass = final_hists['qcd']['pass'][test].Clone('full_bkg_pass_'+test)
        data_minus_nonqcd_fail = final_hists['data_obs']['fail'][test].Clone('data_minus_nonqcd_fail_'+test)
        data_minus_nonqcd_pass = final_hists['data_obs']['pass'][test].Clone('data_minus_nonqcd_pass_'+test)

        # Subtract away nonqcd from data and add nonqcd to the full bkg estimate
        for nonqcd in inputConfig['PROCESS'].keys():
            if nonqcd == 'HELP':
                continue
            elif inputConfig['PROCESS'][nonqcd]['CODE'] == 2 or inputConfig['PROCESS'][nonqcd]['CODE'] == 3:
                full_bkg_fail.Add(final_hists[nonqcd]['fail'][test])
                full_bkg_pass.Add(final_hists[nonqcd]['pass'][test])
                # print 'Subtracting ' + nonqcd + ' from data to estimate QCD'
                data_minus_nonqcd_fail.Add(final_hists[nonqcd]['fail'][test],-1)
                data_minus_nonqcd_pass.Add(final_hists[nonqcd]['pass'][test],-1)
        
        # It's possible that some of the data_minus_nonqcd bins are < 0 so need to fix that
        for xbin in range(1,data_minus_nonqcd_fail.GetNbinsX()+1):
            for ybin in range(1,data_minus_nonqcd_fail.GetNbinsY()+1):
                if data_minus_nonqcd_fail.GetBinContent(xbin,ybin) < 0:
                    data_minus_nonqcd_fail.SetBinContent(xbin,ybin,0)
                if data_minus_nonqcd_pass.GetBinContent(xbin,ybin) < 0:
                    data_minus_nonqcd_pass.SetBinContent(xbin,ybin,0)

        # Titles
        full_bkg_fail.SetTitle('Full Estimate - '+test+' - Fail')
        full_bkg_pass.SetTitle('Full Estimate - '+test+' - Pass')
        data_minus_nonqcd_fail.SetTitle('True QCD - '+test+' - Fail')
        data_minus_nonqcd_pass.SetTitle('True QCD - '+test+' - Pass')
        final_hists['data_obs']['pass'][test].SetTitle('Data - '+test+' - Pass')
        final_hists['data_obs']['fail'][test].SetTitle('Data - '+test+' - Fail')
        final_hists['qcd']['pass'][test].SetTitle('QCD Estimate - '+test+' - Pass')
        final_hists['qcd']['fail'][test].SetTitle('QCD Estimate - '+test+' - Fail')

        # Plot the full and qcd comparisons in 2D
        makeCan('full_comparison_2D',tag,[final_hists['data_obs']['pass'][test],    final_hists['data_obs']['fail'][test],    full_bkg_pass,              full_bkg_fail])
        makeCan('bkg_comparison_2D',tag, [data_minus_nonqcd_pass,         data_minus_nonqcd_fail,         final_hists['qcd']['pass'][test], final_hists['qcd']['fail'][test]])

        # Now make 1D Projections
        data_minus_nonqcd_fail_1D = data_minus_nonqcd_fail.ProjectionY()
        data_minus_nonqcd_pass_1D = data_minus_nonqcd_pass.ProjectionY()
        data_pass_1D = final_hists['data_obs']['pass'][test].ProjectionY()
        data_fail_1D = final_hists['data_obs']['fail'][test].ProjectionY()
        qcd_pass_1D = final_hists['qcd']['pass'][test].ProjectionY()
        qcd_fail_1D = final_hists['qcd']['fail'][test].ProjectionY()
        

        # For each bkg, create pass and fail y projections
        bkgs_fail_1D = []
        bkgs_pass_1D = []
        for bkg in [bkg for bkg in inputConfig['PROCESS'] if bkg != 'HELP' and (inputConfig['PROCESS'][bkg]['CODE'] == 3 or inputConfig['PROCESS'][bkg]['CODE'] == 2)]+['qcd']:
            bkgs_fail_1D.append(final_hists[bkg]['fail'][test].ProjectionY())
            bkgs_pass_1D.append(final_hists[bkg]['pass'][test].ProjectionY())


        # Plot comparisons in 1D
        makeCan('full_comparison_'+test+'_1D',tag,
                [data_pass_1D, data_fail_1D],
                [bkgs_pass_1D, bkgs_fail_1D])
        makeCan('full_comparison_'+test+'_1D',tag,
                [data_pass_1D, data_fail_1D],
                [bkgs_pass_1D, bkgs_fail_1D],
                True)       # True = semilog y
        makeCan('bkg_comparison_'+test+'_1D',tag, 
                [data_minus_nonqcd_pass_1D, data_minus_nonqcd_fail_1D],
                [[qcd_pass_1D], [qcd_fail_1D]])


        # Plot renormalized backrounds (pass, fail, before, after)
        for bkg in [bkg for bkg in inputConfig['PROCESS'] if bkg != 'HELP' and (inputConfig['PROCESS'][bkg]['CODE'] == 3 or inputConfig['PROCESS'][bkg]['CODE'] == 2)]:
            print 'Doing ' + bkg
            # Get some stuff to make copyHistWithNewXbounds call easier to read
            old_dist_pass = organizedDict[bkg]['pass']['nominal']
            old_dist_pass.SetTitle(bkg + ' - Original - Pass')
            old_dist_fail = organizedDict[bkg]['fail']['nominal']
            old_dist_fail.SetTitle(bkg + ' - Original - Fail')

            final_hists[bkg]['pass'][test].SetTitle(bkg + ' - Renormalized - Pass')
            final_hists[bkg]['fail'][test].SetTitle(bkg + ' - Renormalized - Fail')

            makeCan(bkg+'_'+test+'_distributions',tag,[old_dist_pass,old_dist_fail,final_hists[bkg]['pass'][test],final_hists[bkg]['fail'][test]])

            # 1D Projections
            old_dist_pass_1D = old_dist_pass.ProjectionY()
            old_dist_fail_1D = old_dist_fail.ProjectionY()
            bkg_pass_1D = final_hists[bkg]['pass'][test].ProjectionY()
            bkg_fail_1D = final_hists[bkg]['fail'][test].ProjectionY()

            makeCan(bkg+'_'+test+'_distributions_1D',tag,
                    [old_dist_pass_1D,old_dist_fail_1D],
                    [[bkg_pass_1D],[bkg_fail_1D]])


        # Quick thing to save out qcd estimates in a rootfile
        bkg_out = TFile(tag+'/qcd_estimate_'+test+'_'+tag+'.root',"RECREATE")
        bkg_out.cd()
        qcd_pass_1D.Write()
        qcd_pass_1D_up, qcd_pass_1D_down = Make_up_down(qcd_pass_1D)
        qcd_pass_1D_up.Write()
        qcd_pass_1D_down.Write()
        bkg_out.Close()


