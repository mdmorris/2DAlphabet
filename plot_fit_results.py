import ROOT
from ROOT import *
import header
from header import getRRVs, dictStructureCopy, make4x4Can, copyHistWithNewXbounds, make1x1Can
import pprint
pp = pprint.PrettyPrinter(indent = 2)

gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)

def main(inputConfig, organizedDict, blinded, tag):
    allVars = []

    #####################
    #   Get everything  #
    #####################

    # Vars
    if blinded:
        x_var,x_var_low,x_var_high,y_var = getRRVs(inputConfig,blinded)
        var_list = RooArgList(x_var,y_var)
        low_var_list = RooArgList(x_var_low,y_var)
        high_var_list = RooArgList(x_var_high,y_var)
        allVars.extend([x_var_low,x_var_high,low_var_list,high_var_list])
        
        signal_region_start = inputConfig['BINNING']['X']['SIGSTART']
        signal_region_end = inputConfig['BINNING']['X']['SIGEND']

    else:
        x_var,y_var = getRRVs(inputConfig,blinded)
        var_list = RooArgList(x_var,y_var)

    allVars.extend([x_var,y_var,var_list])

    # Binning
    x_low = inputConfig['BINNING']['X']['LOW']
    x_high = inputConfig['BINNING']['X']['HIGH']
    x_nbins = inputConfig['BINNING']['X']['NBINS']
    y_low = inputConfig['BINNING']['Y']['LOW']
    y_high = inputConfig['BINNING']['Y']['HIGH']
    y_nbins = inputConfig['BINNING']['Y']['NBINS']
    x_binWidth = float(x_high-x_low)/float(x_nbins)

    if blinded:
        x_nbins_low = int((signal_region_start - x_low)/x_binWidth)
        x_nbins_high = int((x_high - signal_region_end)/x_binWidth)

    # Open up our files and workspaces
    new_file = TFile.Open('MaxLikelihoodFitResult.root')
    new_w = new_file.Get('MaxLikelihoodFitResult')

    old_file = TFile.Open(tag+'/base_'+tag+'.root')    # Need to get the data_obs pass and fail
    old_w = old_file.Get('w_2D')                # which is not stored in the output of Combine
                                                # (the sum of the two is stored there)

    # Build another dictionary that can sort and save our pdfs and RDHs
    # based upon info from the config
    new_roo_dict = dictStructureCopy(organizedDict)

    if blinded:
        catagories = ['passLow','passHigh','failLow', 'failHigh']
    else:
        catagories = ['pass','fail']

    # Need to add qcd to this
    new_roo_dict['qcd'] = {}
    for cat in catagories:
        new_roo_dict['qcd'][cat] = 0

    # Grab the RDH for data and PDFs and normalizations for all else
    # For each process, catagory, and distribution...
    for proc in new_roo_dict.keys():
        for cat in new_roo_dict[proc].keys():
            if cat not in catagories:                   # Remove any keys that aren't part of catagories
                del new_roo_dict[proc][cat]             # so we don't accidentally access them later
                continue

            # Empty out the old stuff
            new_roo_dict[proc][cat] = {}

            # Grab the correct Combine output

            # Check for qcd first
            if proc == 'qcd':
                new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeBkg_'+proc+'_'+cat)
                new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_qcd')

            # Now the rest
            elif inputConfig['PROCESS'][proc]['CODE'] == 0:       # If signal
                new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeSig_'+cat+'_'+proc+'_morph')
                new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_'+proc)     # normalization

            elif inputConfig['PROCESS'][proc]['CODE'] == 1:     # If data
                new_roo_dict[proc][cat]['RDH'] = new_w.data('data_obs')

            elif inputConfig['PROCESS'][proc]['CODE'] == 2:     # If unchanged MC bkg
                new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeBkg_'+cat+'_'+proc+'_morph')
                new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_'+proc)

            elif inputConfig['PROCESS'][proc]['CODE'] == 3:     # If renormalized MC bkg
                new_roo_dict[proc][cat]['PDF'] = new_w.pdf('shapeBkg_'+cat+'_'+proc+'_morph')
                new_roo_dict[proc][cat]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_'+proc)     # normalization

            else: 
                print 'Process ' + proc + ' has code ' + str(inputConfig['PROCESS'][proc]['CODE']) + ' in the input configuration which is not valid. Quitting...'
                quit()

    # Get fit parameters
    # Parameters of the fit (store them in python list immediately)
    RAS_rpfParams = new_w.allVars().selectByName('polyCoeff_x*y*',True)
    iter_params = RAS_rpfParams.createIterator()
    PolyCoeffs = {}
    RPV_par = iter_params.Next()
    while RPV_par:
        coeffName = RPV_par.GetName()[RPV_par.GetName().find('x'):] # returns "x#y#"
        PolyCoeffs[coeffName] = RPV_par
        # print coeffName + ': ',
        RPV_par.Print()
        allVars.append(RPV_par)
        RPV_par = iter_params.Next()


    ############################################################
    #   Do some rebuilding of the polynomial and make a histo  #
    ############################################################

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

    # Rebuild the RooPolyVar (can't just grab since we have one in each bin stored! Need something over whole space)
    xPolyList = RooArgList()
    for yCoeff in range(polYO+1):
        xCoeffList = RooArgList()

        # Get each x coefficient for this y
        for xCoeff in range(polXO+1):                    
            xCoeffList.add(PolyCoeffs['x'+str(xCoeff)+'y'+str(yCoeff)])

        # Make the RooPolyVar in x and save it to the list of x polynomials
        thisXPolyVarLabel = "xPol_y_"+str(yCoeff)
        xPolyVar = RooPolyVar(thisXPolyVarLabel,thisXPolyVarLabel,x_var,xCoeffList)
        xPolyList.add(xPolyVar)
        allVars.append(xPolyVar)

    # Now make a RooPolyVar out of the x polynomials
    RPV_rpf_func = RooPolyVar("FullPol","FullPol",y_var,xPolyList)
    allVars.append(RPV_rpf_func)

    # And make a histogram from that
    # VERY IMPORTANT NOTE: You need to call RooFit.Scaling(False) here otherwise it will scale each bin by the xBinWidth*yBinWidth and you'll get huge values
    TH2_rpf_func = RPV_rpf_func.createHistogram("Rpf_func",x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)),RooFit.Scaling(False))


    ############################
    #   Start making the rest  #
    ############################

    final_hists = {}

    for proc in new_roo_dict.keys():
        final_hists[proc] = {}
        # First make the histograms
        for cat in new_roo_dict[proc].keys():
            thisDist = new_roo_dict[proc][cat]
            if cat.find('Low') != -1:
                this_x_var = x_var_low
                this_x_nbins = x_nbins_low
            elif cat.find('High') != -1:
                this_x_var = x_var_high
                this_x_nbins = x_nbins_high
            else:
                this_x_var = x_var
                this_x_nbins = x_nbins

            if proc == 'data_obs':
                # thisDist['TH2'] = thisDist['RDH'].createHistogram('data_obs_'+cat,x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)))
                thisDist['RDH'].Print()
                thisDist['TH2'] = thisDist['RDH'].createHistogram(this_x_var,y_var,this_x_nbins,y_nbins,'','data_obs_'+cat)
                make1x1Can('data_obs_'+cat,thisDist['TH2'])

            # PDFs need to be scaled
            else:
                thisDist['TH2'] = thisDist['PDF'].createHistogram(proc +'_'+cat,this_x_var,RooFit.Binning(x_nbins,x_low,x_high),RooFit.YVar(y_var,RooFit.Binning(y_nbins,y_low,y_high)))
                make1x1Can(proc +'_'+cat,thisDist['TH2'])
                if abs(1.0-thisDist['TH2'].Integral()) > 0.001:
                    print 'ERROR: Double check PDF ' + thisDist['PDF'].GetName() + '. It integrated to ' + str(thisDist['TH2'].Integral()) + ' instead of 1'
                thisDist['TH2'].Scale(thisDist['NORM'].getValV())
                make1x1Can(proc +'_'+cat+'_scaled',thisDist['TH2'])

        # If blinded, add the lows and highs and store in final_hists
        if blinded:
            for reg in ['pass','fail']:
                print proc + ' ' + reg
                # Get some stuff to make copyHistWithNewXbounds call easier to read
                new_xbin_width = organizedDict['data_obs']['passLow']['nominal'].GetXaxis().GetBinWidth(1)      # Doesn't really matter what we grab here. Anything should work
                new_x_low = new_roo_dict[proc][reg + 'Low']['TH2'].GetXaxis().GetXmin()
                new_x_high = new_roo_dict[proc][reg + 'High']['TH2'].GetXaxis().GetXmax()

                # Grab the low and high for the region (pass or fail)
                thisTH2_low = new_roo_dict[proc][reg + 'Low']['TH2']
                thisTH2_high = new_roo_dict[proc][reg + 'High']['TH2']

                # Expand the binning to the full range
                thisTH2_low_rebinned = copyHistWithNewXbounds(thisTH2_low, proc + '_'+reg+'Low_original', new_xbin_width, new_x_low, new_x_high)
                thisTH2_high_rebinned = copyHistWithNewXbounds(thisTH2_high, proc + '_'+reg+'High_original', new_xbin_width, new_x_low, new_x_high)
                make1x1Can(proc + '_'+reg+'Low_original',thisTH2_low_rebinned)
                make1x1Can(proc + '_'+reg+'High_original',thisTH2_high_rebinned)

                # Now add them together
                fullTH2 = thisTH2_low_rebinned.Clone(proc + '_'+reg)
                fullTH2.Add(thisTH2_high_rebinned)
                make1x1Can(proc + '_'+reg+'_full',fullTH2)
                final_hists[proc][reg] = fullTH2

        # If not blinded, just copy this hists over to final_hist
        else:
            for reg in ['pass','fail']:
                final_hists[proc][reg] = new_roo_dict[proc][reg]['TH2']

    #########################
    #   Plot on canvases    #
    #########################

    # Simultaneously build qcd_estimate + non_qcd and data - non_qcd
    full_bkg_fail = final_hists['qcd']['fail'].Clone()
    full_bkg_pass = final_hists['qcd']['pass'].Clone()
    data_minus_nonqcd_fail = final_hists['data_obs']['fail'].Clone()
    data_minus_nonqcd_pass = final_hists['data_obs']['pass'].Clone()
    for nonqcd in [nonqcd for nonqcd in final_hists.keys() if nonqcd != 'qcd' and nonqcd != 'data_obs']:
        full_bkg_fail.Add(final_hists[nonqcd]['fail'])
        full_bkg_pass.Add(final_hists[nonqcd]['pass'])
        data_minus_nonqcd_fail.Add(final_hists[nonqcd]['fail'],-1)
        data_minus_nonqcd_pass.Add(final_hists[nonqcd]['pass'],-1)
    
    # Titles
    full_bkg_fail.SetTitle('Full Estimate - Fail')
    full_bkg_pass.SetTitle('Full Estimate - Pass')
    data_minus_nonqcd_fail.SetTitle('True QCD - Fail')
    data_minus_nonqcd_pass.SetTitle('True QCD - Pass')
    final_hists['data_obs']['pass'].SetTitle('Data - Pass')
    final_hists['data_obs']['fail'].SetTitle('Data - Fail')
    final_hists['qcd']['pass'].SetTitle('QCD Estimate - Pass')
    final_hists['qcd']['fail'].SetTitle('QCD Estimate - Fail')

    make4x4Can('full_comparison_2D',final_hists['data_obs']['pass'],    final_hists['data_obs']['fail'],    full_bkg_pass,              full_bkg_fail)
    make4x4Can('bkg_comparison_2D', data_minus_nonqcd_pass,         data_minus_nonqcd_fail,         final_hists['qcd']['pass'], final_hists['qcd']['fail'])

    # 1D Projections
    full_bkg_fail_1D = full_bkg_fail.ProjectionY()
    full_bkg_pass_1D = full_bkg_pass.ProjectionY()
    data_minus_nonqcd_fail_1D = data_minus_nonqcd_fail.ProjectionY()
    data_minus_nonqcd_pass_1D = data_minus_nonqcd_pass.ProjectionY()
    data_pass_1D = final_hists['data_obs']['pass'].ProjectionY()
    data_fail_1D = final_hists['data_obs']['fail'].ProjectionY()
    qcd_pass_1D = final_hists['qcd']['pass'].ProjectionY()
    qcd_fail_1D = final_hists['qcd']['fail'].ProjectionY()

    make4x4Can('full_comparison_1D',data_pass_1D,               data_fail_1D,               full_bkg_pass_1D,  full_bkg_fail_1D)
    make4x4Can('bkg_comparison_1D', data_minus_nonqcd_pass_1D,  data_minus_nonqcd_fail_1D,  qcd_pass_1D,       qcd_fail_1D)

    # Rp/f
    derived_rpf = TCanvas('derived_rpf','derived_rpf',800,700)
    derived_rpf.cd()
    TH2_rpf_func.SetTitle('Derived R_{P/F}')
    TH2_rpf_func.Draw('surf')
    derived_rpf.Print('derived_rpf.pdf','pdf')

    # Plot renormalized backrounds (pass, fail, before, after)
    for renorm_bkg in [renorm_bkg for renorm_bkg in inputConfig['PROCESS'] if renorm_bkg != 'HELP' and inputConfig['PROCESS'][renorm_bkg]['CODE'] == 3]:
        # Get some stuff to make copyHistWithNewXbounds call easier to read
        new_xbin_width = organizedDict[renorm_bkg]['passLow']['nominal'].GetXaxis().GetBinWidth(1)
        new_x_low = final_hists[renorm_bkg]['pass'].GetXaxis().GetXmin()
        new_x_high = final_hists[renorm_bkg]['pass'].GetXaxis().GetXmax()

        old_dist_passLow = copyHistWithNewXbounds(organizedDict[renorm_bkg]['passLow']['nominal'], renorm_bkg + '_passLos_original', new_xbin_width, new_x_low, new_x_high)
        old_dist_passHigh = copyHistWithNewXbounds(organizedDict[renorm_bkg]['passHigh']['nominal'], renorm_bkg + '_passHigh_original', new_xbin_width, new_x_low, new_x_high)
        old_dist_pass = old_dist_passLow.Clone()
        old_dist_pass.Add(old_dist_passHigh)
        old_dist_pass.SetTitle(renorm_bkg + ' - Original - Pass')

        old_dist_failLow = copyHistWithNewXbounds(organizedDict[renorm_bkg]['failLow']['nominal'], renorm_bkg + '_failLow_original', new_xbin_width, new_x_low, new_x_high)
        old_dist_failHigh = copyHistWithNewXbounds(organizedDict[renorm_bkg]['failHigh']['nominal'], renorm_bkg + '_failHigh_original', new_xbin_width, new_x_low, new_x_high)
        old_dist_fail = old_dist_failLow.Clone()
        old_dist_fail.Add(old_dist_failHigh)
        old_dist_fail.SetTitle(renorm_bkg + ' - Original - Fail')

        final_hists[renorm_bkg]['pass'].SetTitle(renorm_bkg + ' - Renormalized - Pass')
        final_hists[renorm_bkg]['fail'].SetTitle(renorm_bkg + ' - Renormalized - Fail')

        make4x4Can(renorm_bkg+'_distributions',old_dist_pass,old_dist_fail,final_hists[renorm_bkg]['pass'],final_hists[renorm_bkg]['fail'])

        # 1D Projections
        old_dist_pass_1D = old_dist_pass.ProjectionY()
        old_dist_fail_1D = old_dist_fail.ProjectionY()
        renorm_bkg_pass_1D = final_hists[renorm_bkg]['pass'].ProjectionY()
        renorm_bkg_fail_1D = final_hists[renorm_bkg]['fail'].ProjectionY()

        make4x4Can(renorm_bkg+'_distributions_1D',old_dist_pass,old_dist_fail,renorm_bkg_pass_1D,renorm_bkg_fail_1D)

