import ROOT
from ROOT import *

import sys

import pprint 
pp = pprint.PrettyPrinter(indent=4)

import header
from header import makeCan, FindCommonString

def main(inputConfig, blindData, globalDir, fittype='s', suffix='',procAddString=''):
    allVars = []

    subdir = ''
    batch = False
    if len(globalDir.split('/')) > 1:
        if globalDir.split('/')[1] != '':
            subdir = '_'+globalDir.split('/')[1]
            batch = True


    #####################
    #   Get everything  #
    #####################

    # File with histograms and RooFitResult parameters
    if not batch:
        post_file = TFile.Open(globalDir+'/postfitshapes_'+fittype+'.root')
        fd_file = TFile.Open(globalDir+'/fitDiagnostics.root')


    else:
        post_file = TFile.Open(globalDir.split('/')[0]+'/postfitshapes_'+fittype+'.root')
        fd_file = TFile.Open(globalDir.split('/')[0]+'/fitDiagnostics.root')

    fit_result = fd_file.Get('fit_'+fittype)

    # Binning
    x_low = inputConfig['BINNING']['X']['LOW']
    x_high = inputConfig['BINNING']['X']['HIGH']
    x_nbins = inputConfig['BINNING']['X']['NBINS']
    x_name = inputConfig['BINNING']['X']['NAME']
    x_binWidth = float(x_high-x_low)/float(x_nbins)
    try:
        x_title = inputConfig['BINNING']['X']['TITLE']
    except:
        x_title = ''

    sigstart = inputConfig['BINNING']['X']['SIGSTART']
    sigend = inputConfig['BINNING']['X']['SIGEND']

    # Get the corresponding bin numbers with sigstart and end edges
    for xbin in range(1,x_nbins+1):
        if post_file.Get('pass'+subdir+'_prefit/data_obs').GetXaxis().GetBinLowEdge(xbin) == sigstart:
            sigstart_bin = xbin
        if post_file.Get('pass'+subdir+'_prefit/data_obs').GetXaxis().GetBinUpEdge(xbin) == sigend:
            sigend_bin = xbin

    y_low = inputConfig['BINNING']['Y']['LOW']
    y_high = inputConfig['BINNING']['Y']['HIGH']
    y_nbins = inputConfig['BINNING']['Y']['NBINS']
    y_name = inputConfig['BINNING']['Y']['NAME']
    try:
        y_title = inputConfig['BINNING']['Y']['TITLE']
    except:
        y_title = ''

    # Define low, middle, high projection regions for y (x regions defined already via signal region bounds)
    y_turnon_endBin = post_file.Get('pass'+subdir+'_prefit/data_obs').ProjectionY().GetMaximumBin()
    y_tail_beginningBin = int((y_nbins - y_turnon_endBin)/2.0 + y_turnon_endBin)

    if y_turnon_endBin > y_nbins/2.0:  # in case this isn't a distribution with a turn-on
        y_turnon_endBin = int(round(y_nbins/3.0))
        y_tail_beginningBin = 2*y_turnon_endBin
    y_turnon_endVal = str(int(post_file.Get('pass'+subdir+'_prefit/data_obs').GetYaxis().GetBinUpEdge(y_turnon_endBin)))
    y_tail_beginningVal = str(int(post_file.Get('pass'+subdir+'_prefit/data_obs').GetYaxis().GetBinLowEdge(y_tail_beginningBin)))
 

    # Final fit signal strength
    tree_fit_sb = fd_file.Get('tree_fit_sb')
    tree_fit_sb.GetEntry(0)
    signal_strength = tree_fit_sb.r

    #####################
    #    Data vs Bkg    #
    #####################

    hist_dict = {}

    # Organize and make any projections, etc
    for process in [process for process in inputConfig['PROCESS'] if process != 'HELP']+['qcd']:
        hist_dict[process] = {}
        for cat in ['fail','pass']:
            hist_dict[process][cat] = {}

            hist_dict[process][cat]['prefit_2D'] = post_file.Get(cat +subdir+ '_prefit/'+process).Clone()
            hist_dict[process][cat]['postfit_2D'] = post_file.Get(cat +subdir+ '_postfit/'+process).Clone()

            hist_dict[process][cat]['prefit_2D'].SetMinimum(0)
            hist_dict[process][cat]['postfit_2D'].SetMinimum(0)
            hist_dict[process][cat]['prefit_2D'].SetTitle(process + ', ' + cat +subdir+ ', pre-fit')
            hist_dict[process][cat]['postfit_2D'].SetTitle(process + ', ' + cat +subdir+ ', post-fit')

            base_proj_name_pre = process+'_'+cat+subdir+'_pre_'
            base_proj_name_post = process+'_'+cat+subdir+'_post_'

            hist_dict[process][cat]['prefit_projx1'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+str(y_low)+'-'+y_turnon_endVal,              1,                      y_turnon_endBin, 'e')
            hist_dict[process][cat]['prefit_projx2'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+y_turnon_endVal+'-'+y_tail_beginningVal,     y_turnon_endBin+1,      y_tail_beginningBin,'e')
            hist_dict[process][cat]['prefit_projx3'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+y_tail_beginningVal+'-'+str(y_high),         y_tail_beginningBin+1,  y_nbins,'e')

            hist_dict[process][cat]['prefit_projy1'] = hist_dict[process][cat]['prefit_2D'].ProjectionY(base_proj_name_pre+'projy_'+str(x_low)+'-'+str(sigstart),                1,                      sigstart_bin-1,'e')
            hist_dict[process][cat]['prefit_projy2'] = hist_dict[process][cat]['prefit_2D'].ProjectionY(base_proj_name_pre+'projy_'+str(sigstart)+'-'+str(sigend),               sigstart_bin,           sigend_bin,'e')
            hist_dict[process][cat]['prefit_projy3'] = hist_dict[process][cat]['prefit_2D'].ProjectionY(base_proj_name_pre+'projy_'+str(sigend)+'-'+str(x_high),                 sigend_bin+1,           x_nbins,'e')

            hist_dict[process][cat]['postfit_projx1'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+str(y_low)+'-'+y_turnon_endVal,           1,                      y_turnon_endBin, 'e')
            hist_dict[process][cat]['postfit_projx2'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+y_turnon_endVal+'-'+y_tail_beginningVal,  y_turnon_endBin+1,      y_tail_beginningBin,'e')
            hist_dict[process][cat]['postfit_projx3'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+y_tail_beginningVal+'-'+str(y_high),      y_tail_beginningBin+1,  y_nbins,'e')

            hist_dict[process][cat]['postfit_projy1'] = hist_dict[process][cat]['postfit_2D'].ProjectionY(base_proj_name_post+'projy_'+str(x_low)+'-'+str(sigstart),             1,                      sigstart_bin-1,'e')
            hist_dict[process][cat]['postfit_projy2'] = hist_dict[process][cat]['postfit_2D'].ProjectionY(base_proj_name_post+'projy_'+str(sigstart)+'-'+str(sigend),            sigstart_bin,           sigend_bin,'e')
            hist_dict[process][cat]['postfit_projy3'] = hist_dict[process][cat]['postfit_2D'].ProjectionY(base_proj_name_post+'projy_'+str(sigend)+'-'+str(x_high),              sigend_bin+1,           x_nbins,'e')


            x_edges = [x_low,sigstart,sigend,x_high]
            y_edges = [y_low,y_turnon_endVal,y_tail_beginningVal,y_high]

            for z in ['x','y']:
                for i in range(1,4):
                    hist_dict[process][cat]['postfit_proj'+z+str(i)].SetMinimum(0)
                    if z == 'x':
                        hist_dict[process][cat]['postfit_proj'+z+str(i)].SetTitle(process + ', ' + cat+subdir + ', ' +str(y_edges[i-1]) +'-'+ str(y_edges[i]))
                    elif z == 'y':
                        hist_dict[process][cat]['postfit_proj'+z+str(i)].SetTitle(process + ', ' + cat+subdir + ', ' +str(x_edges[i-1]) +'-'+ str(x_edges[i]))

    post_file.Close()

    # Add together processes that we want to see as one
    process_list = hist_dict.keys() # save this real quick
    colors = []
    if procAddString != '':
        for summation in procAddString.split(','):  # For each set we're add together
            totalProcName = FindCommonString(summation.split('+'))  # Get the common name

            process_list.append(totalProcName)   # add the name to a list so we can keep track
            hist_dict[totalProcName] = {}
            inputConfig["PROCESS"][totalProcName] = {"COLOR":inputConfig["PROCESS"][summation.split('+')[0]]["COLOR"],
                                                     "CODE":inputConfig["PROCESS"][summation.split('+')[0]]["CODE"] }

            for cat in hist_dict[summation.split('+')[0]].keys():   # for each pass/fail
                hist_dict[totalProcName][cat] = {}
                for reg in hist_dict[summation.split('+')[0]][cat].keys():  # and each region
                    firstHist = hist_dict[summation.split('+')[0]][cat][reg].Clone(totalProcName+'_'+cat+'_'+reg)    # Clone the "first" histogram and give it the totalProcName
                    for proc in summation.split('+'):   # for each process in the list of ones we're adding together
                        if proc != summation.split('+')[0]: # not the first one since we've cloned that
                            firstHist.Add(hist_dict[proc][cat][reg])    # add it
                    hist_dict[totalProcName][cat][reg] = firstHist    # Put it in the hist_dict


            for proc in summation.split('+'):
                process_list.remove(proc)

    # Create lists for the 2D projections (ones you want to see together)
    for process in hist_dict.keys():    # Canvas
        isSignal = (process != 'qcd' and inputConfig['PROCESS'][process]['CODE'] == 0)
        twoDList = []         
        for cat in ['fail','pass']:
            for fit in ['prefit', 'postfit']:
                if isSignal and fittype == 's':
                    hist_dict[process][cat][fit+'_2D'].Scale(signal_strength)
                    twoDList.append(hist_dict[process][cat][fit+'_2D'])
                else:
                    twoDList.append(hist_dict[process][cat][fit+'_2D'])

        if isSignal and fittype != 's':
            continue
        else:
            makeCan(process+'_fit'+fittype+'_2D',globalDir+'/',twoDList,xtitle=x_title,ytitle=y_title)

    process_list[-1],process_list[-2] = process_list[-2],process_list[-1]

    # Get the colors
    for process in process_list:
        if process != 'data_obs':
            if process not in inputConfig['PROCESS'].keys():
                continue
            if (process != 'qcd' and inputConfig['PROCESS'][process]['CODE'] != 0):
                if (process in inputConfig['PROCESS'].keys()) and ('COLOR' in inputConfig['PROCESS'][process].keys()):
                    colors.append(inputConfig['PROCESS'][process]["COLOR"])
                else:
                    colors.append(None)

    # Create lists for makeCan of the projections
    for plotType in ['postfit_projx','postfit_projy']:   # Canvases
        bkgList = []
        dataList = []
        signal_list = []
        
        for cat in ['fail','pass']: # Row 
            for regionNum in range(1,4):    # Column
                bkg_process_list = []
                for process in process_list:
                    if process != 'data_obs':
                        if (process != 'qcd' and inputConfig['PROCESS'][process]['CODE'] != 0):
                            bkg_process_list.append(hist_dict[process][cat][plotType+str(regionNum)])
                        elif (process != 'qcd' and inputConfig['PROCESS'][process]['CODE'] == 0):
                            hist_dict[process][cat][plotType+str(regionNum)].Scale(signal_strength)
                            signal_list.append(hist_dict[process][cat][plotType+str(regionNum)])
                        
                    else:
                        dataList.append(hist_dict[process][cat][plotType+str(regionNum)])
                # print colors
                # print bkg_process_list
                # raw_input('waiting')

                bkg_process_list.append(hist_dict['qcd'][cat][plotType+str(regionNum)]) # QCD goes last
                colors.append(kYellow)
                bkgList.append(bkg_process_list)


        if 'x' in plotType:
            makeCan(plotType+'_fit'+fittype,globalDir+'/',dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=x_title)
            makeCan(plotType+'_fit'+fittype+'_log',globalDir+'/',dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=x_title,logy=True)
        elif 'y' in plotType:
            makeCan(plotType+'_fit'+fittype,globalDir+'/',dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=y_title)
            makeCan(plotType+'_fit'+fittype+'_log',globalDir+'/',dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=y_title,logy=True)

    ##############
    #    Rp/f    #
    ##############

    # Need to sample the space to get the Rp/f with proper errors (1000 samples)
    rpf_xnbins = x_nbins
    rpf_ynbins = y_nbins
    rpf_samples = TH3F('rpf_samples','rpf_samples',rpf_xnbins, x_low, x_high, rpf_ynbins, y_low, y_high, 1000,0,1)# TH3 to store samples
    sample_size = 100
    # First figure out the names of parameters and their form so we can evaluate them as we sample
    if 'FORM' in inputConfig['FIT'].keys():
        polXO = 0
        polYO = 0
        for param_name in [key for key in inputConfig['FIT'].keys() if key != 'HELP' and key != 'FORM']:
            # Assuming poly order is a single digit (pretty reasonable I think...)
            tempXorder = int(param_name[param_name.find('X')+1])
            tempYorder = int(param_name[param_name.find('Y')+1])
            if tempXorder > polXO:
                polXO = tempXorder
            if tempYorder > polYO:
                polYO = tempYorder


        # Grab and save final coefficient values
        param_final = fit_result.floatParsFinal()
        coeffs_final = param_final.selectByName('polyCoeff_*'+suffix)        # Another trick here - if suffix='', this will grab everything including those
        if suffix == '':                                                        # polyCoeffs with the suffix. Need to remove those by explicitely grabbing them
            coeffsToRemove_final = param_final.selectByName('polyCoeff_*_*') # and using .remove(collection)
            coeffs_final.remove(coeffsToRemove_final)

        param_out = open(globalDir+'/rpf_params'+suffix+'.txt','w')
        coeffIter_final = coeffs_final.createIterator()
        coeff_final = coeffIter_final.Next()
        while coeff_final:
            param_out.write(coeff_final.GetName()+': ' + str(coeff_final.getValV()) + ' +/- ' + str(coeff_final.getError())+'\n')
            coeff_final = coeffIter_final.Next()
        param_out.close()


        for i in range(sample_size):
            sys.stdout.write('\rSampling '+str(100*float(i)/float(sample_size)) + '%')
            sys.stdout.flush()
            param_sample = fit_result.randomizePars()

            for xbin in range(1,rpf_xnbins+1):
                for ybin in range(1,rpf_ynbins+1):
                    bin_val = 0

                    thisXCenter = rpf_samples.GetXaxis().GetBinCenter(xbin)
                    thisYCenter = rpf_samples.GetYaxis().GetBinCenter(ybin)

                    thisXMapped = (thisXCenter - inputConfig['BINNING']['X']['LOW'])/(inputConfig['BINNING']['X']['HIGH'] - inputConfig['BINNING']['X']['LOW'])
                    thisYMapped = (thisYCenter - inputConfig['BINNING']['Y']['LOW'])/(inputConfig['BINNING']['Y']['HIGH'] - inputConfig['BINNING']['Y']['LOW'])

                    for iy in range(polYO+1):
                        for ix in range(polXO+1):
                            coeff = param_sample.find('polyCoeff_'+'x'+str(ix)+'y'+str(iy)+suffix).getValV()
                            bin_val += coeff*(thisXMapped**ix)*(thisYMapped**iy)

                    # print '['+str(thisXCenter)+','+str(thisYCenter)+'] -> ['+str(thisXMapped)+','+str(thisYMapped)+']: ' +str(bin_val)

                    rpf_samples.Fill(thisXCenter,thisYCenter,bin_val)


    elif 'XFORM' in inputConfig['FIT'].keys() and 'YFORM' in inputConfig['FIT'].keys():
        nxparams = max([int(param[1:]) for param in inputConfig['FIT'].keys() if param.find('X') != -1 and param != 'XFORM'])
        nyparams = max([int(param[1:]) for param in inputConfig['FIT'].keys() if param.find('Y') != -1 and param != 'YFORM'])


        # Grab and save final coefficient values
        param_final = fit_result.floatParsFinal()
        coeffs_final = param_final.selectByName('fitParam*_*'+suffix)        # Another trick here - if suffix='', this will grab everything including those
        if suffix == '':                                                        # polyCoeffs with the suffix. Need to remove those by explicitely grabbing them
            coeffsToRemove_final = param_final.selectByName('fitParam*_*_*') # and using .remove(collection)
            coeffs_final.remove(coeffsToRemove_final)

        param_out = open(globalDir+'/rpf_params'+suffix+'.txt','w')
        coeffIter_final = coeffs_final.createIterator()
        coeff_final = coeffIter_final.Next()
        while coeff_final:
            param_out.write(coeff_final.GetName()+': ' + str(coeff_final.getValV()) + ' +/- ' + str(coeff_final.getError())+'\n')
            coeff_final = coeffIter_final.Next()
        param_out.close()



        # Generate samples
        for i in range(sample_size):
            sys.stdout.write('\rSampling '+str(100*float(i)/float(sample_size)) + '%')
            sys.stdout.flush()
            param_sample = fit_result.randomizePars()

            for xbin in range(1,rpf_xnbins+1):
                for ybin in range(1,rpf_ynbins+1):
                    bin_valx = 0
                    bin_valy = 0

                    thisXCenter = rpf_samples.GetXaxis().GetBinCenter(xbin)
                    thisYCenter = rpf_samples.GetYaxis().GetBinCenter(ybin)

                    thisXMapped = (thisXCenter - inputConfig['BINNING']['X']['LOW'])/(inputConfig['BINNING']['X']['HIGH'] - inputConfig['BINNING']['X']['LOW'])
                    thisYMapped = (thisYCenter - inputConfig['BINNING']['Y']['LOW'])/(inputConfig['BINNING']['Y']['HIGH'] - inputConfig['BINNING']['Y']['LOW'])

                    # Construct the equation from RooRealVar sample
                    for ix in range(0,nxparams):
                        paramName = 'fitParamX_'+str(ix+1)+suffix
                        bin_valx += param_sample.find(paramName).getValV()*thisXMapped**ix

                    for iy in range(0,nyparams):
                        paramName = 'fitParamY_'+str(iy+1)+suffix
                        bin_valy += param_sample.find(paramName).getValV()*thisYMapped**iy

                    bin_val = bin_valx*bin_valy

                    rpf_samples.Fill(thisXCenter,thisYCenter,bin_val)


    elif 'CHEBYSHEV' in inputConfig['FIT'].keys():
        # Import the basis shapes
        cheb_shapes = TFile.Open(globalDir+'basis_plots/basis_shapes.root')
        first_shape_name = cheb_shapes.GetListOfKeys().First().GetName()
        first_shape = cheb_shapes.Get(first_shape_name) # just used to grab binning and such
        cheb_xnbins = first_shape.GetNbinsX()
        cheb_xmin = first_shape.GetXaxis().GetXmin()
        cheb_xmax = first_shape.GetXaxis().GetXmax()
        cheb_ynbins = first_shape.GetNbinsY()
        cheb_ymin = first_shape.GetYaxis().GetXmin()
        cheb_ymax = first_shape.GetYaxis().GetXmax()


        # Grab and save final coefficient values
        param_final = fit_result.floatParsFinal()
        chebCoeffs_final = param_final.selectByName('ChebCoeff_*x*y*'+suffix)        # Another trick here - if suffix='', this will grab everything including those
        if suffix == '':                                                        # polyCoeffs with the suffix. Need to remove those by explicitely grabbing them
            chebCoeffsToRemove_final = param_final.selectByName('ChebCoeff_*x*y*_*') # and using .remove(collection)
            chebCoeffs_final.remove(chebCoeffsToRemove_final)

        param_out = open(globalDir+'/rpf_params'+suffix+'.txt','w')
        chebIter_final = chebCoeffs_final.createIterator()
        chebCoeff_final = chebIter_final.Next()
        while chebCoeff_final:
            param_out.write(chebCoeff_final.GetName()+': ' + str(chebCoeff_final.getValV()) + ' +/- ' + str(chebCoeff_final.getError())+'\n')
            chebCoeff_final = chebIter_final.Next()        
        param_out.close()


        # Make TH3F to store the samples
        rpf_samples = TH3F('rpf_samples','rpf_samples',cheb_xnbins,cheb_xmin,cheb_xmax,cheb_ynbins,cheb_ymin,cheb_ymax, 1000,0,1)# TH3 to store samples

        # Loop over samples
        for i in range(sample_size):
            sys.stdout.write('\rSampling '+str(100*float(i)/float(sample_size)) + '%')
            sys.stdout.flush()

            # Randomize the parameters
            param_sample = fit_result.randomizePars()

            # Make TH2 for this sample
            chebSum = TH2F('chebSum','chebSum',cheb_xnbins,cheb_xmin,cheb_xmax,cheb_ynbins,cheb_ymin,cheb_ymax)

            
            # Grab relevant coefficients and loop over them to sum over the shapes
            chebCoeffs = param_sample.selectByName('ChebCoeff_*x*y*'+suffix)        # Another trick here - if suffix='', this will grab everything including those
            if suffix == '':                                                        # polyCoeffs with the suffix. Need to remove those by explicitely grabbing them
                chebCoeffsToRemove = param_sample.selectByName('ChebCoeff_*x*y*_*') # and using .remove(collection)
                chebCoeffs.remove(chebCoeffsToRemove)

            # Looping...
            chebIter = chebCoeffs.createIterator()
            chebCoeff = chebIter.Next()
            while chebCoeff:
                chebName = chebCoeff.GetName()
                xLabel = chebName[len('chebCoeff_'):len('chebCoeff_')+2] 
                yLabel = chebName[len('chebCoeff_'+xLabel):len('chebCoeff_'+xLabel)+2]

                # Grab and scale the basis shape
                tempScaled = cheb_shapes.Get('cheb_Tx'+xLabel[-1]+'_Ty'+yLabel[-1]).Clone()
                tempScaled.Scale(chebCoeffs.find(chebName).getValV())

                # Add to the sum
                chebSum.Add(tempScaled)
                chebCoeff = chebIter.Next()

            for xbin in range(1,chebSum.GetNbinsX()+1):
                for ybin in range(1,chebSum.GetNbinsY()+1):
                    thisXCenter = rpf_samples.GetXaxis().GetBinCenter(xbin)
                    thisYCenter = rpf_samples.GetYaxis().GetBinCenter(ybin)
                    rpf_samples.Fill(thisXCenter,thisYCenter,chebSum.GetBinContent(xbin,ybin))

            del chebSum

    print '\n'
    rpf_final = TH2F('rpf_final','rpf_final',rpf_samples.GetNbinsX(),rpf_samples.GetXaxis().GetXmin(),rpf_samples.GetXaxis().GetXmax(),rpf_samples.GetNbinsY(),rpf_samples.GetYaxis().GetXmin(),rpf_samples.GetYaxis().GetXmax())
    # Now loop over all x,y bin in rpf_samples, project onto Z axis, 
    # get the mean and RMS and set as the bin content and error in rpf_final
    for xbin in range(1,rpf_final.GetNbinsX()+1):
        for ybin in range(1,rpf_final.GetNbinsY()+1):
            temp_projz = rpf_samples.ProjectionZ('temp_projz',xbin,xbin,ybin,ybin)
            rpf_final.SetBinContent(xbin,ybin,temp_projz.GetMean())
            rpf_final.SetBinError(xbin,ybin,temp_projz.GetRMS())

    rpf_file = TFile.Open(globalDir+'/plots/rpf_comparisons.root','UPDATE')
    

    rpf_c = TCanvas('rpf_c','Post-fit R_{P/F}',800,700)
    rpf_final.Draw('lego')
    rpf_c.Print(globalDir+'/plots/fit_'+fittype+'/postfit_rpf_lego.pdf','pdf')
    rpf_final.Draw('pe')
    rpf_c.Print(globalDir+'/plots/fit_'+fittype+'/postfit_rpf_errs.pdf','pdf')

    # Do a ratio and diff with pre-fit
    prefit_rpf = rpf_file.Get('rebinnedRpf')

    # Ratio
    rpf_ratio = rpf_final.Clone('rpf_ratio_fit'+fittype)
    rpf_ratio.Divide(prefit_rpf)

    # Difference
    rpf_diff = rpf_final.Clone('rpf_diff_fit'+fittype)
    rpf_diff.Add(prefit_rpf,-1)


    # rpf_ratio_c = TCanvas('rpf_ratio_c','Ratio of post-fit to pre-fit R_{P/F}',800,700)
    # rpf_ratio.Draw('surf')
    # rpf_ratio_c.Print(globalDir+'/plots/rpf_post-to-pre_ratio.pdf','pdf')

    rpf_file.cd()
    rpf_final.Write()
    rpf_ratio.Write()
    rpf_diff.Write()

    rpf_file.Close()
