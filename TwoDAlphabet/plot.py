import ROOT, os, warnings, pandas, collections, math
from TwoDAlphabet.helpers import get_hist_maximum, set_hist_maximums, executeCmd
from TwoDAlphabet.binning import stitch_hists_in_x, convert_to_events_per_unit, get_min_bin_width
from TwoDAlphabet.ext import tdrstyle, CMS_lumi

class Plotter(object):
    def __init__(self,twoD,fittag,loadExisting=False):
        columns = [
            'process','region','type','title',
            'prefit_2D','postfit_2D',
            'prefit_projx1', 'prefit_projx2', 'prefit_projx3',
            'prefit_projy1', 'prefit_projy2', 'prefit_projy3',
            'postfit_projx1', 'postfit_projx2', 'postfit_projx3',
            'postfit_projy1', 'postfit_projy2', 'postfit_projy3'
        ]

        fd_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
        self.fittag = fittag
        self.fit_results = fd_file.Get('fit_'+fittag)
        self.signal_strength = _get_signal_strength(self.fd_file)
        self.twoD = twoD
        self.yaxis1D_title = 'Events / bin'
        self.df = pandas.DataFrame(columns=columns)
        self.dir = '{t}/plots_fit_{f}'.format(t=self.twoD.tag,f=self.fittag)
        self.slices = {'x': {}, 'y': {}}

        if not loadExisting:
            self._make(fittag)
        else:
            self._load()

    def __del__(self):
        self.df.to_csv('{t}/plots_fit_{f}/df.csv'.format(t=self.twoD.tag,f=self.fittag))
        self.root_out.Close()

    class _entryTracker(object):
        def __init__(self,plotter,p,r,t):
            self.entry = {'process': p, 'region': r, 'type': t, 'title': self.plotter.twoD.GetProcessTitle(p)}
            self.plotter = plotter

        def track(self,k,h):
            self.entry[k] = h.GetName()
            self.plotter.root_out.WriteTObject(h.GetName(),h)

    def _load(self):
        # open pickled dataframe and root_out
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name)
        self.df = pandas.read_csv('%s/df.csv'%self.dir)

    def _make(self):
        # make root_out and fill dataframe
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name,'RECREATE')

        shapes_file = ROOT.TFile.Open('postfitshapes_%s.root'%self.fittag)
        loc_base = '{r}_{c}_{t}/{p}'

        for region in self.twoD.GetRegions():
            binning,_ = self.twoD.GetBinningFor(region)
            blinding = [1] if region in self.twoD.options.blindedPlots else []
            self.slices['x'][region] = {'vals': binning.xSlices,'idxs':binning.xSlices_idxs}
            self.slices['y'][region] = {'vals': binning.ySlices,'idxs':binning.ySlices_idxs}
            
            for process in self.twoD.GetProcesses()+['TotalBkg']:
                if process != 'TotalBkg':
                    color = self.twoD.GetProcessColor(process)
                    proc_type = self.twoD.GetProcessType(process)
                    proc_title = self.twoD.GetProcessTitle(process)
                else:
                    color = ROOT.kBlack
                    proc_type = 'TOTAL'
                    proc_title = 'TotalBkg'

                entry = self._entryTracker(self, process, region, proc_type)
                for time in ['prefit','postfit']:
                    # 2D distributions first
                    out2d_name = '%s_%s_%s_2D'%(process,region,time)
                    low  = shapes_file.Get(loc_base.format(r=region, c='LOW', t=time, p=process))
                    sig  = shapes_file.Get(loc_base.format(r=region, c='SIG', t=time, p=process))
                    high = shapes_file.Get(loc_base.format(r=region, c='HIGH',t=time, p=process))

                    full = stitch_hists_in_x(out2d_name, binning, [low,sig,high], blinded=blinding)
                    full.SetMinimum(0)
                    full.SetTitle('%s, %s, %s'%(process,region,time))
                    # if self.twoD.IsSignal(process) and fittag == 's' and time == 'postfit':
                    #     full.Scale(self.signal_strength)
                    entry.track('%s_2D'%time, full)

                    # Now do projections using the 2D
                    out_proj_name = '{p}_{r}_{t}_proj{x}{i}'
                    for proj in ['X','Y']:
                        slices = self.slices[proj.lower()][region]

                        for islice in range(3):
                            hname = out_proj_name.format(p=process,r=region,t=time,x=proj,i=islice)
                            start,stop = self._get_start_stop(islice,slices['idx'])
                            
                            hslice = getattr(full,'Projection'+proj)(hname,start,stop,'e')
                            #hslice_title = '%s, %s, %s, %s-%s'%(process,region,time,slices['val'][islice],slices['val'][islice+1])
                            self._format_hist(
                                hslice, proc_title,
                                binning.xtitle, binning.ytitle,
                                color, proc_type)

                            # TRACK
                            col_name = '{t}_proj{x}{i}'.format(t=time,x=proj.lower(),i=islice+1)
                            entry.track(col_name, hslice)

                    self.df = self.df.append(entry.entry,ignore_index=True)

        shapes_file.Close()

    def _get_start_stop(i,slice_idxs):
        start = slice_idxs[i]+1
        stop  = slice_idxs[i+1]
        return start, stop

    def _format_hist(self, hslice, title, xtitle, ytitle, color, proc_type):
        if not isinstance(hslice,ROOT.TH2):
            ytitle = self.yaxis1D_title
            if self.twoD.config.options.plotEvtsPerUnit:
                hslice = convert_to_events_per_unit(hslice)
                ytitle = 'Events / %s GeV' % get_min_bin_width(hslice)
        hslice.SetMinimum(0)
        hslice.SetTitle(title)
        hslice.GetXaxis().SetTitle(xtitle)
        hslice.GetYaxis().SetTitle(ytitle)

        if proc_type == 'BKG':
            hslice.SetFillColor(color)
        elif proc_type == 'SIGNAL' or proc_type == 'TOTAL':
            hslice.SetLineColor(color)
        elif proc_type == 'DATA':
            hslice.SetLineColor(color)
            hslice.SetMarkerColor(color)

    def _order_df_on_proc_list(self,df,proclist=[],alphaBottom=True):
        if self.proclist == []:
            if alphaBottom:
                process_order = self.twoD.GetProcesses(onlyNonConfig=True) + self.twoD.GetProcesses(includeNonConfig=False)
            else:
                process_order = self.twoD.GetProcesses(includeNonConfig=False) + self.twoD.GetProcesses(onlyNonConfig=True)
        else:
            process_order = proclist

        process_order_df = pandas.DataFrame({'process':process_order})
        return process_order_df.merge(df,on='process',how='outer')

    def plot_2D_distributions(self):
        for process, group in self.df.groupby('process'):
            out_file_name = 'plots_fit_{f}/{p}_2D'.format(f=self.fittag,p=process)
            hist_list = []
            for _,subgroup in group.groupby('region'):
                hist_list.append(subgroup.iloc[0,'prefit_2D'])
                hist_list.append(subgroup.iloc[0,'postfit_2D'])
            
            out_file_name = 'plots_fit_{f}/{p}_2D'.format(f=self.fittag,p=process)
            MakeCan(out_file_name,hist_list,year=self.twoD.options.year)

    def plot_projections(self, logyFlag):
        for proj in ['postfit_projx','postfit_projy']:
            for _, group in self.df.groupby('region'):
                data, bkgs, totalbkgs, signals = [], [], [], []
                ordered_bkgs = self._order_df_on_proc_list(
                                    group.loc[group.type.eq('BKG')],
                                    alphaBottom=(not logyFlag))
                ordered_signals = self._order_df_on_proc_list(
                                    group.loc[group.type.eq('SIGNAL')],
                                    alphaBottom=(not logyFlag))
                bkgNames =  ordered_bkgs.title.to_list()
                signalNames = ordered_signals.title.to_list()
                for islice in range(3):
                    projn = proj+str(islice)
                    this_data = group.loc[group.type.eq('DATA')][projn].iloc[0]
                    this_totalbkg = group.loc[group.type.eq('TOTAL')][projn].iloc[0]
                    these_bkgs = ordered_bkgs[projn].to_list()
                    if self.twoD.config.options.plotPrefitSigInFitB and self.fittag == 'b':
                        projn = projn.replace('postfit','prefit')
                    these_signals = ordered_signals[projn].to_list()

                    # Add everything to track
                    data.append(this_data)
                    totalbkgs.append(this_totalbkg)
                    bkgs.append(these_bkgs)
                    signals.append(these_signals)
                
                slices = [h.GetName().split(',')[-1].lstrip() for h in data]
                titles = [', '.join(h.GetName().split(', ')[1:]) for h in data]

                out_file_name = 'plots_fit_{f}/{proj}%s'.format(f=self.fittag,proj=proj)
                log_str = '' if logyFlag == False else '_logy'
                MakeCan(
                    out_file_name%log_str,
                    histlist=data,      bkglist=bkgs,
                    totalBkg=totalbkgs, signals=signals,
                    titles=titles,      subtitles=slices,
                    bkgNames=bkgNames,  signalNames=signalNames,
                    addSignals=self.twoD.config.options.haddSignals
                )

    def plot_pre_vs_post(self):
        # Make comparisons for each background process of pre and post fit projections
        for proj in ['projx','projy']:
            for process, group in self.df[~self.df.process.isin(['data_obs','TotalBkg'])].groupby('process'):
                pre_list, post_list, titleList = [], [], []
                for region, subgroup in group.groupby('region'):
                    for islice in range(3):
                        projn = proj+str(islice)
                        slices = self.slices[proj.replace('proj','')]
                        post = subgroup['postfit_'+projn].iloc[0]
                        post.SetLineColor(ROOT.kBlack)
                        post.SetLineStyle(9)
                        pre_list.append([subgroup['prefit_'+projn].iloc[0]])  # in terms of makeCan these are "bkg hists"
                        post_list.append(post)   # and these are "data hists"
                        titleList.append('Pre vs Postfit - %s - %s - [%s,%s]'%(process,region,slices['val'][islice],slices['val'][islice+1]))

                MakeCan(
                    'plots_fit_{f}/{p}_{proj}%s'.format(f=self.fittag,p=process,proj=proj),
                    histlist=post_list, bkglist=pre_list,
                    totalBkg=[b[0] for b in pre_list],
                    titles=titleList,   bkgNames=['Prefit, '+process],
                    dataName='Postfit, '+process,
                    datastyle='histe', year=self.twoD.options.year
                )

    def plot_transfer_funcs(self):
        raise NotImplementedError()
        # # Need to sample the space to get the Rp/f with proper errors (1000 samples)
        # rpf_xnbins = len(self.fullXbins)-1
        # rpf_ynbins = len(self.newYbins)-1
        # if self.rpfRatio == False: rpf_zbins = [i/1000000. for i in range(0,1000001)]
        # else: rpf_zbins = [i/1000. for i in range(0,5001)]
        # rpf_samples = TH3F('rpf_samples','rpf_samples',rpf_xnbins, array.array('d',self.fullXbins), rpf_ynbins, array.array('d',self.newYbins), len(rpf_zbins)-1, array.array('d',rpf_zbins))# TH3 to store samples
        # sample_size = 500

        # # Collect all final parameter values
        # param_final = fit_result.floatParsFinal()
        # coeffs_final = RooArgSet()
        # for v in self.rpf.funcVars.keys():
        #     coeffs_final.add(param_final.find(v))

        # # Now sample to generate the Rpf distribution
        # for i in range(sample_size):
        #     sys.stdout.write('\rSampling '+str(100*float(i)/float(sample_size)) + '%')
        #     sys.stdout.flush()
        #     param_sample = fit_result.randomizePars()

        #     # Set params of the Rpf object
        #     coeffIter_sample = param_sample.createIterator()
        #     coeff_sample = coeffIter_sample.Next()
        #     while coeff_sample:
        #         # Set the rpf parameter to the sample value
        #         if coeff_sample.GetName() in self.rpf.funcVars.keys():
        #             self.rpf.setFuncParam(coeff_sample.GetName(), coeff_sample.getValV())
        #         coeff_sample = coeffIter_sample.Next()

        #     # Loop over bins and fill
        #     for xbin in range(1,rpf_xnbins+1):
        #         for ybin in range(1,rpf_ynbins+1):
        #             bin_val = 0

        #             thisXCenter = rpf_samples.GetXaxis().GetBinCenter(xbin)
        #             thisYCenter = rpf_samples.GetYaxis().GetBinCenter(ybin)

        #             if self.recycleAll:
        #                 # Remap to [-1,1]
        #                 x_center_mapped = (thisXCenter - self.newXbins['LOW'][0])/(self.newXbins['HIGH'][-1] - self.newXbins['LOW'][0])
        #                 y_center_mapped = (thisYCenter - self.newYbins[0])/(self.newYbins[-1] - self.newYbins[0])

        #                 # And assign it to a RooConstVar 
        #                 x_const = RooConstVar("ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,x_center_mapped)
        #                 y_const = RooConstVar("ConstVar_y_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,y_center_mapped)
                        
        #                 # Now get the Rpf function value for this bin 
        #                 self.allVars.append(x_const)
        #                 self.allVars.append(y_const)
        #                 self.rpf.evalRpf(x_const, y_const,xbin,ybin)

        #             # Determine the category
        #             if thisXCenter > self.newXbins['LOW'][0] and thisXCenter < self.newXbins['LOW'][-1]: # in the LOW category
        #                 thisxcat = 'LOW'
        #             elif thisXCenter > self.newXbins['SIG'][0] and thisXCenter < self.newXbins['SIG'][-1]: # in the SIG category
        #                 thisxcat = 'SIG'
        #             elif thisXCenter > self.newXbins['HIGH'][0] and thisXCenter < self.newXbins['HIGH'][-1]: # in the HIGH category
        #                 thisxcat = 'HIGH'

        #             bin_val = self.rpf.getFuncBinVal(thisxcat,xbin,ybin)

        #             rpf_samples.Fill(thisXCenter,thisYCenter,bin_val)

        # print ('\n')
        # rpf_final = TH2F('rpf_final','rpf_final',rpf_xnbins, array.array('d',self.fullXbins), rpf_ynbins, array.array('d',self.newYbins))
        # # Now loop over all x,y bin in rpf_samples, project onto Z axis, 
        # # get the mean and RMS and set as the bin content and error in rpf_final
        # for xbin in range(1,rpf_final.GetNbinsX()+1):
        #     for ybin in range(1,rpf_final.GetNbinsY()+1):
        #         temp_projz = rpf_samples.ProjectionZ('temp_projz',xbin,xbin,ybin,ybin)
        #         rpf_final.SetBinContent(xbin,ybin,temp_projz.GetMean())
        #         rpf_final.SetBinError(xbin,ybin,temp_projz.GetRMS())

        # rpf_final.SetTitle('')
        # rpf_final.GetXaxis().SetTitle(self.xVarTitle)
        # rpf_final.GetYaxis().SetTitle(self.yVarTitle)
        # rpf_final.GetZaxis().SetTitle('R_{P/F}' if self.rpfRatio == False else 'R_{Ratio}')
        # rpf_final.GetXaxis().SetTitleSize(0.045)
        # rpf_final.GetYaxis().SetTitleSize(0.045)
        # rpf_final.GetZaxis().SetTitleSize(0.045)
        # rpf_final.GetXaxis().SetTitleOffset(1.2)
        # rpf_final.GetYaxis().SetTitleOffset(1.5)
        # rpf_final.GetZaxis().SetTitleOffset(1.3)

        # rpf_c = TCanvas('rpf_c','Post-fit R_{P/F}',1000,700)
        # CMS_lumi.lumiTextSize = 0.75
        # CMS_lumi.cmsTextSize = 0.85
        # CMS_lumi.extraText = 'Preliminary'
        # CMS_lumi.CMS_lumi(rpf_c, self.year, 11)
        # rpf_c.SetRightMargin(0.2)
        # rpf_final.Draw('colz')
        # rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_colz.pdf','pdf')
        # rpf_final.Draw('surf')
        # rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_surf.pdf','pdf')
        # rpf_final.Draw('pe')
        # rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_errs.pdf','pdf')

        # rpf_file = TFile.Open(self.projPath+'/plots/postfit_rpf_fit'+fittag+'.root','RECREATE')
        # rpf_file.cd()
        # rpf_final.Write()
        # rpf_file.Close()

def plotAllFitResults(twoD,fittag):
    plotter = Plotter(twoD,fittag)
    plotter.plot_2D_distributions()
    plotter.plot_projections(logyFlag=False)
    plotter.plot_projections(logyFlag=True)
    plotter.plot_pre_vs_post()
    # plotter.plot_transfer_funcs()


def MakeCan(name, tag, histlist, bkglist=[],totalBkg=None,signals=[],
            titles=[],subtitles=[],sliceVar='X',dataName='Data',bkgNames=[],
            signalNames=[],logy=False,rootfile=False,dataOff=False,
            datastyle='pe',year=1, addSignals=True, extraText=''):
    '''histlist is just the generic list but if bkglist is specified (non-empty)
    then this function will stack the backgrounds and compare against histlist as if 
    it is data. The imporant bit is that bkglist is a list of lists. The first index
    of bkglist corresponds to the index in histlist (the corresponding data). 
    For example you could have:
      histlist = [data1, data2]
      bkglist = [[bkg1_1,bkg2_1],[bkg1_2,bkg2_2]]

    Args:
        name ([type]): [description]
        tag ([type]): [description]
        histlist ([type]): [description]
        bkglist (list, optional): [description]. Defaults to [].
        totalBkg ([type], optional): [description]. Defaults to None.
        signals (list, optional): [description]. Defaults to [].
        titles (list, optional): [description]. Defaults to [].
        subtitles (list, optional): [description]. Defaults to [].
        sliceVar (str, optional): [description]. Defaults to 'X'.
        dataName (str, optional): [description]. Defaults to 'Data'.
        bkgNames (list, optional): [description]. Defaults to [].
        signalNames (list, optional): [description]. Defaults to [].
        logy (bool, optional): [description]. Defaults to False.
        rootfile (bool, optional): [description]. Defaults to False.
        dataOff (bool, optional): [description]. Defaults to False.
        datastyle (str, optional): [description]. Defaults to 'pe'.
        year (int, optional): [description]. Defaults to 1.
        addSignals (bool, optional): [description]. Defaults to True.
        extraText (str, optional): [description]. Defaults to ''.

    Raises:
        RuntimeError: [description]
        TypeError: [description]
    '''
    if len(histlist) == 1:
        width = 800; height = 700; padx = 1; pady = 1
    elif len(histlist) == 2:
        width = 1200; height = 700; padx = 2; pady = 1
    elif len(histlist) == 3:
        width = 1800; height = 600; padx = 3; pady = 1
    elif len(histlist) == 4:
        width = 1200; height = 1000; padx = 2; pady = 2
    elif len(histlist) <= 6:
        height = 1600; width = 1200; padx = 2; pady = 3
    elif len(histlist) <= 9:
        height = 1600; width = 1600; padx = 3; pady = 3
    else:
        raise RuntimeError('histlist of size %s not currently supported'%(len(histlist),histlist))

    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetTitleAlign(33)
    ROOT.gStyle.SetTitleX(.77)
    myCan = ROOT.TCanvas(name,name,width,height)
    myCan.Divide(padx,pady)

    # A bunch of empty lists for storage if needed
    stacks, tot_hists_err, tot_hists, legends, mains, subs, pulls, tot_sigs = [], [], [], [], [], [], [], []

    # For each hist/data distribution
    for hist_index, hist in enumerate(histlist):
        # Grab the pad we want to draw in
        myCan.cd(hist_index+1)
        thisPad = myCan.GetPrimitive(name+'_'+str(hist_index+1))
        thisPad.cd(); thisPad.SetRightMargin(0.0); thisPad.SetTopMargin(0.0); thisPad.SetBottomMargin(0.0)

        # If this is a TH2, just draw the lego
        if hist.ClassName().find('TH2') != -1:
            if len(bkglist) > 0:
                raise TypeError('It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.')

            ROOT.gPad.SetLeftMargin(0.15)
            ROOT.gPad.SetRightMargin(0.2)
            ROOT.gPad.SetBottomMargin(0.12)
            ROOT.gPad.SetTopMargin(0.1)
            if logy: ROOT.gPad.SetLogz()
            hist.GetXaxis().SetTitleOffset(1.15); hist.GetXaxis().SetLabelSize(0.05); hist.GetXaxis().SetTitleSize(0.05)
            hist.GetYaxis().SetTitleOffset(1.5);  hist.GetYaxis().SetLabelSize(0.05); hist.GetYaxis().SetTitleSize(0.05)
            hist.GetZaxis().SetTitleOffset(1.5);  hist.GetZaxis().SetLabelSize(0.05); hist.GetZaxis().SetTitleSize(0.05)
            hist.GetXaxis().SetNdivisions(505)
            if 'lego' in datastyle.lower():
                hist.GetZaxis().SetTitleOffset(1.4)
            if len(titles) > 0:
                hist.SetTitle(titles[hist_index])

            if datastyle != 'pe': hist.Draw(datastyle)
            else: hist.Draw('colz')
            
            CMS_lumi.extraText = extraText
            CMS_lumi.CMS_lumi(thisPad, year, 11, sim=False if 'data' in name.lower() else True)
        
        # Otherwise it's a TH1 hopefully
        else:
            hist.SetLineColorAlpha(ROOT.kBlack, 0 if dataOff else 1)
            if 'pe' in datastyle.lower():
                hist.SetMarkerColorAlpha(ROOT.kBlack,0 if dataOff else 1)
                hist.SetMarkerStyle(8)
            if 'hist' in datastyle.lower():
                hist.SetFillColorAlpha(0,0)

            # If there are no backgrounds, only plot the data (semilog if desired)
            if len(bkglist) == 0:
                hist.SetMaximum(1.13*hist.GetMaximum())
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.SetTitleOffset(1.1)
                hist.Draw(datastyle)
                CMS_lumi.CMS_lumi(thisPad, year, 11)
            
            # Otherwise...
            else:
                # Create some subpads, a legend, a stack, and a total bkg hist that we'll use for the error bars
                if not dataOff:
                    mains.append(ROOT.TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.35, 1, 1))
                    subs.append(ROOT.TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 1, 0.35))

                else:
                    mains.append(ROOT.TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.1, 1, 1))
                    subs.append(ROOT.TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 0, 0))

                if len(signals) == 0:
                    nsignals = 0
                elif addSignals:
                    nsignals = 1
                else:
                    nsignals = len(signals[0])
                legend_topY = 0.73-0.03*(min(len(bkglist[0]),6)+nsignals+1)
                # legend_bottomY = 0.2+0.02*(len(bkglist[0])+nsignals+1)

                legends.append(ROOT.TLegend(0.65,legend_topY,0.90,0.88))
                if not dataOff: legends[hist_index].AddEntry(hist,dataName,datastyle)

                stacks.append(ROOT.THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
                if totalBkg == None:
                    tot_hist = hist.Clone(hist.GetName()+'_tot')
                    tot_hist.Reset()
                else:
                    tot_hist = totalBkg[hist_index]

                tot_hist.SetTitle(hist.GetName()+'_tot')
                tot_hist.SetMarkerStyle(0)
                tot_hists.append(tot_hist)
                tot_hists_err.append(tot_hist.Clone())
                tot_hists[hist_index].SetLineColor(ROOT.kBlack)
                tot_hists_err[hist_index].SetLineColor(ROOT.kBlack)
                tot_hists_err[hist_index].SetLineWidth(0)
                tot_hists_err[hist_index].SetFillColor(ROOT.kBlack)
                tot_hists_err[hist_index].SetFillStyle(3354)
                legends[hist_index].AddEntry(tot_hists_err[hist_index],'Total bkg unc.','F')

                # Set margins and make these two pads primitives of the division, thisPad
                mains[hist_index].SetBottomMargin(0.04)
                mains[hist_index].SetLeftMargin(0.17)
                mains[hist_index].SetRightMargin(0.05)
                mains[hist_index].SetTopMargin(0.1)

                subs[hist_index].SetLeftMargin(0.17)
                subs[hist_index].SetRightMargin(0.05)
                subs[hist_index].SetTopMargin(0)
                subs[hist_index].SetBottomMargin(0.35)
                mains[hist_index].Draw()
                subs[hist_index].Draw()
                
                # Build the stack
                legend_info = collections.OrderedDict()
                for bkg in bkglist[hist_index]:     # Won't loop if bkglist is empty
                    if totalBkg == None:
                        tot_hists[hist_index].Add(bkg)
                    if logy:
                        bkg.SetMinimum(1e-3)

                    stacks[hist_index].Add(bkg)
                    legend_info[bkg.GetTitle()] = bkg

                # Deal with legend which needs ordering reversed from stack build
                legend_duplicates = []
                for bname in reversed(legend_info.keys()):
                    if bname not in legend_duplicates:
                        legends[hist_index].AddEntry(legend_info[bname],bname,'f')
                        legend_duplicates.append(bname)
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()

                # Set y max of all hists to be the same to accommodate the tallest
                histList = [stacks[hist_index],tot_hists[hist_index],hist]

                yMax = get_hist_maximum(histList)
                for h in histList:
                    h.SetMaximum(yMax*(2.5-legend_topY+0.03))
                    if logy == True:
                        h.SetMaximum(yMax*10**(2.5-legend_topY+0.1))

                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.SetTitleOffset(1.15,"xy")
                hist.GetYaxis().SetTitleOffset(1.04)
                hist.GetYaxis().SetLabelSize(0.07)
                hist.GetYaxis().SetTitleSize(0.09)
                hist.GetXaxis().SetLabelSize(0.07)
                hist.GetXaxis().SetTitleSize(0.09)
                hist.GetXaxis().SetLabelOffset(0.05)
                if logy == True:
                    hist.SetMinimum(1e-3)
                
                hist.GetYaxis().SetNdivisions(508)
                hist.Draw(datastyle+' X0')

                if logy == True:stacks[hist_index].SetMinimum(1e-3) 
                
                stacks[hist_index].Draw('same hist') # need to draw twice because the axis doesn't exist for modification until drawing
                try:
                    stacks[hist_index].GetYaxis().SetNdivisions(508)
                except:
                    stacks[hist_index].GetYaxis().SetNdivisions(8,5,0)
                stacks[hist_index].Draw('same hist')

                # Do the signals
                sigs_to_plot = []
                if len(signals) > 0: 
                    # Can add together for total signal
                    if addSignals:
                        totsig = signals[hist_index][0].Clone()
                        for isig in range(1,len(signals[hist_index])):
                            totsig.Add(signals[hist_index][isig])
                        sigs_to_plot = [totsig]
                    # or treat separately
                    else:
                        for sig in signals[hist_index]:
                            sigs_to_plot.append(sig)

                    # Plot either way
                    tot_sigs.append(sigs_to_plot)
                    for isig,sig in enumerate(sigs_to_plot):
                        sig.SetLineWidth(2)
                        if logy == True:
                            sig.SetMinimum(1e-3)

                        legends[hist_index].AddEntry(sig,sig.GetTitle(),'L')
                        sig.Draw('hist same')

                # Draw total hist and error
                if logy: 
                    tot_hists[hist_index].SetMinimum(1e-3)
                    tot_hists_err[hist_index].SetMinimum(1e-3)
                tot_hists[hist_index].Draw('hist same')
                tot_hists_err[hist_index].Draw('e2 same')

                legends[hist_index].SetBorderSize(0)
                legends[hist_index].Draw()

                if not dataOff:
                    hist.Draw(datastyle+' X0 same')

                ROOT.gPad.RedrawAxis()

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull
                pulls.append(MakePullPlot(hist,tot_hists[hist_index]))
                pulls[hist_index].SetFillColor(ROOT.kBlue)
                pulls[hist_index].SetTitle(";"+hist.GetXaxis().GetTitle()+";(Data-Bkg)/#sigma")
                pulls[hist_index].SetStats(0)

                pulls[hist_index].GetYaxis().SetRangeUser(-2.9,2.9)
                pulls[hist_index].GetYaxis().SetTitleOffset(0.4)                             
                pulls[hist_index].GetYaxis().SetLabelSize(0.13)
                pulls[hist_index].GetYaxis().SetTitleSize(0.12)
                pulls[hist_index].GetYaxis().SetNdivisions(306)

                pulls[hist_index].GetXaxis().SetLabelSize(0.13)
                pulls[hist_index].GetXaxis().SetTitleSize(0.15)
                pulls[hist_index].GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
                pulls[hist_index].GetYaxis().SetTitle("(Data-Bkg)/#sigma")
                pulls[hist_index].Draw('hist')

                if logy == True:
                    mains[hist_index].SetLogy()

                CMS_lumi.extraText = extraText
                CMS_lumi.cmsTextSize = 0.9
                CMS_lumi.cmsTextOffset = 2
                CMS_lumi.lumiTextSize = 0.9
                
                CMS_lumi.CMS_lumi(mains[hist_index], year, 11)
                mains[hist_index].cd()
                lumiE = ROOT.TLatex()
                lumiE.SetNDC()
                lumiE.SetTextAngle(0)
                lumiE.SetTextColor(ROOT.kBlack)
                lumiE.SetTextFont(42)
                lumiE.SetTextAlign(31) 
                lumiE.SetTextSize(0.7*0.1)
                lumiE.DrawLatex(1-0.05,1-0.1+0.2*0.1,"137 fb^{-1} (13 TeV)")
                
                if isinstance(subtitles,list) and len(subtitles) > 0:
                    subtitle = ROOT.TLatex()
                    subtitle.SetNDC()
                    subtitle.SetTextAngle(0)
                    subtitle.SetTextColor(ROOT.kBlack)
                    subtitle.SetTextFont(42)
                    subtitle.SetTextAlign(12) 
                    subtitle.SetTextSize(0.06)
                    # print (subtitles[hist_index])
                    subtitle_string = '%s < %s < %s %s'%(subtitles[hist_index].split('-')[0], sliceVar.split(' ')[0], subtitles[hist_index].split('-')[1], 'GeV')
                    subtitle.DrawLatex(0.208,0.74,subtitle_string)

    if rootfile:
        myCan.Print(tag+'/'+name+'.root','root')
    else:
        myCan.Print(tag+'/'+name+'.pdf','pdf')
        myCan.Print(tag+'/'+name+'.png','png')

def MakeSystematicPlots(configObj):
    '''Make plots of the systematic shape variations of each process based on those
    processes and systematic shapes specified in the config. Shapes are presented 
    as projections onto 1D axis where no selection has been made on the axis not
    being plotted. Plots are saved to UncertPlots/.
    '''
    c = ROOT.TCanvas('c','c',800,700)

    for pair, group in configObj.df.groupby(['process','region']):
        p,r = pair
        nominal_full = configObj.organizedHists.Get(process=p,region=r,systematic='')
        binning = configObj.BinningLookup(nominal_full.GetName())
        for axis in ['X','Y']:
            nominal = getattr(nominal_full,'Projection'+axis)('%s_%s_%s_%s'%(p,r,'nom','proj'+axis))
            for s, _ in group.loc[group.variation.ne('nominal')].groupby('variation'):
                up = getattr(configObj.organizedHists.Get(process=p,region=r,systematic=s+'Up'),'Projection'+axis)('%s_%s_%s_%s'%(p,r,s+'Up','proj'+axis))
                down = getattr(configObj.organizedHists.Get(process=p,region=r,systematic=s+'Down'),'Projection'+axis)('%s_%s_%s_%s'%(p,r,s+'Down','proj'+axis))

                c.cd()
                nominal.SetLineColor(ROOT.kBlack)
                nominal.SetFillColor(ROOT.kYellow-9)
                up.SetLineColor(ROOT.kRed)
                down.SetLineColor(ROOT.kBlue)

                up.SetLineStyle(9)
                down.SetLineStyle(9)
                up.SetLineWidth(2)
                down.SetLineWidth(2)

                nominal,up,down = set_hist_maximums([nominal,up,down])
                nominal.SetXTitle(getattr(binning,axis.lower()+'title'))

                nominal.SetTitle('')
                nominal.GetXaxis().SetTitleOffset(1.0)
                nominal.GetXaxis().SetTitleSize(0.05)
                c.SetRightMargin(0.16)

                nominal.Draw('hist')
                up.Draw('same hist')
                down.Draw('same hist')

                c.Print(configObj.projPath+'/UncertPlots/Uncertainty_%s_%s_%s_%s.png'%(p,r,s,'proj'+axis),'png')

def _get_signal_strength(fd_file):
    tree_fit = fd_file.Get('tree_fit_s')
    tree_fit.GetEntry(0)
    return tree_fit.r

def reducedCorrMatrixHist(fit_result,varsOfInterest=[]):
    ROOT.gStyle.SetOptStat(0)
    # ROOT.gStyle.SetPaintTextFormat('.3f')
    CM = fit_result.correlationMatrix()
    finalPars = fit_result.floatParsFinal()

    nParams = CM.GetNcols()
    finalParamsDict = {}
    for cm_index in range(nParams):
        if varsOfInterest == []:
            if 'Fail_' not in finalPars.at(cm_index).GetName():
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index
        else:
            if finalPars.at(cm_index).GetName() in varsOfInterest:
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index

    nFinalParams = len(finalParamsDict.keys())
    out = ROOT.TH2D('correlation_matrix','correlation_matrix',nFinalParams,0,nFinalParams,nFinalParams,0,nFinalParams)
    out_txt = open('correlation_matrix.txt','w')

    for out_x_index, paramXName in enumerate(sorted(finalParamsDict.keys())):
        cm_index_x = finalParamsDict[paramXName]
        for out_y_index, paramYName in enumerate(sorted(finalParamsDict.keys())):
            cm_index_y = finalParamsDict[paramYName]
            if cm_index_x > cm_index_y:
                out_txt.write('%s:%s = %s\n'%(paramXName,paramYName,CM[cm_index_x][cm_index_y]))
            out.Fill(out_x_index+0.5,out_y_index+0.5,CM[cm_index_x][cm_index_y])

        out.GetXaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
        out.GetYaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
    out.SetMinimum(-1)
    out.SetMaximum(+1)

    return out

def ColorCodeSortedIndices(colors):
    possible_colors = []
    for c in colors:
        if c not in possible_colors:
            possible_colors.append(c)

    new_color_order = []
    for c in possible_colors:
        for idx in [idx for idx,color in enumerate(colors) if color == c]:
            if idx not in new_color_order:
                new_color_order.append(idx)
    return new_color_order

def MakePullPlot( DATA,BKG):
    BKGUP, BKGDOWN = MakeUpDown(BKG)
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
        if FSerr != None:
            sigma = math.sqrt(FSerr*FSerr + BKGerr*BKGerr)
        else:
            sigma = math.sqrt(BKGerr*BKGerr)
        if FScont == 0.0:
            pull.SetBinContent(ibin, 0.0 )  
        else:
            if sigma != 0 :
                pullcont = (pull.GetBinContent(ibin))/sigma
                pull.SetBinContent(ibin, pullcont)
            else :
                pull.SetBinContent(ibin, 0.0 )
    return pull

def MakeUpDown(hist):
    hist_up = hist.Clone(hist.GetName()+'_up')
    hist_down = hist.Clone(hist.GetName()+'_down')

    for xbin in range(1,hist.GetNbinsX()+1):
        errup = hist.GetBinErrorUp(xbin)
        errdown = hist.GetBinErrorLow(xbin)
        nom = hist.GetBinContent(xbin)

        hist_up.SetBinContent(xbin,nom+errup)
        hist_down.SetBinContent(xbin,nom-errdown)

    return hist_up,hist_down

def ReorderHists(histlist):
    if len(histlist) != 6:
        raise IndexError('reorderHists() only built to rearrange list of six hists from 2x3 to 3x2')

    outlist = []
    outlist.append(histlist[0])
    outlist.append(histlist[3])
    outlist.append(histlist[1])
    outlist.append(histlist[4])
    outlist.append(histlist[2])
    outlist.append(histlist[5])

    return outlist

def NuisPulls():
    diffnuis_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnostics.root --abs -g nuisance_pulls.root'
    executeCmd(diffnuis_cmd)
    # Make a PDF of the nuisance_pulls.root
    if os.path.exists('nuisance_pulls.root'):
        nuis_file = ROOT.TFile.Open('nuisance_pulls.root')
        nuis_can = nuis_file.Get('nuisances')
        nuis_can.Print('nuisance_pulls.pdf','pdf')
        nuis_file.Close()

def _getGoodFitResults(tfile):
    successful_fits = []
    for fittag in ['b','s']:
        if 'fit_'+fittag not in [k.GetName() for k in tfile.GetListOfKeys()]:
            warnings.warn('Unable to find result fit_%s...'%fittag,RuntimeWarning)
        else:
            successful_fits.append(fittag)

def SavePostFitParametricFuncVals():
    fit_result_file = ROOT.TFile.Open('fitDiagnostics.root')
    goodFitTags = _getGoodFitResults(fit_result_file)
    for fittag in goodFitTags:
        coeffs_final = fit_result_file.Get('fit_'+fittag).floatParsFinal()
        all_par_names = [] # get names of anything matching *_par<N>
        for i in range(coeffs_final.getSize):
            name = coeffs_final.at(i).GetName()
            if name.split('_')[-1].startswith('par'):
                all_par_names.append(name)
        all_par_names.sort()
        # Get unique prefixes to _par<N>
        all_obj_names = set(['_'.join(name.split('_')[:-1]) for name in all_par_names])

        for obj_name in all_obj_names:
            with open('{tag}/rpf_params_%s_fit%s.txt'%(self.tag,obj_name,fittag),'w') as param_out:
                for par_name in [p for p in all_par_names if p.startswith(obj_name)]:
                    var = coeffs_final.find(par_name)
                    param_out.write('%s: %s +/- %s\n'%(par_name, var.getValV(), var.getError()))

def GenPostFitShapes():
    fit_result_file = ROOT.TFile.Open('fitDiagnostics.root')
    goodFitTags = _getGoodFitResults(fit_result_file)
    for t in goodFitTags:
        shapes_cmd = 'PostFit2DShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -o postfitshapes_{0}.root -f fitDiagnostics.root:fit_{0} --postfit --sampling --samples 100 --print 2> PostFitShapes2D_stderr_{0}.txt'.format(t)
        executeCmd(shapes_cmd)

    if 'b' in goodFitTags:
        fit_result = fit_result_file.Get("fit_b")
        if hasattr(fit_result,'correlationMatrix'):
            corrMtrx = reducedCorrMatrixHist(fit_result)
            corrMtrxCan = ROOT.TCanvas('c','c',1400,1000)
            corrMtrxCan.cd()
            corrMtrxCan.SetBottomMargin(0.22)
            corrMtrxCan.SetLeftMargin(0.17)
            corrMtrxCan.SetTopMargin(0.06)

            corrMtrx.GetXaxis().SetLabelSize(0.01)
            corrMtrx.GetYaxis().SetLabelSize(0.01)
            corrMtrx.Draw('colz text')
            corrMtrxCan.Print('correlation_matrix.png','png')
        else:
            warnings.warn('Not able to produce correlation matrix.',RuntimeWarning)