from collections import OrderedDict
import ROOT, os, warnings, pandas, collections, math
from TwoDAlphabet.helpers import get_hist_maximum, set_hist_maximums, execute_cmd
from TwoDAlphabet.binning import stitch_hists_in_x, convert_to_events_per_unit, get_min_bin_width
from TwoDAlphabet.ext import tdrstyle, CMS_lumi

class Plotter(object):
    '''Class to manage output distributions, manipulate them, and provide access to plotting
    standard groups of distributions.

    Attributes:
        fit_tag (str): Either 's' or 'b'.
        fit_results (RooFitResult): The RooFitResult corresponding to the fit_tag.
        signal_strength (float): The post-fit signal strength.
        twoD (TwoDAlphabet): TwoDAlphabet object storing various meta information needed for access.
        yaxis1D_title (str): Title for "counts" axis of 1D plots. Defaults to 'Events / bin' but can change if plotting events per unit.
        df (pandas.DataFrame): DataFrame organizing the post-fit and pre-fit plots and their 1D projections.
        dir (str): Directory path to save final images.
        slices (dict): Stores edges to slice "x" and "y" axes. 
        root_out (ROOT.TFile): File storing all histograms that are made.
    '''
    def __init__(self,twoD,fittag,loadExisting=False):
        '''Constructor.

        Args:
            twoD (TwoDAlphabet): Object with meta information about the run.
            fittag (str): Either 's' or 'b'.
            loadExisting (bool, optional): Flag to load existing projections instead of remaking everything. Defaults to False.
        '''
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
        self.signal_strength = _get_signal_strength(fd_file)
        self.twoD = twoD
        self.yaxis1D_title = 'Events / bin'
        self.df = pandas.DataFrame(columns=columns)
        self.dir = 'plots_fit_{f}'.format(f=self.fittag)
        self.slices = {'x': {}, 'y': {}}
        self.root_out = None

        if not loadExisting:
            self._make()
        else:
            self._load()

    def __del__(self):
        '''On deletion, save DataFrame to csv and close ROOT file.'''
        self.df.to_csv('%s/df.csv'%self.dir)
        self.root_out.Close()

    class _entryTracker(object):
        '''Subclass to track entries for the output ROOT file and DataFrame.
        The entry itself is a row in the DataFrame with columns holding meta
        information and histogram objects. The meta information describes the
        process, process type, a "pretty" process title, and the selectioin region.
        Since all histogram objects share this information, the various interpretations
        (1D vs 2D, projections, pre-fit vs post-fit, etc) are stored in the same row.
        These histograms can be added to the row via the `track` method. Only the string
        storing the name is tracked and the actual histogram is saved to the output ROOT file
        storing all histograms for all rows.
        '''
        def __init__(self,plotter,p,r,t,title):
            '''Constructor.

            Args:
                plotter (Plotter): Parent plotter object for access to open ROOT file.
                p (str): Process name of entry.
                r (str): Region name of entry.
                t (str): Process type of entry (ex. "BKG").
                title (str): Pretty process title for the legend.
            '''
            self.plotter = plotter
            self.entry = {'process': p, 'region': r, 'type': t, 'title': title}

        def track(self,k,h):
            '''Add a key to the `self.entry` dict to track a histogram of name `h.GetName()`.
            Also write the histogram to output ROOT file.

            Args:
                k (str): Key name (ex. 'prefit_2D').
                h (ROOT.TH1): ROOT histogram to track.
            '''
            self.entry[k] = h.GetName()
            self.plotter.root_out.WriteTObject(h,h.GetName())

    def _load(self):
        '''Open pickled DataFrame and output ROOT file
        and reference with `self.df` and `self.root_out` attributes.'''
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name)
        self.df = pandas.read_csv('%s/df.csv'%self.dir)

    def _make(self):
        '''Make the DataFrame and output ROOT file from scratch
        and reference with `self.df` and `self.root_out` attributes.
        
        Loops over all regions and processes from the pre-fit and post-fit shapes
        and tracks/constructs the 2D histograms and six projections (three each for "x" and "y"). 
        '''
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name,'RECREATE')

        shapes_file = ROOT.TFile.Open('postfitshapes_%s.root'%self.fittag)
        loc_base = '{r}_{c}_{t}/{p}'

        for region in self.twoD.GetRegions():
            binning,_ = self.twoD.GetBinningFor(region)
            if self.twoD.config.options.blindedPlots != None:
                if region in self.twoD.config.options.blindedPlots:
                    blinding = [1]
                else: blinding = []
            else: blinding = []

            self.slices['x'][region] = {'vals': binning.xSlices,'idxs':binning.xSliceIdx}
            self.slices['y'][region] = {'vals': binning.ySlices,'idxs':binning.ySliceIdx}
            
            for process in self.twoD.GetProcesses()+['TotalBkg']:
                if process != 'TotalBkg':
                    color = self.twoD.GetProcessColor(process)
                    proc_type = self.twoD.GetProcessType(process)
                    proc_title = self.twoD.GetProcessTitle(process)
                else:
                    color = ROOT.kBlack
                    proc_type = 'TOTAL'
                    proc_title = 'TotalBkg'

                entry = self._entryTracker(self, process, region, proc_type, proc_title)
                for time in ['prefit','postfit']:
                    # 2D distributions first
                    out2d_name = '%s_%s_%s_2D'%(process,region,time)
                    low  = shapes_file.Get(loc_base.format(r=region, c='LOW', t=time, p=process))
                    sig  = shapes_file.Get(loc_base.format(r=region, c='SIG', t=time, p=process))
                    high = shapes_file.Get(loc_base.format(r=region, c='HIGH',t=time, p=process))

                    full = stitch_hists_in_x(out2d_name, binning, [low,sig,high], blinded=blinding if process == 'data_obs' else [])
                    full.SetMinimum(0)
                    full.SetTitle('%s, %s, %s'%(process,region,time))

                    entry.track('%s_2D'%time, full)

                    # Now do projections using the 2D
                    out_proj_name = '{p}_{r}_{t}_proj{x}{i}'
                    for proj in ['X','Y']:
                        slices = self.slices['x' if proj == 'Y' else 'y'][region]

                        for islice in range(3):
                            hname = out_proj_name.format(p=process,r=region,t=time,x=proj,i=islice)
                            start,stop = _get_start_stop(islice,slices['idxs'])
                            
                            hslice = getattr(full,'Projection'+proj)(hname,start,stop,'e')
                            hslice_title = '%s, %s, %s, %s-%s'%(proc_title,region,time,slices['vals'][islice],slices['vals'][islice+1])
                            hslice = self._format_1Dhist(
                                hslice, hslice_title,
                                binning.xtitle if proj == 'X' else binning.ytitle,
                                self.yaxis1D_title,
                                color, proc_type)

                            col_name = '{t}_proj{x}{i}'.format(t=time,x=proj.lower(),i=islice)
                            entry.track(col_name, hslice)

                self.df = self.df.append(entry.entry,ignore_index=True)

        shapes_file.Close()
        # self.root_out.Close()
        # self.root_out = ROOT.TFile.Open(root_out_name)

    def Get(self,hname):
        '''Get a histogram by name from the master ROOT file.
        Does quick check for if the histogram is saved.
        
        Args:
            hname (str): Histogram name.

        Raises:
            LookupError: If histogram cannot be found.
        '''
        if hname not in [k.GetName() for k in self.root_out.GetListOfKeys()]:
            raise LookupError('Histogram %s not found in %s'%(hname,self.root_out.GetName()))
        return self.root_out.Get(hname).Clone()

    def _format_1Dhist(self, hslice, title, xtitle, ytitle, color, proc_type):
        '''Perform some basic formatting of a 1D histogram so that the ROOT.TH1
        already has some of the meta information set like line and fill colors.
        Also renormalizes bins to the minimum bin width if the y-axis is requested
        to be plotted as events per unit rather than events per bin.

        Args:
            hslice (ROOT.TH1): Histogram.
            title (str): Title for the output histogram.
            xtitle (str): X-axis title for the output histogram.
            ytitle (str): Y-axis title for the output histogram. Will be modified if requesting to plot events/unit.
            color (int): ROOT color code.
            proc_type (str): Process type. Either 'BKG', 'SIGNAL', or 'DATA'.

        Raises:
            NameError: If proc_type is not 'BKG', 'SIGNAL', or 'DATA'.

        Returns:
            TH1: Formatted histogram.
        '''
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
            hslice.SetLineColorAlpha(0,0)
        elif proc_type == 'SIGNAL' or proc_type == 'TOTAL':
            hslice.SetLineColor(color)
        elif proc_type == 'DATA':
            hslice.SetLineColor(color)
            hslice.SetMarkerColor(color)
        else:
            raise NameError('Process type "%s" is not supported.'%proc_type)

        return hslice

    def _order_df_on_proc_list(self,df,proc_type,proclist=[],alphaBottom=True):
        '''Re-order input dataframe (`df`) based on the ordered list of process names (`proclist`).
        Useful for pre-ordering process before trying to construct a THStack.

        Args:
            df (pandas.DataFrame): Input DataFrame to manipulate.
            proc_type (str): Process type.  Either 'BKG', 'SIGNAL', or 'DATA'.
            proclist (list, optional): Ordered list of processes to order by. Defaults to [] in which case order is determined based on alphaBottom.
            alphaBottom (bool, optional): Only matters if proclist == []. Defaults to True in which case parametric Alphabet objects included first (thus, on the bottom of the THStack).

        Returns:
            pandas.DataFrame: Ordered DataFrame.
        '''
        if proclist == []:
            if alphaBottom:
                process_order = self.twoD.GetProcesses(ptype=proc_type,onlyNonConfig=True) + self.twoD.GetProcesses(ptype=proc_type,includeNonConfig=False)
            else:
                process_order = self.twoD.GetProcesses(ptype=proc_type,includeNonConfig=False) + self.twoD.GetProcesses(ptype=proc_type,onlyNonConfig=True)
        else:
            process_order = proclist

        process_order_df = pandas.DataFrame({'process':process_order})
        return process_order_df.merge(df,on='process',how='outer')

    def plot_2D_distributions(self):
        '''Take the saved 2D distributions and plot them together on sub-pads
        based on process and region groupings.

        Plots are grouped based on process and then the regions and pre-fit/post-fit
        plots share the same canvas as sub-pads.

        Returns:
            None
        '''
        for process, group in self.df.groupby('process'):
            pads = []
            for region,subgroup in group.groupby('region'):           
                out_file_name = 'plots_fit_{f}/%s_{p}_{r}_2D'.format(f=self.fittag,p=process,r=region)
                pads.append(self.MakePad2D(outname=out_file_name%('prefit'), hist=self.Get(subgroup.prefit_2D.iloc[0]),
                               year=self.twoD.config.options.year, savePDF=False, savePNG=False, saveROOT=True))
                pads.append(self.MakePad2D(outname=out_file_name%('postfit'), hist=self.Get(subgroup.postfit_2D.iloc[0]),
                               year=self.twoD.config.options.year, savePDF=False, savePNG=False, saveROOT=True))

                self.MakeCan('plots_fit_{f}/{p}_{r}_2D'.format(f=self.fittag,p=process,r=region), pads=pads)


    def plot_projections(self,logyFlag):
        '''Plot comparisons of data and the post-fit background model and signal
        using the 1D projections. Canvases are grouped based on projection axis.
        The canvas rows are separate selection regions while the columns 
        are the different slices of the un-plotted axis.

        Args:
            logyFlag (bool): If True, set the y-axis to be log scale.

        Returns:
            None
        '''
        
        for proj in ['postfit_projx','postfit_projy']:
            pads = OrderedDict()
            for region, group in self.df.groupby('region'):
                binning,_ = self.twoD.GetBinningFor(region)
                ordered_bkgs = self._order_df_on_proc_list(
                                    group.loc[group.type.eq('BKG')],
                                    proc_type='BKG',
                                    alphaBottom=(not logyFlag))
                ordered_signals = self._order_df_on_proc_list(
                                    group.loc[group.type.eq('SIGNAL')],
                                    proc_type='SIGNAL',
                                    alphaBottom=(not logyFlag))

                for islice in range(3):
                    projn = proj+str(islice)
                    this_data = self.Get( group.loc[group.type.eq('DATA')][projn].iloc[0] )
                    this_totalbkg = self.Get( group.loc[group.type.eq('TOTAL')][projn].iloc[0] )
                    these_bkgs = [self.Get(hname) for hname in ordered_bkgs[projn].to_list()]
                    if self.twoD.config.options.plotPrefitSigInFitB and self.fittag == 'b':
                        these_signals = [self.Get(hname) for hname in ordered_signals[projn.replace('postfit','prefit')].to_list()]
                    else:
                        these_signals = [self.Get(hname) for hname in ordered_signals[projn].to_list()]

                    slice_edges = (
                        self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice],
                        binning.xtitle if 'y' in proj else binning.ytitle,
                        self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice+1],
                        'GeV'
                    )
                    slice_str = '%s < %s < %s %s'%slice_edges

                    out_file_name = 'plots_fit_{f}/{projn}_{reg}{logy}'.format(
                                        f=self.fittag, projn=projn, reg=region,
                                        logy='' if logyFlag == False else '_logy')
                    print ('DEBUG: %s' %out_file_name)
                    pads['%s_%s'%(projn,region)] = self.MakePad1D(out_file_name, this_data, these_bkgs, these_signals,
                                   subtitle=slice_str, totalBkg=this_totalbkg,
                                   logyFlag=logyFlag, year=self.twoD.config.options.year,
                                   extraText='Preliminary', savePDF=False, savePNG=False, saveROOT=True)
            
            out_can_name = 'plots_fit_{f}/{proj}{logy}'.format(
                                        f=self.fittag, proj=proj,
                                        logy='' if logyFlag == False else '_logy')
            print ('DEBUG: %s'%out_can_name)
            # TODO: Change order by unpacking pads with an ordered list
            self.MakeCan(out_can_name, pads.values())

    def plot_pre_vs_post(self):
        '''Make comparisons for each background process of pre and post fit projections.
        '''
        for proj in ['projx','projy']:
            for process, group in self.df[~self.df.process.isin(['data_obs','TotalBkg'])].groupby('process'):
                pads = OrderedDict()
                for region, subgroup in group.groupby('region'):
                    binning,_ = self.twoD.GetBinningFor(region)
                    for islice in range(3):
                        projn = proj+str(islice)

                        post = self.Get( subgroup['postfit_'+projn].iloc[0] )
                        post.SetLineColor(ROOT.kBlack)
                        post.SetTitle('          Postfit, '+process)

                        pre = self.Get( subgroup['prefit_'+projn].iloc[0] )
                        pre.SetTitle('Prefit, '+process)

                        slice_edges = (
                            self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice],
                            binning.xtitle if 'y' in proj else binning.ytitle,
                            self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice+1],
                            'GeV'
                        )
                        slice_str = '%s < %s < %s %s'%slice_edges

                        key = '{p}_{projn}_{reg}'.format(p=process,projn=projn, reg=region)
                        print ('DEBUG: %s'%key)
                        pads[key] = self.MakePad1D(
                            'plots_fit_{f}/{k}'.format(f=self.fittag, k=key),
                            post, [pre], totalBkg=pre, subtitle=slice_str, savePNG=False, savePDF=False, saveROOT=True,
                            datastyle='histe', year=self.twoD.config.options.year, extraText='Preliminary'
                        )

                print ('DEBUG: %s'%'plots_fit_{f}/{p}_{proj}'.format(f=self.fittag, p=process,proj=proj))
                self.MakeCan(
                    'plots_fit_{f}/{p}_{proj}'.format(f=self.fittag, p=process,proj=proj),
                    pads.values()
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

    def _save_pad_generic(self, outname, pad, saveROOT, savePDF, savePNG):
        if saveROOT:
            self.root_out.WriteTObject(pad,outname)
        if savePDF:
            pad.Print(outname+'.pdf','pdf')
        if savePNG:
            pad.Print(outname+'.png','png')

    def _make_pad_gen(self, name):
        tdrstyle.setTDRStyle()
        ROOT.gStyle.SetLegendFont(42)
        ROOT.gStyle.SetTitleBorderSize(0)
        ROOT.gStyle.SetTitleAlign(33)
        ROOT.gStyle.SetTitleX(.77)

        pad = ROOT.TCanvas(name, name, 800, 700)
        pad.cd(); pad.SetRightMargin(0.0); pad.SetTopMargin(0.0); pad.SetBottomMargin(0.0)
        return pad

    def MakePad2D(self, outname, hist, style='lego', logzFlag=False,
                  saveROOT=True, savePDF=False, savePNG=False, year=1, extraText='Preliminary'):
        '''Make a pad holding a 2D plot with standardized formatting conventions.

        Args:
            outname (str): Output file path name.
            hist (TH2): Histogram to draw on the pad.
            style (str, optional): ROOT drawing style. Defaults to 'lego'.
            logzFlag (bool, optional): Make log z-axis. Defaults to False.
            saveROOT (bool, optional): Save to master ROOT file. Defaults to True.
            savePDF (bool, optional): Save to PDF. Defaults to False.
            savePNG (bool, optional): Save to PNG. Defaults to False.
            year (int, optional): Luminosity formatting. Options are 16, 17, 18, 1 (full Run 2), 2 (16+17+18). Defaults to 1.
            extraText (str, optional): Prepended to the CMS subtext. Defaults to 'Preliminary'.

        Returns:
            [type]: [description]
        '''
        pad = self._make_pad_gen(outname)
        pad.SetLeftMargin(0.15)
        pad.SetRightMargin(0.2)
        pad.SetBottomMargin(0.12)
        pad.SetTopMargin(0.1)
        if logzFlag: pad.SetLogz()

        hist.GetXaxis().SetTitleOffset(1.15); hist.GetXaxis().SetLabelSize(0.05); hist.GetXaxis().SetTitleSize(0.05)
        hist.GetYaxis().SetTitleOffset(1.5);  hist.GetYaxis().SetLabelSize(0.05); hist.GetYaxis().SetTitleSize(0.05)
        hist.GetZaxis().SetTitleOffset(1.5);  hist.GetZaxis().SetLabelSize(0.05); hist.GetZaxis().SetTitleSize(0.05)
        hist.GetXaxis().SetNdivisions(505)
        
        if 'lego' in style.lower():
            hist.GetZaxis().SetTitleOffset(1.4)

        hist.Draw(style)
        
        CMS_lumi.extraText = extraText
        CMS_lumi.CMS_lumi(pad, year, 11, sim=False if 'data' in hist.GetName().lower() else True)

        self._save_pad_generic(outname, pad, saveROOT, savePDF, savePNG)

        return pad

    def MakePad1D(self, outname, data, bkgs=[], signals=[], title='', subtitle='',
                totalBkg=None, logyFlag=False, saveROOT=True, savePDF=False, savePNG=False,
                dataOff=False, datastyle='pe X0', year=1, addSignals=True, extraText='Preliminary'):
        '''Make a pad holding a 1D plot with standardized formatting conventions.

        Args:
            outname (str): Output file path name.
            data (TH1): Data histogram.
            bkgs ([TH1]): List of background histograms (will be stacked).
            signals ([TH1]): List of signal histograms.
            title (str, optional): Title of plot. Only applicable if bkgs is empty. Defaults to ''.
            subtitle (str, optional): Subtitle text for physics information (like slice ranges). Defaults to ''.
            totalBkg (TH1, optional): Total background estimate from fit. Used to get total background uncertianty. Defaults to None.
            logyFlag (bool, optional): Make log y-axis. Defaults to False.
            saveROOT (bool, optional): Save to master ROOT file. Defaults to False.
            savePDF (bool, optional): Save to PDF. Defaults to True.
            savePNG (bool, optional): Save to PNG. Defaults to True.
            dataOff (bool, optional): Turn off the data from plotting. Defaults to False.
            datastyle (str, optional): ROOT drawing style for the data. Defaults to 'pe X0'.
            year (int, optional): Luminosity formatting. Options are 16, 17, 18, 1 (full Run 2), 2 (16+17+18). Defaults to 1.
            addSignals (bool, optional): If True, multiple signals will be added together and plotted as one. If False, signals are plotted individually. Defaults to True.
            extraText (str, optional): Prepended to the CMS subtext. Defaults to 'Preliminary'.

        Returns:
            ROOT.TPad: Output pad.
        '''
        def _format_data():
            data.SetBinErrorOption(ROOT.TH1.kPoisson)
            data.SetLineColorAlpha(ROOT.kBlack, 0 if dataOff else 1)
            if 'pe' in datastyle.lower():
                data.SetMarkerColorAlpha(ROOT.kBlack,0 if dataOff else 1)
                data.SetMarkerStyle(8)
            if 'hist' in datastyle.lower():
                data.SetFillColorAlpha(0,0)

            data.SetTitleOffset(1.15,"xy")
            data.GetYaxis().SetTitleOffset(1.04)
            data.GetYaxis().SetLabelSize(0.07)
            data.GetYaxis().SetTitleSize(0.09)
            data.GetXaxis().SetLabelSize(0.07)
            data.GetXaxis().SetTitleSize(0.09)
            data.GetXaxis().SetLabelOffset(0.05)
            if logyFlag == True:
                data.SetMinimum(1e-3)
            
            data.GetYaxis().SetNdivisions(508)
            
            return data

        def _make_sub_pads():
            if not dataOff:
                main_pad = ROOT.TPad(data.GetName()+'_main',data.GetName()+'_main',0, 0.35, 1, 1)
                sub_pad  = ROOT.TPad(data.GetName()+'_sub',data.GetName()+'_sub',0, 0, 1, 0.35)
            else:
                main_pad = ROOT.TPad(data.GetName()+'_main',data.GetName()+'_main',0, 0.1, 1, 1)
                sub_pad  = ROOT.TPad(data.GetName()+'_sub',data.GetName()+'_sub',0, 0, 0, 0)

            main_pad.SetBottomMargin(0.04)
            main_pad.SetLeftMargin(0.17)
            main_pad.SetRightMargin(0.05)
            main_pad.SetTopMargin(0.1)

            sub_pad.SetLeftMargin(0.17)
            sub_pad.SetRightMargin(0.05)
            sub_pad.SetTopMargin(0)
            sub_pad.SetBottomMargin(0.35)

            if logyFlag == True:
                main_pad.SetLogy()
            
            main_pad.Draw()
            sub_pad.Draw()

            return main_pad, sub_pad

        def _make_legend():
            legend = ROOT.TLegend(0.65,legend_topY,0.90,0.88)
            legend.SetBorderSize(0)
            if not dataOff: legend.AddEntry(data,data.GetName().split(',')[0],datastyle)
            return legend

        def _make_totalBkg():           
            totalBkg.SetMarkerStyle(0)
            totalBkg_err = totalBkg.Clone()
            totalBkg.SetLineColor(ROOT.kBlack)
            totalBkg_err.SetLineColor(ROOT.kBlack)
            totalBkg_err.SetLineWidth(0)
            totalBkg_err.SetFillColor(ROOT.kBlack)
            totalBkg_err.SetFillStyle(3354)
            if logyFlag: 
                totalBkg.SetMinimum(1e-3)
                totalBkg_err.SetMinimum(1e-3)
            legend.AddEntry(totalBkg_err,'Total bkg unc.','F')

            return totalBkg, totalBkg_err

        def _make_stack():
            stack = ROOT.THStack(data.GetName()+'_stack',data.GetName()+'_stack')
            # Build the stack
            legend_info = collections.OrderedDict()
            for bkg in bkgs:     # Won't loop if bkglist is empty
                if logyFlag:
                    bkg.SetMinimum(1e-3)

                stack.Add(bkg)
                legend_info[bkg.GetName().split(',')[0]] = bkg

            if logyFlag:
                stack.SetMinimum(1e-3) 

            # Deal with legend which needs ordering reversed from stack build
            legend_duplicates = []
            for bname in reversed(legend_info.keys()):
                if bname not in legend_duplicates:
                    legend.AddEntry(legend_info[bname],bname,'f')
                    legend_duplicates.append(bname)
            
            return stack

        def _make_signals():
            sigs_to_plot = signals
            # Can add together for total signal
            if addSignals:
                totsig = signals[0].Clone()
                for isig in range(1,len(signals)):
                    totsig.Add(signals[isig])
                sigs_to_plot = [totsig]

            # Plot either way
            for isig,sig in enumerate(sigs_to_plot):
                sig.SetLineWidth(2)
                if logyFlag == True:
                    sig.SetMinimum(1e-3)

                legend.AddEntry(sig,sig.GetTitle().split(',')[0],'L')
            
            return sigs_to_plot

        def _draw_subtitle_tex():
            subtitle_tex = ROOT.TLatex()
            subtitle_tex.SetNDC()
            subtitle_tex.SetTextAngle(0)
            subtitle_tex.SetTextColor(ROOT.kBlack)
            subtitle_tex.SetTextFont(42)
            subtitle_tex.SetTextAlign(12) 
            subtitle_tex.SetTextSize(0.06)
            subtitle_tex.DrawLatex(0.208,0.74,subtitle)

        def _draw_extralumi_tex():
            lumiE = ROOT.TLatex()
            lumiE.SetNDC()
            lumiE.SetTextAngle(0)
            lumiE.SetTextColor(ROOT.kBlack)
            lumiE.SetTextFont(42)
            lumiE.SetTextAlign(31) 
            lumiE.SetTextSize(0.7*0.1)
            lumiE.DrawLatex(1-0.05,1-0.1+0.2*0.1,"137 fb^{-1} (13 TeV)")

        def _draw_just_data():
            data.SetMaximum(1.13*data.GetMaximum())
            data.SetTitle(title)
            data.SetTitleOffset(1.1)
            data.Draw(datastyle)
            CMS_lumi.CMS_lumi(pad, year, 11)

        def _get_nsignals():
            if len(signals) == 0: nsignals = 0
            elif addSignals:      nsignals = 1
            else:                 nsignals = len(signals[0])
            return nsignals

        pad = self._make_pad_gen(outname)
        data = _format_data()

        if len(bkgs) == 0:
            _draw_just_data()
        else:
            main_pad, sub_pad = _make_sub_pads()
            nsignals = _get_nsignals()
            legend_topY = 0.73-0.03*(min(len(bkgs[0]),6)+nsignals+1)
            legend = _make_legend()
            totalBkg, totalBkg_err = _make_totalBkg()
            signals = signals if len(signals) == 0 else _make_signals()        
            stack = _make_stack()
            set_hist_maximums([stack, totalBkg, data], 2.5-legend_topY+0.03)

            # Go to main pad and draw
            main_pad.cd()
            data.Draw(datastyle)
            stack.Draw('same hist') # need to draw twice because the axis doesn't exist for modification until drawing
            try:    stack.GetYaxis().SetNdivisions(508)
            except: stack.GetYaxis().SetNdivisions(8,5,0)
            stack.Draw('same hist')
            for sig in signals:
                sig.Draw('hist same')

            # Draw total hist and error
            totalBkg.Draw('hist same')
            totalBkg_err.Draw('e2 same')
            legend.Draw()
            if not dataOff:
                data.Draw(datastyle+' same')

            ROOT.gPad.RedrawAxis()

            sub_pad.cd()
            pull = MakePullPlot(data,totalBkg)
            pull.Draw('hist')

            CMS_lumi.extraText = extraText
            CMS_lumi.cmsTextSize = 0.9
            CMS_lumi.cmsTextOffset = 2
            CMS_lumi.lumiTextSize = 0.9
            CMS_lumi.CMS_lumi(main_pad, year, 11)
            main_pad.cd()           
            _draw_subtitle_tex()
            
        self._save_pad_generic(outname, pad, saveROOT, savePDF, savePNG)
        return pad

    def MakeCan(self, outname, pads, saveROOT=False, savePDF=True, savePNG=True, padx=0, pady=0):
        '''Combine multiple pads/canvases into one canvas for convenience of viewing.
        Input pad order matters.

        Args:
            outname (str): Output file path name.
            pads ([TCanvas,TPad]): List of canvases/pads to plot together on one canvas.
            saveROOT (bool, optional): Save to master ROOT file. Defaults to False.
            savePDF (bool, optional): Save to PDF. Defaults to True.
            savePNG (bool, optional): Save to PNG. Defaults to True.

        Raises:
            RuntimeError: If 10 or more subdivisions are requested.

        Returns:
            ROOT.TCanvas: Output canvas.
        '''
        if padx == 0 and pady == 0:
            if len(pads) == 1:
                width = 800; height = 700; padx = 1; pady = 1
            elif len(pads) == 2:
                width = 1200; height = 700; padx = 2; pady = 1
            elif len(pads) == 3:
                width = 1800; height = 600; padx = 3; pady = 1
            elif len(pads) == 4:
                width = 1200; height = 1000; padx = 2; pady = 2
            elif len(pads) <= 6:
                width = 1600; height = 1000; padx = 3; pady = 2
            elif len(pads) <= 9:
                width = 1600; height = 1600; padx = 3; pady = 3
            else:
                raise RuntimeError('histlist of size %s not currently supported'%(len(pads),pads))

        canvas = ROOT.TCanvas(outname, outname, width, height)
        canvas.cd()
        canvas.Divide(padx, pady)
        for i in range(len(pads)):
            canvas.cd(i+1)
            pads[i].DrawClonePad()

        self._save_pad_generic(outname, canvas, saveROOT, savePDF, savePNG)
        return canvas

def _get_start_stop(i,slice_idxs):
    start = slice_idxs[i]+1
    stop  = slice_idxs[i+1]
    return start, stop

def plotAllFitResults(twoD,fittag):
    plotter = Plotter(twoD,fittag)
    plotter.plot_2D_distributions()
    plotter.plot_projections(logyFlag=False)
    plotter.plot_projections(logyFlag=True)
    plotter.plot_pre_vs_post()
    # plotter.plot_transfer_funcs()

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
    tree_fit = fd_file.Get('tree_fit_sb')
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

def MakePullPlot(data, bkg):
    bkgup, bkgdown = MakeUpDown(bkg)
    pull = data.Clone(data.GetName()+"_pull")
    pull.Add(bkg,-1)
    
    sigma = 0.0
    for ibin in range(1,pull.GetNbinsX()+1):
        d = data.GetBinContent(ibin)
        b = bkg.GetBinContent(ibin)
        if d >= b:
            derr = data.GetBinErrorLow(ibin)
            berr = abs(bkgup.GetBinContent(ibin)-bkg.GetBinContent(ibin))
        elif d < b:
            derr = data.GetBinErrorUp(ibin)
            berr = abs(bkgdown.GetBinContent(ibin)-bkg.GetBinContent(ibin))
        
        if d == 0:
            derr = 1

        sigma = math.sqrt(derr*derr + berr*berr)
        if sigma != 0:
            pull.SetBinContent(ibin, (pull.GetBinContent(ibin))/sigma)
        else:
            pull.SetBinContent(ibin, 0.0 )

    pull.SetFillColor(ROOT.kBlue)
    pull.SetTitle(";"+data.GetXaxis().GetTitle()+";(Data-Bkg)/#sigma")
    pull.SetStats(0)

    pull.GetYaxis().SetRangeUser(-2.9,2.9)
    pull.GetYaxis().SetTitleOffset(0.4)                             
    pull.GetYaxis().SetLabelSize(0.13)
    pull.GetYaxis().SetTitleSize(0.12)
    pull.GetYaxis().SetNdivisions(306)

    pull.GetXaxis().SetLabelSize(0.13)
    pull.GetXaxis().SetTitleSize(0.15)
    pull.GetXaxis().SetTitle(data.GetXaxis().GetTitle())
    pull.GetYaxis().SetTitle("(Data-Bkg)/#sigma")
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
    diffnuis_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnosticsTest.root --abs -g nuisance_pulls.root'
    execute_cmd(diffnuis_cmd)
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
    return successful_fits

def SavePostFitParametricFuncVals():
    fit_result_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
    goodFitTags = _getGoodFitResults(fit_result_file)
    for fittag in goodFitTags:
        coeffs_final = fit_result_file.Get('fit_'+fittag).floatParsFinal()
        all_par_names = [] # get names of anything matching *_par<N>
        for i in range(coeffs_final.getSize()):
            name = coeffs_final.at(i).GetName()
            if name.split('_')[-1].startswith('par'):
                all_par_names.append(name)
        all_par_names.sort()
        # Get unique prefixes to _par<N>
        all_obj_names = set(['_'.join(name.split('_')[:-1]) for name in all_par_names])

        for obj_name in all_obj_names:
            with open('rpf_params_%s_fit%s.txt'%(obj_name,fittag),'w') as param_out:
                for par_name in [p for p in all_par_names if p.startswith(obj_name)]:
                    var = coeffs_final.find(par_name)
                    param_out.write('%s: %s +/- %s\n'%(par_name, var.getValV(), var.getError()))

def GenPostFitShapes():
    fit_result_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
    goodFitTags = _getGoodFitResults(fit_result_file)
    for t in goodFitTags:
        shapes_cmd = 'PostFit2DShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -o postfitshapes_{0}.root -f fitDiagnosticsTest.root:fit_{0} --postfit --sampling --samples 100 --print 2> PostFitShapes2D_stderr_{0}.txt'.format(t)
        execute_cmd(shapes_cmd)

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