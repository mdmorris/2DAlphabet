import itertools, ROOT
from TwoDAlphabet.src.helpers import nested_dict, set_hist_maximums

def MakeSystematicPlots(configObj):
    '''Make plots of the systematic shape variations of each process based on those
    processes and systematic shapes specified in the config. Shapes are presented 
    as projections onto 1D axis where no selection has been made on the axis not
    being plotted. Plots are saved to UncertPlots/.
    '''
    regions = ['pass','fail']
    variations = ['nom','up','down']
    axes = ['X','Y']
    for proc in configObj.Section('PROCESS'):
        for syst in configObj.processes[proc]['SYSTEMATICS']:
            if configObj.systematics[syst]['CODE'] < 2: continue

            tracking_dict = nested_dict(3,None)
            for r,v,x in itertools.product(regions,variations,axes):
                if v == 'nom': h = configObj.orgFile.Get(configObj.organizedDict[proc][r+'_FULL']['nominal'])
                else:          h = configObj.orgFile.Get(configObj.organizedDict[proc][r+'_FULL'][syst+v.capitalize()])
                if x == 'Y':   tracking_dict[r][x][v] = h.ProjectionY(proc +'_'+r+ '_'+syst+'_'+x+'_'+v)
                elif x == 'X': tracking_dict[r][x][v] = h.ProjectionX(proc +'_'+r+ '_'+syst+'_'+x+'_'+v)

            for r,x in itertools.product(regions,axes):
                thisCan = ROOT.TCanvas('canvas_'+proc+'_'+syst,'canvas_'+proc+'_'+syst,800,700)

                nom, up, down = (tracking_dict[r][x]['nom'],
                                    tracking_dict[r][x]['up'],
                                    tracking_dict[r][x]['down'])

                nom.SetLineColor(ROOT.kBlack)
                nom.SetFillColor(ROOT.kYellow-9)
                up.SetLineColor(ROOT.kRed)
                down.SetLineColor(ROOT.kBlue)

                up.SetLineStyle(9)
                down.SetLineStyle(9)
                up.SetLineWidth(2)
                down.SetLineWidth(2)

                nom,up,down = set_hist_maximums([nom,up,down])
                
                if x == 'X': nom.SetXTitle(configObj.inputConfig['BINNING']['X']['TITLE'])
                elif x == 'Y': nom.SetXTitle(configObj.inputConfig['BINNING']['Y']['TITLE'])

                nom.SetTitle('')
                nom.GetXaxis().SetTitleOffset(1.0)
                nom.GetXaxis().SetTitleSize(0.05)
                thisCan.SetRightMargin(0.16)

                nom.Draw('hist'); up.Draw('same hist'); down.Draw('same hist')
                thisCan.Print(configObj.projPath+'/UncertPlots/Uncertainty_'+proc+'_'+syst+r+x+'.png','png')

def plotFitResults(self,fittag):#,simfit=False): # fittag means 'b' or 's'
    allVars = []

    #####################
    #   Get everything  #
    #####################

    # File with histograms and RooFitResult parameters
    # if simfit == False:
    #     post_file = TFile.Open(self.projPath+'/postfitshapes_'+fittag+'.root')
    #     fd_file = TFile.Open(self.projPath+'/fitDiagnosticsTest.root')
    # else:
    post_file = TFile.Open(self.tag+'/postfitshapes_'+fittag+'.root')
    axis_hist = post_file.Get('pass_LOW_'+self.name+'_prefit/data_obs')
    fd_file = TFile.Open(self.tag+'/fitDiagnosticsTest.root')

    if 'RunII_' in fittag:
        runII = True
        fittag = fittag.replace('RunII_','')
    else:
        runII = False

    fit_result = fd_file.Get('fit_'+fittag)

    x_low = self.newXbins['LOW'][0]
    x_high = self.newXbins['HIGH'][-1]

    print ('Finding start and end bin indexes of signal range. Looking for '+str(self.sigStart)+', '+str(self.sigEnd))
    for ix,xwall in enumerate(self.fullXbins):
        if xwall == self.sigStart:
            print ('Assigning start bin as '+str(ix+1))
            x_sigstart_bin = ix+1
        if xwall == self.sigEnd:
            print ('Assigning end bin as '+str(ix))
            x_sigend_bin = ix

    y_low = self.newYbins[0]
    y_high = self.newYbins[-1]
    y_nbins = len(self.newYbins)-1

    if self.ySlices == False:
        # Define low, middle, high projection regions for y (x regions defined already via signal region bounds)
        y_turnon_endBin = axis_hist.ProjectionY().GetMaximumBin()
        y_turnon_endVal = int(axis_hist.GetYaxis().GetBinUpEdge(y_turnon_endBin))
        y_tail_beginningBin = axis_hist.GetYaxis().FindBin((y_high - y_turnon_endVal)/3.0 + y_turnon_endVal)
        
        if y_turnon_endBin > y_nbins/2.0:  # in case this isn't a distribution with a turn-on
            y_turnon_endBin = int(round(y_nbins/3.0))
            y_turnon_endVal = int(axis_hist.GetYaxis().GetBinUpEdge(y_turnon_endBin))
            y_tail_beginningBin = 2*y_turnon_endBin

        y_tail_beginningVal = str(int(axis_hist.GetYaxis().GetBinLowEdge(y_tail_beginningBin)))
    
    else: # custom slices
        y_low = self.ySlices[0]
        y_turnon_endVal = self.ySlices[1]
        y_turnon_endBin = axis_hist.GetYaxis().FindBin(y_turnon_endVal)-1
        y_tail_beginningVal = self.ySlices[2]
        y_tail_beginningBin = axis_hist.GetYaxis().FindBin(y_tail_beginningVal)
        y_high = self.ySlices[3]
    
    y_turnon_endVal = str(y_turnon_endVal)
    y_tail_beginningVal = str(y_tail_beginningVal)
    
    # Save out slice edges in case they aren't plotted in title (don't do for x since those are defined by user in config)
    y_slices_out = open(self.projPath+'/y_slices.txt','w')
    y_slices_out.write('Bin edge values:  %i %s %s %i\n'%(y_low, y_turnon_endVal, y_tail_beginningVal, y_high))
    y_slices_out.write('Bin edge numbers: %i %i %i %i'%(1,     y_turnon_endBin+1, y_tail_beginningBin, len(self.newYbins)-1))
    y_slices_out.close()

    # Final fit signal strength
    if fittag == 's':
        tree_fit_sb = fd_file.Get('tree_fit_sb')
        tree_fit_sb.GetEntry(0)
        signal_strength = tree_fit_sb.r
    else:
        tree_fit_b = fd_file.Get('tree_fit_b')
        tree_fit_b.GetEntry(0)
        signal_strength = tree_fit_b.r

    #####################
    #    Data vs Bkg    #
    #####################

    hist_dict = {}

    # Organize and make any projections or 2D distributions
    for process in [process for process in self.inputConfig['PROCESS'] if process != 'HELP']+['qcd','TotalBkg']:
        hist_dict[process] = {}
        for cat in ['fail','pass']:
            hist_dict[process][cat] = {'LOW':{},'SIG':{},'HIGH':{}}
            x_slice_list_pre = []
            x_slice_list_post = []
            # Grab everything and put clones in a dictionary
            for c in ['LOW','SIG','HIGH']:
                file_dir = cat+'_'+c+'_'+self.name
                hist_dict[process][cat]['prefit_'+c] = post_file.Get(file_dir+'_prefit/'+process).Clone()
                hist_dict[process][cat]['postfit_'+c] = post_file.Get(file_dir+'_postfit/'+process).Clone()
                x_slice_list_pre.append(hist_dict[process][cat]['prefit_'+c])    # lists for 2D making
                x_slice_list_post.append(hist_dict[process][cat]['postfit_'+c])

            # First rebuild the 2D distributions
            if self.blindedPlots and process == 'data_obs':
                hist_dict[process][cat]['prefit_2D'] = header.stitchHistsInX(process+'_'+cat+'_prefit2D',self.fullXbins,self.newYbins,x_slice_list_pre,blinded=[1])
                hist_dict[process][cat]['postfit_2D'] = header.stitchHistsInX(process+'_'+cat+'_postfit2D',self.fullXbins,self.newYbins,x_slice_list_post,blinded=[1])

            else:
                hist_dict[process][cat]['prefit_2D'] = header.stitchHistsInX(process+'_'+cat+'_prefit2D',self.fullXbins,self.newYbins,x_slice_list_pre,blinded=[])
                hist_dict[process][cat]['postfit_2D'] = header.stitchHistsInX(process+'_'+cat+'_postfit2D',self.fullXbins,self.newYbins,x_slice_list_post,blinded=[])

            hist_dict[process][cat]['prefit_2D'].SetMinimum(0)
            hist_dict[process][cat]['postfit_2D'].SetMinimum(0)
            hist_dict[process][cat]['prefit_2D'].SetTitle(process + ', ' + cat +', '+self.name+ ', pre-fit')
            hist_dict[process][cat]['postfit_2D'].SetTitle(process + ', ' + cat +', '+self.name+ ', post-fit')

            # Now projections
            base_proj_name_pre = process+'_'+cat+'_'+self.name+'_pre_'
            base_proj_name_post = process+'_'+cat+'_'+self.name+'_post_'

            hist_dict[process][cat]['prefit_projx1'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+str(y_low)+'-'+y_turnon_endVal,              1,                      y_turnon_endBin, 'e')
            hist_dict[process][cat]['prefit_projx2'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+y_turnon_endVal+'-'+y_tail_beginningVal,     y_turnon_endBin+1,      y_tail_beginningBin-1,'e')
            hist_dict[process][cat]['prefit_projx3'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+y_tail_beginningVal+'-'+str(y_high),         y_tail_beginningBin,  y_nbins,'e')

            hist_dict[process][cat]['prefit_projy1'] = hist_dict[process][cat]['prefit_LOW'].ProjectionY(base_proj_name_pre+'projy_'+str(x_low)+'-'+str(self.sigStart),          1,                      hist_dict[process][cat]['prefit_LOW'].GetNbinsX(),'e')
            if self.blindedPlots:
                hist_dict[process][cat]['prefit_projy2'] = hist_dict[process][cat]['prefit_2D'].ProjectionY(base_proj_name_pre+'projy_'+str(self.sigStart)+'-'+str(self.sigEnd), x_sigstart_bin,           x_sigend_bin,'e')
            else:
                hist_dict[process][cat]['prefit_projy2'] = hist_dict[process][cat]['prefit_SIG'].ProjectionY(base_proj_name_pre+'projy_'+str(self.sigStart)+'-'+str(self.sigEnd),    1,                      hist_dict[process][cat]['prefit_SIG'].GetNbinsX(),'e')
            hist_dict[process][cat]['prefit_projy3'] = hist_dict[process][cat]['prefit_HIGH'].ProjectionY(base_proj_name_pre+'projy_'+str(self.sigEnd)+'-'+str(x_high),          1,                      hist_dict[process][cat]['prefit_HIGH'].GetNbinsX(),'e')

            hist_dict[process][cat]['postfit_projx1'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+str(y_low)+'-'+y_turnon_endVal,           1,                      y_turnon_endBin, 'e')
            hist_dict[process][cat]['postfit_projx2'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+y_turnon_endVal+'-'+y_tail_beginningVal,  y_turnon_endBin+1,      y_tail_beginningBin-1,'e')
            hist_dict[process][cat]['postfit_projx3'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+y_tail_beginningVal+'-'+str(y_high),      y_tail_beginningBin,  y_nbins,'e')

            hist_dict[process][cat]['postfit_projy1'] = hist_dict[process][cat]['postfit_LOW'].ProjectionY(base_proj_name_post+'projy_'+str(x_low)+'-'+str(self.sigStart),       1,                      hist_dict[process][cat]['postfit_LOW'].GetNbinsX(),'e')
            if self.blindedPlots:
                hist_dict[process][cat]['postfit_projy2'] = hist_dict[process][cat]['postfit_2D'].ProjectionY(base_proj_name_pre+'projy_'+str(self.sigStart)+'-'+str(self.sigEnd), x_sigstart_bin,           x_sigend_bin,'e')
            else:
                hist_dict[process][cat]['postfit_projy2'] = hist_dict[process][cat]['postfit_SIG'].ProjectionY(base_proj_name_post+'projy_'+str(self.sigStart)+'-'+str(self.sigEnd), 1,                      hist_dict[process][cat]['postfit_SIG'].GetNbinsX(),'e')
            hist_dict[process][cat]['postfit_projy3'] = hist_dict[process][cat]['postfit_HIGH'].ProjectionY(base_proj_name_post+'projy_'+str(self.sigEnd)+'-'+str(x_high),       1,                      hist_dict[process][cat]['postfit_HIGH'].GetNbinsX(),'e')

            x_edges = [x_low,self.sigStart,self.sigEnd,x_high]
            y_edges = [y_low,y_turnon_endVal,y_tail_beginningVal,y_high]

            for z in ['x','y']:
                for i in range(1,4):
                    hist_dict[process][cat]['postfit_proj'+z+str(i)].SetMinimum(0)
                    if z == 'x':
                        hist_dict[process][cat]['postfit_proj'+z+str(i)].SetTitle(process + ', ' + cat+', '+self.name+ ', ' +str(y_edges[i-1]) +'-'+ str(y_edges[i]))
                    elif z == 'y':
                        hist_dict[process][cat]['postfit_proj'+z+str(i)].SetTitle(process + ', ' + cat+', '+self.name+ ', ' +str(x_edges[i-1]) +'-'+ str(x_edges[i]))

    post_file.Close()
        
    process_list = hist_dict.keys()

    # Create lists for the 2D projections (ones you want to see together)
    for process in hist_dict.keys():    # Canvas
        isSignal = (process != 'qcd' and process != 'TotalBkg' and self.inputConfig['PROCESS'][process]['CODE'] == 0)
        twoDList = []
        twoDtitles = []   
        for cat in ['fail','pass']:
            for fit in ['prefit', 'postfit']:
                if isSignal and fittag == 's' and fit == 'postfit':
                    hist_dict[process][cat][fit+'_2D'].Scale(signal_strength)

                twoDList.append(hist_dict[process][cat][fit+'_2D'])
                twoDtitles.append(process + ', '+cat+', '+fit)

        if isSignal and fittag != 's':
            continue
        else:
            header.makeCan('plots/fit_'+fittag+'/'+process+'_fit'+fittag+'_2D',self.projPath,twoDList,titles=twoDtitles,xtitle=self.xVarTitle,ytitle=self.yVarTitle,year=self.year)

    # Invert the last two items (unique to b*) - customize as needed
    process_list[-1],process_list[-2] = process_list[-2],process_list[-1]

    # Get the colors
    colors = []
    for process in process_list:
        if process != 'data_obs':
            if process not in self.inputConfig['PROCESS'].keys():
                continue
            if (process != 'qcd' and process !='TotalBkg' and self.inputConfig['PROCESS'][process]['CODE'] != 0):
                if (process in self.inputConfig['PROCESS'].keys()) and ('COLOR' in self.inputConfig['PROCESS'][process].keys()):
                    colors.append(self.inputConfig['PROCESS'][process]["COLOR"])
                else:
                    colors.append(None)

    # Put QCD on bottom of stack since it's smooth
    colors_logy = colors+[kYellow]
    colors = [kYellow]+colors

    # Create lists for makeCan of the projections
    for plotType in ['postfit_projx','postfit_projy']:   # Canvases
        bkgList = []
        bkgNameList = []
        bkgList_logy = []
        bkgNameList_logy = []
        dataList = []
        signalList = []
        titleList = []
        totalBkgs = []
        for cat in ['fail','pass']: # Row 
            for regionNum in range(1,4):    # Column
                bkg_process_list = []
                bkg_process_names = []
                signal_process_list = []
                signal_names = []
                sliceranges = []
                this_totalbkg = hist_dict['TotalBkg'][cat][plotType+str(regionNum)]
                totalBkgs.append(this_totalbkg)
                if 'y' in plotType:
                    if regionNum == 1: low_str,high_str = str(x_low),str(self.sigStart)
                    elif regionNum == 2: low_str,high_str = str(self.sigStart),str(self.sigEnd)
                    elif regionNum == 3: low_str,high_str = str(self.sigEnd),str(x_high)
                elif 'x' in plotType:
                    if regionNum == 1: low_str,high_str = str(y_low),str(y_turnon_endVal)
                    elif regionNum == 2: low_str,high_str = str(y_turnon_endVal),str(y_tail_beginningVal)
                    elif regionNum == 3: low_str,high_str = str(y_tail_beginningVal),str(y_high)
                for process in process_list:
                    if process != 'data_obs':
                        if (process != 'qcd' and process != 'TotalBkg' and self.inputConfig['PROCESS'][process]['CODE'] != 0):
                            process_name = process if 'TITLE' not in self.inputConfig['PROCESS'][process].keys() else self.inputConfig['PROCESS'][process]['TITLE']
                            bkg_process_list.append(hist_dict[process][cat][plotType+str(regionNum)])
                            bkg_process_names.append(process_name)
                        elif (process != 'qcd' and process != 'TotalBkg' and self.inputConfig['PROCESS'][process]['CODE'] == 0):
                            process_name = process if 'TITLE' not in self.inputConfig['PROCESS'][process].keys() else self.inputConfig['PROCESS'][process]['TITLE']
                            if self.plotPrefitSigInFitB and fittag == 'b':
                                signal_process_list.append(hist_dict[process][cat][plotType.replace("postfit","prefit")+str(regionNum)])
                            else:
                                hist_dict[process][cat][plotType+str(regionNum)].Scale(signal_strength)
                                signal_process_list.append(hist_dict[process][cat][plotType+str(regionNum)])
                            signal_names.append(process_name)
                            
                    else:
                        this_data = hist_dict[process][cat][plotType+str(regionNum)]
                        dataList.append(this_data)

                # Put QCD on bottom of stack since it's smooth
                bkg_process_list = [hist_dict['qcd'][cat][plotType+str(regionNum)]]+bkg_process_list
                bkgNameList.append(['Multijet']+bkg_process_names)
                bkgList.append(bkg_process_list)

                # Put QCD on top of logy
                bkg_process_list_logy = bkg_process_list[1:]+[hist_dict['qcd'][cat][plotType+str(regionNum)]]
                bkgNameList_logy.append(bkg_process_names+['Multijet'])
                bkgList_logy.append(bkg_process_list_logy)

                # Do same for signal
                signalList.append(signal_process_list)

                if self.plotTitles: titleList.append('Data vs bkg - %s - [%s,%s]'%(cat,low_str,high_str))
                else: titleList.append('')

                sliceranges = [h.GetName().split('_')[-1] for h in dataList]

        if self.plotEvtsPerUnit:
            new_dataList = []
            new_bkgList = []
            new_bkgList_logy = []
            new_signalList = []
            new_totalBkgs = [] 
            for i,h in enumerate(dataList):
                new_dataList.append(ConvertToEvtsPerUnit(h))
                new_totalBkgs.append(ConvertToEvtsPerUnit(totalBkgs[i]))
                new_sub_bkgList = []
                new_sub_bkgList_logy = []
                new_sub_signalList = []
                for b in bkgList[i]:
                    new_sub_bkgList.append(ConvertToEvtsPerUnit(b))
                for blog in bkgList_logy[i]:
                    new_sub_bkgList_logy.append(ConvertToEvtsPerUnit(blog))
                for s in signalList[i]:
                    new_sub_signalList.append(ConvertToEvtsPerUnit(s))

                new_bkgList.append(new_sub_bkgList)
                new_bkgList_logy.append(new_sub_bkgList_logy)
                new_signalList.append(new_sub_signalList)
            dataList = new_dataList
            bkgList = new_bkgList
            bkgList_logy = new_bkgList_logy
            signalList = new_signalList
            totalBkgs = new_totalBkgs

            # bstar specific (GeV specifically)
            yAxisTitle = 'Events / %s GeV' % header.GetMinWidth(dataList[0])

        else:
            yAxisTitle = 'Events / bin'

        root_out = TFile.Open(self.projPath+'/plots/fit_'+fittag+'/'+plotType+'_fit'+fittag+'.root','RECREATE')
        for h in dataList+bkgList+totalBkgs+signalList:
            if isinstance(h,list):
                for i in h:
                    root_out.WriteObject(i,i.GetName())
            else:
                root_out.WriteObject(h,h.GetName())
        root_out.Close()

        if 'x' in plotType:
            header.makeCan('plots/fit_'+fittag+'/'+plotType+'_fit'+fittag,self.projPath,
                dataList,bkglist=bkgList,subtitles=sliceranges,sliceVar=self.yVarTitle,totalBkg=totalBkgs,signals=signalList,
                bkgNames=bkgNameList,signalNames=signal_names,titles=titleList,
                colors=colors,xtitle=self.xVarTitle,ytitle=yAxisTitle,year=self.year,addSignals=self.addSignals)
            header.makeCan('plots/fit_'+fittag+'/'+plotType+'_fit'+fittag+'_log',self.projPath,
                dataList,bkglist=bkgList_logy,subtitles=sliceranges,sliceVar=self.yVarTitle,totalBkg=totalBkgs,signals=signalList,
                bkgNames=bkgNameList_logy,signalNames=signal_names,titles=titleList,
                colors=colors_logy,xtitle=self.xVarTitle,ytitle=yAxisTitle,logy=True,year=self.year,addSignals=self.addSignals)
        elif 'y' in plotType:
            header.makeCan('plots/fit_'+fittag+'/'+plotType+'_fit'+fittag,self.projPath,
                dataList,bkglist=bkgList,subtitles=sliceranges,sliceVar=self.xVarTitle,totalBkg=totalBkgs,signals=signalList,
                bkgNames=bkgNameList,signalNames=signal_names,titles=titleList,
                colors=colors,xtitle=self.yVarTitle,ytitle=yAxisTitle,year=self.year,addSignals=self.addSignals)
            header.makeCan('plots/fit_'+fittag+'/'+plotType+'_fit'+fittag+'_log',self.projPath,
                dataList,bkglist=bkgList_logy,subtitles=sliceranges,sliceVar=self.xVarTitle,totalBkg=totalBkgs,signals=signalList,
                bkgNames=bkgNameList_logy,signalNames=signal_names,titles=titleList,
                colors=colors_logy,xtitle=self.yVarTitle,ytitle=yAxisTitle,logy=True,year=self.year,addSignals=self.addSignals)

    # Make comparisons for each background process of pre and post fit projections
    for plotType in ['projx','projy']:
        for process in process_list:
            if process != 'data_obs' and process != 'TotalBkg':
                pre_list = []
                post_list = []
                titleList = []
                for cat in ['fail','pass']: # Row 
                    for regionNum in range(1,4):    # Column
                        if 'y' in plotType:
                            if regionNum == 1: low_str,high_str = str(x_low),str(self.sigStart)
                            elif regionNum == 2: low_str,high_str = str(self.sigStart),str(self.sigEnd)
                            elif regionNum == 3: low_str,high_str = str(self.sigEnd),str(x_high)
                        elif 'x' in plotType:
                            if regionNum == 1: low_str,high_str = str(y_low),str(y_turnon_endVal)
                            elif regionNum == 2: low_str,high_str = str(y_turnon_endVal),str(y_tail_beginningVal)
                            elif regionNum == 3: low_str,high_str = str(y_tail_beginningVal),str(y_high)

                        pre_list.append([hist_dict[process][cat]['prefit_'+plotType+str(regionNum)]])  # in terms of makeCan these are "bkg hists"
                        post_list.append(hist_dict[process][cat]['postfit_'+plotType+str(regionNum)])   # and these are "data hists"
                        if process != 'qcd' and process != 'TotalBkg':
                            if 'COLOR' in self.inputConfig['PROCESS'][process].keys():
                                prepostcolors = [self.inputConfig['PROCESS'][process]['COLOR']]
                            else:
                                prepostcolors = [0]
                        else:
                            prepostcolors = [kYellow]

                        titleList.append('Pre vs Postfit - %s - %s - [%s,%s]'%(process,cat,low_str,high_str))

                if 'x' in plotType: header.makeCan('plots/fit_'+fittag+'/'+process+'_'+plotType+'_fit'+fittag,self.projPath,
                    post_list,bkglist=pre_list,totalBkg=[b[0] for b in pre_list],
                    titles=titleList,bkgNames=['Prefit, '+process],dataName='        Postfit, '+process,
                    colors=prepostcolors,xtitle=self.xVarTitle,datastyle='histpe',year=self.year)
                if 'y' in plotType: header.makeCan('plots/fit_'+fittag+'/'+process+'_'+plotType+'_fit'+fittag,self.projPath,
                    post_list,bkglist=pre_list,totalBkg=[b[0] for b in pre_list],
                    titles=titleList,bkgNames=['Prefit, '+process],dataName='        Postfit, '+process,
                    colors=prepostcolors,xtitle=self.yVarTitle,datastyle='histpe',year=self.year)


    ##############
    #    Rp/f    #
    ##############      
    # Don't run Rp/f if this is just summed plots 
    if runII: return 0

    # Need to sample the space to get the Rp/f with proper errors (1000 samples)
    rpf_xnbins = len(self.fullXbins)-1
    rpf_ynbins = len(self.newYbins)-1
    if self.rpfRatio == False: rpf_zbins = [i/1000000. for i in range(0,1000001)]
    else: rpf_zbins = [i/1000. for i in range(0,5001)]
    rpf_samples = TH3F('rpf_samples','rpf_samples',rpf_xnbins, array.array('d',self.fullXbins), rpf_ynbins, array.array('d',self.newYbins), len(rpf_zbins)-1, array.array('d',rpf_zbins))# TH3 to store samples
    sample_size = 500

    # Collect all final parameter values
    param_final = fit_result.floatParsFinal()
    coeffs_final = RooArgSet()
    for v in self.rpf.funcVars.keys():
        coeffs_final.add(param_final.find(v))

    # Now sample to generate the Rpf distribution
    for i in range(sample_size):
        sys.stdout.write('\rSampling '+str(100*float(i)/float(sample_size)) + '%')
        sys.stdout.flush()
        param_sample = fit_result.randomizePars()

        # Set params of the Rpf object
        coeffIter_sample = param_sample.createIterator()
        coeff_sample = coeffIter_sample.Next()
        while coeff_sample:
            # Set the rpf parameter to the sample value
            if coeff_sample.GetName() in self.rpf.funcVars.keys():
                self.rpf.setFuncParam(coeff_sample.GetName(), coeff_sample.getValV())
            coeff_sample = coeffIter_sample.Next()

        # Loop over bins and fill
        for xbin in range(1,rpf_xnbins+1):
            for ybin in range(1,rpf_ynbins+1):
                bin_val = 0

                thisXCenter = rpf_samples.GetXaxis().GetBinCenter(xbin)
                thisYCenter = rpf_samples.GetYaxis().GetBinCenter(ybin)

                if self.recycleAll:
                    # Remap to [-1,1]
                    x_center_mapped = (thisXCenter - self.newXbins['LOW'][0])/(self.newXbins['HIGH'][-1] - self.newXbins['LOW'][0])
                    y_center_mapped = (thisYCenter - self.newYbins[0])/(self.newYbins[-1] - self.newYbins[0])

                    # And assign it to a RooConstVar 
                    x_const = RooConstVar("ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,x_center_mapped)
                    y_const = RooConstVar("ConstVar_y_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,y_center_mapped)
                    
                    # Now get the Rpf function value for this bin 
                    self.allVars.append(x_const)
                    self.allVars.append(y_const)
                    self.rpf.evalRpf(x_const, y_const,xbin,ybin)

                # Determine the category
                if thisXCenter > self.newXbins['LOW'][0] and thisXCenter < self.newXbins['LOW'][-1]: # in the LOW category
                    thisxcat = 'LOW'
                elif thisXCenter > self.newXbins['SIG'][0] and thisXCenter < self.newXbins['SIG'][-1]: # in the SIG category
                    thisxcat = 'SIG'
                elif thisXCenter > self.newXbins['HIGH'][0] and thisXCenter < self.newXbins['HIGH'][-1]: # in the HIGH category
                    thisxcat = 'HIGH'

                bin_val = self.rpf.getFuncBinVal(thisxcat,xbin,ybin)

                rpf_samples.Fill(thisXCenter,thisYCenter,bin_val)

    print ('\n')
    rpf_final = TH2F('rpf_final','rpf_final',rpf_xnbins, array.array('d',self.fullXbins), rpf_ynbins, array.array('d',self.newYbins))
    # Now loop over all x,y bin in rpf_samples, project onto Z axis, 
    # get the mean and RMS and set as the bin content and error in rpf_final
    for xbin in range(1,rpf_final.GetNbinsX()+1):
        for ybin in range(1,rpf_final.GetNbinsY()+1):
            temp_projz = rpf_samples.ProjectionZ('temp_projz',xbin,xbin,ybin,ybin)
            rpf_final.SetBinContent(xbin,ybin,temp_projz.GetMean())
            rpf_final.SetBinError(xbin,ybin,temp_projz.GetRMS())

    rpf_final.SetTitle('')
    rpf_final.GetXaxis().SetTitle(self.xVarTitle)
    rpf_final.GetYaxis().SetTitle(self.yVarTitle)
    rpf_final.GetZaxis().SetTitle('R_{P/F}' if self.rpfRatio == False else 'R_{Ratio}')
    rpf_final.GetXaxis().SetTitleSize(0.045)
    rpf_final.GetYaxis().SetTitleSize(0.045)
    rpf_final.GetZaxis().SetTitleSize(0.045)
    rpf_final.GetXaxis().SetTitleOffset(1.2)
    rpf_final.GetYaxis().SetTitleOffset(1.5)
    rpf_final.GetZaxis().SetTitleOffset(1.3)

    rpf_c = TCanvas('rpf_c','Post-fit R_{P/F}',1000,700)
    CMS_lumi.lumiTextSize = 0.75
    CMS_lumi.cmsTextSize = 0.85
    CMS_lumi.extraText = 'Preliminary'
    CMS_lumi.CMS_lumi(rpf_c, self.year, 11)
    rpf_c.SetRightMargin(0.2)
    rpf_final.Draw('colz')
    rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_colz.pdf','pdf')
    rpf_final.Draw('surf')
    rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_surf.pdf','pdf')
    rpf_final.Draw('pe')
    rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_errs.pdf','pdf')

    rpf_file = TFile.Open(self.projPath+'/plots/postfit_rpf_fit'+fittag+'.root','RECREATE')
    rpf_file.cd()
    rpf_final.Write()
    rpf_file.Close()


def makeCan(name, tag, histlist, bkglist=[],totalBkg=None,signals=[],colors=[],
            titles=[],subtitles=[],sliceVar='X',dataName='Data',bkgNames=[],signalNames=[],logy=False,
            rootfile=False,xtitle='',ytitle='',ztitle='',dataOff=False,
            datastyle='pe',year=1, addSignals=True, extraText=''):  
    # histlist is just the generic list but if bkglist is specified (non-empty)
    # then this function will stack the backgrounds and compare against histlist as if 
    # it is data. The imporant bit is that bkglist is a list of lists. The first index
    # of bkglist corresponds to the index in histlist (the corresponding data). 
    # For example you could have:
    #   histlist = [data1, data2]
    #   bkglist = [[bkg1_1,bkg2_1],[bkg1_2,bkg2_2]]

    if len(histlist) == 1:
        width = 800
        height = 700
        padx = 1
        pady = 1
    elif len(histlist) == 2:
        width = 1200
        height = 700
        padx = 2
        pady = 1
    elif len(histlist) == 3:
        width = 1800
        height = 600
        padx = 3
        pady = 1
    elif len(histlist) == 4:
        width = 1200
        height = 1000
        padx = 2
        pady = 2
    elif len(histlist) == 6 or len(histlist) == 5:
        height = 1600
        width = 1200
        padx = 2
        pady = 3
        histlist = reorderHists(histlist)
        if bkglist != []: bkglist = reorderHists(bkglist)
        if signals != []: signals = reorderHists(signals)
        if totalBkg != None: totalBkg = reorderHists(totalBkg)
        if titles != []: titles = reorderHists(titles)
        if subtitles != []: subtitles = reorderHists(subtitles)
    else:
        print ('histlist of size ' + str(len(histlist)) + ' not currently supported')
        print (histlist)
        return 0

    tdrstyle.setTDRStyle()
    gStyle.SetLegendFont(42)
    gStyle.SetTitleBorderSize(0)
    gStyle.SetTitleAlign(33)
    gStyle.SetTitleX(.77)
        
    myCan = TCanvas(name,name,width,height)
    myCan.Divide(padx,pady)

    # Just some colors that I think work well together and a bunch of empty lists for storage if needed
    default_colors = [kRed,kMagenta,kGreen,kCyan,kBlue]
    if len(colors) == 0:   
        colors = default_colors
    color_idx_order = None
    stacks = []
    tot_hists_err = []
    tot_hists = []
    legends = []
    mains = []
    subs = []
    pulls = []
    logString = ''
    tot_sigs = []

    # For each hist/data distribution
    for hist_index, hist in enumerate(histlist):
        # Grab the pad we want to draw in
        myCan.cd(hist_index+1)
        # if len(histlist) > 1:
        thisPad = myCan.GetPrimitive(name+'_'+str(hist_index+1))
        thisPad.cd()        
        thisPad.SetRightMargin(0.0)
        thisPad.SetTopMargin(0.0)
        thisPad.SetBottomMargin(0.0)

        # If this is a TH2, just draw the lego
        if hist.ClassName().find('TH2') != -1:
            gPad.SetLeftMargin(0.15)
            gPad.SetRightMargin(0.2)
            gPad.SetBottomMargin(0.12)
            gPad.SetTopMargin(0.1)
            if logy: gPad.SetLogz()
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)
            hist.GetZaxis().SetTitle(ztitle)
            hist.GetXaxis().SetTitleOffset(1.15)
            hist.GetYaxis().SetTitleOffset(1.5)
            hist.GetZaxis().SetTitleOffset(1.5)
            hist.GetYaxis().SetLabelSize(0.05)
            hist.GetYaxis().SetTitleSize(0.05)
            hist.GetXaxis().SetLabelSize(0.05)
            hist.GetXaxis().SetTitleSize(0.05)
            hist.GetZaxis().SetLabelSize(0.05)
            hist.GetZaxis().SetTitleSize(0.05)
            hist.GetXaxis().SetNdivisions(505)
            # hist.GetXaxis().SetLabelOffset(0.02)
            if 'lego' in datastyle.lower(): hist.GetZaxis().SetTitleOffset(1.4)
            if len(titles) > 0:
                hist.SetTitle(titles[hist_index])

            if datastyle != 'pe': hist.Draw(datastyle)
            else: hist.Draw('colz')
            if len(bkglist) > 0:
                raise TypeError('ERROR: It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.')

            CMS_lumi.extraText = extraText
            CMS_lumi.CMS_lumi(thisPad, year, 11, sim=False if 'data' in name.lower() else True)
        
        # Otherwise it's a TH1 hopefully
        else:
            titleSize = 0.09
            alpha = 1
            if dataOff:
                alpha = 0
            hist.SetLineColorAlpha(kBlack,alpha)
            if 'pe' in datastyle.lower():
                hist.SetMarkerColorAlpha(kBlack,alpha)
                hist.SetMarkerStyle(8)
            if 'hist' in datastyle.lower():
                hist.SetFillColorAlpha(0,0)
            
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)

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
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.35, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 1, 0.35))

                else:
                    mains.append(TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.1, 1, 1))
                    subs.append(TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 0, 0))

                if len(signals) == 0:
                    nsignals = 0
                elif addSignals:
                    nsignals = 1
                else:
                    nsignals = len(signals[0])
                legend_topY = 0.73-0.03*(min(len(bkglist[0]),6)+nsignals+1)
                # legend_bottomY = 0.2+0.02*(len(bkglist[0])+nsignals+1)

                legends.append(TLegend(0.65,legend_topY,0.90,0.88))
                legend_duplicates = []
                if not dataOff: legends[hist_index].AddEntry(hist,dataName,datastyle)

                stacks.append(THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
                if totalBkg == None:
                    tot_hist = hist.Clone(hist.GetName()+'_tot')
                    tot_hist.Reset()
                else:
                    tot_hist = totalBkg[hist_index]

                tot_hist.SetTitle(hist.GetName()+'_tot')
                tot_hist.SetMarkerStyle(0)
                tot_hists.append(tot_hist)
                tot_hists_err.append(tot_hist.Clone())
                tot_hists[hist_index].SetLineColor(kBlack)
                tot_hists_err[hist_index].SetLineColor(kBlack)
                tot_hists_err[hist_index].SetLineWidth(0)
                tot_hists_err[hist_index].SetFillColor(kBlack)
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

                # If logy, put QCD on top
                # if logy: bkglist[0], bkglist[-1] = bkglist[-1], bkglist[0]

                # Order based on colors
                if color_idx_order == None:
                    color_idx_order = ColorCodeSortedIndices(colors)
                    colors = [colors[i] for i in color_idx_order]
                
                bkglist[hist_index] = [bkglist[hist_index][i] for i in color_idx_order]
                if bkgNames != [] and isinstance(bkgNames[0],list):
                    bkgNames[hist_index] = [bkgNames[hist_index][i] for i in color_idx_order]

                # Build the stack
                legend_info = collections.OrderedDict()
                for bkg_index,bkg in enumerate(bkglist[hist_index]):     # Won't loop if bkglist is empty
                    # bkg.Sumw2()
                    if totalBkg == None: tot_hists[hist_index].Add(bkg)
                    
                    if logy:
                        bkg.SetMinimum(1e-3)

                    if 'qcd' in bkg.GetName():
                        bkg.SetFillColor(kYellow)
                        bkg.SetLineColor(kYellow)
                    else:
                        if colors[bkg_index] != None:
                            bkg.SetFillColor(colors[bkg_index])
                            bkg.SetLineColor(colors[bkg_index] if colors[bkg_index]!=0 else kBlack)
                        else:
                            bkg.SetFillColor(default_colors[bkg_index])
                            bkg.SetLineColor(default_colors[bkg_index] if colors[bkg_index]!=0 else kBlack)

                    stacks[hist_index].Add(bkg)
                    if bkgNames == []: this_bkg_name = bkg.GetName().split('_')[0]
                    elif type(bkgNames[0]) != list: this_bkg_name = bkgNames[bkg_index]
                    else: this_bkg_name = bkgNames[hist_index][bkg_index]
                    
                    legend_info[this_bkg_name] = bkg

                # Deal with legend which needs ordering reversed from stack build
                for bname in reversed(legend_info.keys()):
                    if bname not in legend_duplicates:
                        legends[hist_index].AddEntry(legend_info[bname],bname,'f')
                        legend_duplicates.append(bname)
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()

                # Set y max of all hists to be the same to accommodate the tallest
                histList = [stacks[hist_index],tot_hists[hist_index],hist]

                yMax = histList[0].GetMaximum()
                for h in range(1,len(histList)):
                    if histList[h].GetMaximum() > yMax:
                        yMax = histList[h].GetMaximum()
                for h in histList:
                    h.SetMaximum(yMax*(2.5-legend_topY+0.03))
                    if logy == True:
                        h.SetMaximum(yMax*10**(2.5-legend_topY+0.1))

                
                mLS = 0.08
                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.SetTitleOffset(1.15,"xy")
                hist.GetYaxis().SetTitleOffset(1.04)
                hist.GetYaxis().SetLabelSize(0.07)
                hist.GetYaxis().SetTitleSize(titleSize)
                hist.GetXaxis().SetLabelSize(0.07)
                hist.GetXaxis().SetTitleSize(titleSize)
                hist.GetXaxis().SetLabelOffset(0.05)
                if logy == True:
                    hist.SetMinimum(1e-3)
                
                hist.GetYaxis().SetNdivisions(508)

                hist.Draw(datastyle+' X0')
                #gStyle.SetErrorX(0.5)

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
                        sig.SetLineColor(kBlack)
                        sig.SetLineWidth(2)
                        if logy == True:
                            sig.SetMinimum(1e-3)
                        # if signalNames == []: this_sig_name = sig.GetName().split('_')[0]
                        if type(signalNames) == str: this_sig_name = signalNames
                        else: this_sig_name = signalNames[isig]

                        legends[hist_index].AddEntry(sig,this_sig_name,'L')
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
                    #gStyle.SetErrorX(0)
                    hist.Draw(datastyle+' X0 same')
                    #gStyle.SetErrorX(0.5)

                gPad.RedrawAxis()

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull
                pulls.append(Make_Pull_plot(hist,tot_hists[hist_index]))
                pulls[hist_index].SetFillColor(kBlue)
                pulls[hist_index].SetTitle(";"+hist.GetXaxis().GetTitle()+";(Data-Bkg)/#sigma")
                pulls[hist_index].SetStats(0)

                LS = .16

                pulls[hist_index].GetYaxis().SetRangeUser(-2.9,2.9)
                pulls[hist_index].GetYaxis().SetTitleOffset(0.4)
                # pulls[hist_index].GetXaxis().SetTitleOffset(0.9)
                             
                pulls[hist_index].GetYaxis().SetLabelSize(0.13)
                pulls[hist_index].GetYaxis().SetTitleSize(0.12)
                pulls[hist_index].GetYaxis().SetNdivisions(306)
                pulls[hist_index].GetXaxis().SetLabelSize(0.13)
                pulls[hist_index].GetXaxis().SetTitleSize(0.15)

                pulls[hist_index].GetXaxis().SetTitle(xtitle)
                pulls[hist_index].GetYaxis().SetTitle("(Data-Bkg)/#sigma")
                pulls[hist_index].Draw('hist')

                if logy == True:
                    mains[hist_index].SetLogy()

                CMS_lumi.extraText = extraText#'Preliminary'
                CMS_lumi.cmsTextSize = 0.9
                CMS_lumi.cmsTextOffset = 2
                CMS_lumi.lumiTextSize = 0.9
                
                CMS_lumi.CMS_lumi(mains[hist_index], year, 11)
                mains[hist_index].cd()
                lumiE = TLatex()
                lumiE.SetNDC()
                lumiE.SetTextAngle(0)
                lumiE.SetTextColor(kBlack)
                lumiE.SetTextFont(42)
                lumiE.SetTextAlign(31) 
                lumiE.SetTextSize(0.7*0.1)
                lumiE.DrawLatex(1-0.05,1-0.1+0.2*0.1,"137 fb^{-1} (13 TeV)")
                
                if isinstance(subtitles,list) and len(subtitles) > 0:
                    subtitle = TLatex()
                    subtitle.SetNDC()
                    subtitle.SetTextAngle(0)
                    subtitle.SetTextColor(kBlack)
                    subtitle.SetTextFont(42)
                    subtitle.SetTextAlign(12) 
                    subtitle.SetTextSize(0.06)
                    # print (subtitles[hist_index])
                    subtitle_string = '%s < %s < %s %s'%(subtitles[hist_index].split('-')[0], sliceVar.split(' ')[0], subtitles[hist_index].split('-')[1], 'GeV')
                    subtitle.DrawLatex(0.208,0.74,subtitle_string)
    # CMS_lumi.CMS_lumi(myCan, year, 11)

    if rootfile:
        myCan.Print(tag+'/'+name+'.root','root')
    else:
        myCan.Print(tag+'/'+name+'.pdf','pdf')
        myCan.Print(tag+'/'+name+'.png','png')


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
    out = TH2D('correlation_matrix','correlation_matrix',nFinalParams,0,nFinalParams,nFinalParams,0,nFinalParams)
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
    BKGUP, BKGDOWN = Make_up_down(BKG)
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
            sigma = sqrt(FSerr*FSerr + BKGerr*BKGerr)
        else:
            sigma = sqrt(BKGerr*BKGerr)
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