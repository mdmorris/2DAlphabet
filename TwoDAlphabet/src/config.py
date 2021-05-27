import itertools, ROOT, argparse
from TwoDAlphabet.src.helpers import nested_dict, open_json, parse_arg_dict, set_hist_maximums
from TwoDAlphabet.src.binning import Binning

class Config:
    '''Class to handle the reading and manipulation of data provided 
    in 2DAlphabet JSON configuration files.
    '''
    def __init__(self,json,findreplace={},externalOptions={}):
        '''Initialize a Config object for a given JSON file and perform
        all initial checks and manipulations.

        @param json (str): File name and path.
        @param findreplace (dict, optional): Find-replace pairs. Defaults to {}.
        @param externalOptions (dict, optional): Extra key-value pairs to add to the OPTIONS
            section of the JSON file. Defaults to {}.
        '''
        self.config = open_json(json)
        self._addFindReplace(findreplace)
        self._varReplacement()
        self.name = self.config['NAME']
        self.options = self.GetOptions(externalOptions)

    def Construct(self):
        '''Uses the config as instructions to setup binning, manipulate histograms,
        save outputs, etc. Separated from the constructor so information from the config
        can be extracted for basic sanity checks before it is processed.
        '''
        self.binning = Binning(self.config['BINNING'])
        self.processes = self.Section('PROCESSES')
        self.systematics = self.Section('SYSTEMATICS')
        self.hists = self.organize_hists() # takes config info, makes a TFile with new hists, returns map to the file contents

    def Section(self,key):
        return {k:v for k,v in self.config[key].items() if k != 'HELP'}

    def _addFindReplace(self,findreplace):
        '''Add external find-replace pairs to the "GLOBAL" entry of config.

        @param findreplace (dict): Find-replace pairs in non-nested dictionary.

        Raises:
            ValueError: If a "find" already exists in self.config['GLOBAL'].
        '''
        for s in self.findreplace.keys():
            if s in self.config['GLOBAL'].keys():
                raise ValueError('A command line string replacement (%s) conflicts with one already in the configuration file. Quitting...' %(s))
            self.config['GLOBAL'][s] = self.findreplace[s]

    def _varReplacement(self):
        '''Do string substitution for config entries based on the dictionary of find
        and replace values found in self.config['GLOBAL']. Call self._addFindReplace()
        before running this function to add in external find-replace pairs.

        Args:
            @param findreplace (dict): Non-nested dictionary with key-value pairs to find-replace
                in the internal class configuration dict.

        Returns:
            None.
        '''
        print ("Doing GLOBAL variable replacement in input json...")
        for old_string in self.inputConfig['GLOBAL'].keys():
            if old_string == "HELP":
                print ('WARNING: The HELP entry is deprecated and checking for it will be removed in the future. Please remove it from your config.')
                continue
            new_string = self.inputConfig['GLOBAL'][old_string]
            self.config = self.config_loop_replace(self.config, old_string, new_string)

    def GetOptions(self,externalOpts={}):
        '''Loads options specific to this config file (from 'OPTIONS' section).
        External options (externalOpts) can be provided but will overwrite those in
        the config (and modify the config in-place so that the version later saved
        reflects the conditions under which the config was used).

        @param externalOpts (dict, optional): Option-value pairs. Defaults to {}.

        Returns:
            ArgumentParser.Namespace
        '''
        self.config['OPTIONS'].update(externalOpts)
        parser = argparse.ArgumentParser()

        parser.add_argument('blindedPlots', default=True, type=bool,
            help='Blind plots to SIG region. Does not blind fit. Defaults to True.')
        parser.add_argument('blindedFit', default=True, type=bool,
            help='Blind fit to SIG region. Does not blind plots. Defaults to True.')
        parser.add_argument('haddSignals', default=True, type=bool,
            help='Combine signals into one histogram for the sake of plotting. Still treated as separate in fit. Defaults to True.')

        parser.add_argument('plotTitles', default=False, type=bool,
            help='Include titles in plots. Defaults to False.')
        parser.add_argument('freezeFail', default=False, type=bool,
            help='Freeze the QCD failing bins to their pre-fit values. Defaults to False.')
        # parser.add_argument('fitGuesses', default=False, type=bool,
        #     help='Save settings to file in json format. Ignored in json file. Defaults to False.')
        parser.add_argument('plotUncerts', default=False, type=bool,
            help='Plot comparison of pre-fit uncertainty shape templates in 1D projections. Defaults to False.')
        # parser.add_argument('prerun', default=False, type=bool,
        #     help='Run a fit of just this config by itself to establish. Defaults to False.')
        # parser.add_argument('rpfRatio', default={}, type=dict,
        #     help='. Defaults to empty dict.')
        parser.add_argument('plotPrefitSigInFitB', default=False, type=bool,
            help='In the b-only post-fit plots, plot the signal normalized to its pre-fit value. Defaults to False.')
        # parser.add_argument('parametricFail', default=False, type=bool,
        #     help='. Defaults to False.')
        parser.add_argument('plotEvtsPerUnit', default=False, type=bool,
            help='Post-fit bins are plotted as events per unit rather than events per bin. Defaults to False.')
        
        parser.add_argument('ySlices', default=[], type=list,
            help='Manually define the slices in the y-axis for the sake of plotting. Only needed if the automated algorithm is not working as intended. Defaults to empty list.')
        parser.add_argument('year', default=1, type=int,
            help='Year information used for the sake of plotting text. Defaults to 1 which indicates that the full Run 2 is being analyzed.')
        parser.add_argument('recycle', default=[], type=list,
            help='List of items to recycle. Not currently working.')
        
        return parse_arg_dict(parser,self.config['OPTIONS'])

    def MakeSystematicPlots(self):
        '''Make plots of the systematic shape variations of each process based on those
        processes and systematic shapes specified in the config. Shapes are presented 
        as projections onto 1D axis where no selection has been made on the axis not
        being plotted. Plots are saved to UncertPlots/.
        '''
        regions = ['pass','fail']
        variations = ['nom','up','down']
        axes = ['X','Y']
        for proc in self.Section('PROCESS'):
            for syst in self.processes[proc]['SYSTEMATICS']:
                if self.systematics[syst]['CODE'] < 2: continue

                tracking_dict = nested_dict(3,None)
                for r,v,x in itertools.product(regions,variations,axes):
                    if v == 'nom': h = self.orgFile.Get(self.organizedDict[proc][r+'_FULL']['nominal'])
                    else:          h = self.orgFile.Get(self.organizedDict[proc][r+'_FULL'][syst+v.capitalize()])
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
                    
                    if x == 'X': nom.SetXTitle(self.inputConfig['BINNING']['X']['TITLE'])
                    elif x == 'Y': nom.SetXTitle(self.inputConfig['BINNING']['Y']['TITLE'])

                    nom.SetTitle('')
                    nom.GetXaxis().SetTitleOffset(1.0)
                    nom.GetXaxis().SetTitleSize(0.05)
                    thisCan.SetRightMargin(0.16)

                    nom.Draw('hist'); up.Draw('same hist'); down.Draw('same hist')
                    thisCan.Print(self.projPath+'/UncertPlots/Uncertainty_'+proc+'_'+syst+r+x+'.png','png')

    def _inputOrganizer(self):
        #################################################################################
        # First we need to get histograms from files and store them in a new dictionary #
        #################################################################################
        dict_hists = {}

        # Stores [process,cat] pairs of regions with integral of zero so we can tell the card this
        self.integralZero = []

        # Grab all process names and loop through
        processes = [process for process in self.inputConfig['PROCESS'].keys() if process != "HELP"]
        if self.rpfRatio != False: processes.append('qcdmc')

        for process in processes:
            if process == 'qcdmc': this_process_dict = self.inputConfig['OPTIONS']['rpfRatio']
            else: this_process_dict = self.inputConfig['PROCESS'][process]
            
            dict_hists[process] = {  
                'file': 0,
                'pass': {},
                'fail': {}
            }

            # Grab nominal pass and fail distributions
            file_nominal = TFile.Open(this_process_dict['FILE'])
            hist_pass = file_nominal.Get(this_process_dict['HISTPASS'])
            hist_fail = file_nominal.Get(this_process_dict['HISTFAIL'])

            # DOCUMENT
            # Flat scale
            if "SCALE" in this_process_dict.keys():
                this_proc_scale = this_process_dict["SCALE"]
                hist_pass.Scale(this_proc_scale)
                hist_fail.Scale(this_proc_scale)
            # Scale by another hist or function
            elif "SCALEPASS" in this_process_dict.keys() and "SCALEFAIL" in this_process_dict.keys():
                this_scale_pass_file = TFile.Open(this_process_dict["SCALEPASS"])
                this_scale_fail_file = TFile.Open(this_process_dict["SCALEFAIL"])
                
                this_proc_scale_pass = this_scale_pass_file.Get(this_process_dict["SCALEPASS_HISTNAME"])
                this_proc_scale_fail = this_scale_fail_file.Get(this_process_dict["SCALEFAIL_HISTNAME"])

                hist_pass.Multiply(this_proc_scale_pass)
                hist_fail.Multiply(this_proc_scale_fail)

            # Smooth
            if 'SMOOTH' in this_process_dict.keys() and this_process_dict['SMOOTH']: smooth_this = True# and self.inputConfig['OPTIONS']['rpfRatio'] != False and self.inputConfig['OPTIONS']['rpfRatio']['SMOOTH']: smooth_this = True
            else: smooth_this = False

            if smooth_this:
                if process != 'qcdmc':
                    hist_pass.Smooth(1,"k5a") #= header.smoothHist2D('smooth_'+process+'_pass',hist_pass,renormalize=False,iterate = 1 if process != 'qcdmc' else 1)
                    hist_fail.Smooth(1,"k5a") #= header.smoothHist2D('smooth_'+process+'_fail',hist_fail,renormalize=False,iterate = 1 if process != 'qcdmc' else 1)

            dict_hists[process]['file'] = file_nominal
            dict_hists[process]['pass']['nominal'] = hist_pass
            dict_hists[process]['fail']['nominal'] = hist_fail

            # If there are systematics
            if 'SYSTEMATICS' not in this_process_dict.keys() or len(this_process_dict['SYSTEMATICS']) == 0:
                print 'No systematics for process ' + process
            else:
                # Loop through them and grab info from inputConfig['SYSTEMATIC']
                for syst in this_process_dict['SYSTEMATICS']:
                    try:
                        this_syst_dict = self.inputConfig['SYSTEMATIC'][syst]

                    # Quit if syst does not exist and user does not want to skip
                    except:
                        skip = raw_input('No entry named "' + syst + '" exists in the SYSTEMATIC section of the input JSON. Skip it? (y/n)')
                        if skip == 'y' or skip == 'Y':
                            print 'Skipping ' + syst
                        else: 
                            print 'Quiting'
                            quit()

                    # Handle case where pass and fail are uncorrelated
                    if 'UNCORRELATED' in this_syst_dict and this_syst_dict['UNCORRELATED']:
                        pass_syst = syst+'_pass'
                        fail_syst = syst+'_fail'
                    else:
                        pass_syst = syst
                        fail_syst = syst

                    # Only care about syst if it's a shape (CODE == 2 or 3)
                    if this_syst_dict['CODE'] == 2:   # same file as norm, different hist names

                        dict_hists[process]['pass'][pass_syst+'Up']   = file_nominal.Get(this_syst_dict['HISTPASS_UP'].replace('*',process))
                        dict_hists[process]['pass'][pass_syst+'Down'] = file_nominal.Get(this_syst_dict['HISTPASS_DOWN'].replace('*',process))
                        dict_hists[process]['fail'][fail_syst+'Up']   = file_nominal.Get(this_syst_dict['HISTFAIL_UP'].replace('*',process))
                        dict_hists[process]['fail'][fail_syst+'Down'] = file_nominal.Get(this_syst_dict['HISTFAIL_DOWN'].replace('*',process))

                    if this_syst_dict['CODE'] == 3:   # different file as norm and different files for each process if specified, same hist name if not specified in inputConfig
                        # User will most likely have different file for each process but maybe not so check
                        if 'FILE_UP' in this_syst_dict:
                            file_up = TFile.Open(this_syst_dict['FILE_UP'])
                        # Wild card to replace * with the process name
                        elif 'FILE_UP_*' in this_syst_dict:
                            file_up = TFile.Open(this_syst_dict['FILE_UP_*'].replace('*',process))
                        else:
                            file_up = TFile.Open(this_syst_dict['FILE_UP_'+process])

                        if 'FILE_DOWN' in this_syst_dict:
                            file_down = TFile.Open(this_syst_dict['FILE_DOWN'])
                        elif 'FILE_DOWN_*' in this_syst_dict:
                            file_down = TFile.Open(this_syst_dict['FILE_DOWN_*'].replace('*',process))
                        else:
                            file_down = TFile.Open(this_syst_dict['FILE_DOWN_'+process])

                        dict_hists[process]['file_'+syst+'Up'] = file_up
                        dict_hists[process]['file_'+syst+'Down'] = file_down

                        if 'HISTPASS_UP' in this_syst_dict:
                            dict_hists[process]['pass'][pass_syst+'Up'] = file_up.Get(this_syst_dict['HISTPASS_UP'])            # try to grab hist name from SYSTEMATIC dictionary
                        elif 'HISTPASS' in this_syst_dict:
                            dict_hists[process]['pass'][pass_syst+'Up'] = file_up.Get(this_syst_dict['HISTPASS'])               # else use the same one as nominal distribution
                        elif 'HISTPASS_UP_*' in this_syst_dict:
                            dict_hists[process]['pass'][pass_syst+'Up'] = file_up.Get(this_syst_dict['HISTPASS_UP_*'].replace('*',process))
                        else: 
                            dict_hists[process]['pass'][pass_syst+'Up'] = file_up.Get(this_syst_dict['HISTPASS_UP_'+process])   # or use process specific name

                        if 'HISTPASS_DOWN' in this_syst_dict:
                            dict_hists[process]['pass'][pass_syst+'Down'] = file_down.Get(this_syst_dict['HISTPASS_DOWN'])
                        elif 'HISTPASS' in this_syst_dict:
                            dict_hists[process]['pass'][pass_syst+'Down'] = file_down.Get(this_syst_dict['HISTPASS'])
                        elif 'HISTPASS_DOWN_*' in this_syst_dict:
                            dict_hists[process]['pass'][pass_syst+'Down'] = file_up.Get(this_syst_dict['HISTPASS_DOWN_*'].replace('*',process))
                        else:
                            dict_hists[process]['pass'][pass_syst+'Down'] = file_down.Get(this_syst_dict['HISTPASS_DOWN_' + process])

                        if 'HISTFAIL_UP' in this_syst_dict:
                            dict_hists[process]['fail'][fail_syst+'Up'] = file_up.Get(this_syst_dict['HISTFAIL_UP'])
                        elif 'HISTFAIL' in this_syst_dict:
                            dict_hists[process]['fail'][fail_syst+'Up'] = file_up.Get(this_syst_dict['HISTFAIL'])
                        elif 'HISTFAIL_UP_*' in this_syst_dict:
                            dict_hists[process]['fail'][fail_syst+'Up'] = file_up.Get(this_syst_dict['HISTFAIL_UP_*'].replace('*',process))    
                        else:
                            dict_hists[process]['fail'][fail_syst+'Up'] = file_up.Get(this_syst_dict['HISTFAIL_UP_' + process])

                        if 'HISTFAIL_DOWN' in this_syst_dict:
                            dict_hists[process]['fail'][fail_syst+'Down'] = file_down.Get(this_syst_dict['HISTFAIL_DOWN'])
                        elif 'HISTFAIL' in this_syst_dict:
                            dict_hists[process]['fail'][fail_syst+'Down'] = file_down.Get(this_syst_dict['HISTFAIL'])
                        elif 'HISTFAIL_DOWN_*' in this_syst_dict:
                            dict_hists[process]['fail'][fail_syst+'Down'] = file_up.Get(this_syst_dict['HISTFAIL_DOWN_*'].replace('*',process))
                        else:
                            dict_hists[process]['fail'][fail_syst+'Down'] = file_down.Get(this_syst_dict['HISTFAIL_DOWN_' + process])

                    if this_syst_dict['CODE'] > 1:
                        if smooth_this:
                            dict_hists[process]['pass'][pass_syst+'Up'].Smooth(1,"k5a") #= header.smoothHist2D('smooth_'+process+'_pass_'+syst+'Up',dict_hists[process]['pass'][pass_syst+'Up'],renormalize=False)
                            dict_hists[process]['pass'][pass_syst+'Down'].Smooth(1,"k5a") #= header.smoothHist2D('smooth_'+process+'_pass_'+syst+'Down',dict_hists[process]['pass'][pass_syst+'Down'],renormalize=False)
                            dict_hists[process]['fail'][fail_syst+'Up'].Smooth(1,"k5a")   #= header.smoothHist2D('smooth_'+process+'_fail_'+syst+'Up',dict_hists[process]['fail'][fail_syst+'Up'],renormalize=False)
                            dict_hists[process]['fail'][fail_syst+'Down'].Smooth(1,"k5a") #= header.smoothHist2D('smooth_'+process+'_fail_'+syst+'Down',dict_hists[process]['fail'][fail_syst+'Down'],renormalize=False)

                        if "SCALE" in this_process_dict.keys():
                            dict_hists[process]['pass'][pass_syst+'Up'].Scale(this_proc_scale)
                            dict_hists[process]['pass'][pass_syst+'Down'].Scale(this_proc_scale)
                            dict_hists[process]['fail'][fail_syst+'Up'].Scale(this_proc_scale)
                            dict_hists[process]['fail'][fail_syst+'Down'].Scale(this_proc_scale)
                        elif "SCALEPASS" in this_process_dict.keys() and "SCALEFAIL" in this_process_dict.keys():
                            dict_hists[process]['pass'][pass_syst+'Up'].Multiply(this_proc_scale_pass)
                            dict_hists[process]['pass'][pass_syst+'Down'].Multiply(this_proc_scale_pass)
                            dict_hists[process]['fail'][fail_syst+'Up'].Multiply(this_proc_scale_fail)
                            dict_hists[process]['fail'][fail_syst+'Down'].Multiply(this_proc_scale_fail)

        #####################################################################
        # With dictionary made, we can split around the signal region and   #
        # start renaming to match the format required by Combine. The       #
        # dictionary key names are conveniently named so we can do this     #
        # with minimal pain.                                                #
        #####################################################################
        temp_TH2 = dict_hists['data_obs']['pass']['nominal']
        old_x_min = temp_TH2.GetXaxis().GetXmin()
        old_x_max = temp_TH2.GetXaxis().GetXmax()
        # old_x_nbins = temp_TH2.GetNbinsX()
        # old_x_width = float(old_x_max-old_x_min)/float(old_x_nbins)
        old_y_min = temp_TH2.GetYaxis().GetXmin()
        old_y_max = temp_TH2.GetYaxis().GetXmax()
        old_y_nbins = temp_TH2.GetNbinsY()
        old_y_width = float(old_y_max-old_y_min)/float(old_y_nbins)

        # Print out info
        print "Applying new Y bins: ["+str(old_y_min)+","+str(old_y_max)+"] -> ["+str(self.newYbins[0])+","+str(self.newYbins[-1])+"]"
        print 'Applying new X bins: '
        for c in ['LOW','SIG','HIGH']: 
            print '\t'+c + ': ['+str(old_x_min)+","+str(old_x_max)+"] -> ["+str(self.newXbins[c][0])+","+str(self.newXbins[c][-1])+"]"

        self.orgFile.cd()

        # For each process, category, and dist (nominal, systUp, etc)
        for process in dict_hists.keys():
            self.organizedDict[process] = {'pass_FULL':{}, 'pass_LOW':{}, 'pass_SIG':{}, 'pass_HIGH':{}, 'fail_FULL':{}, 'fail_LOW':{}, 'fail_SIG':{}, 'fail_HIGH':{}}
            for cat in ['pass','fail']:
                for dist in dict_hists[process][cat].keys():
                    print 'Making ' + process +', ' + cat + ', ' + dist

                    # Get new names
                    temp_histname = process + '_' + cat
                    if dist != 'nominal':                           # if not nominal dist
                        temp_histname = temp_histname + '_' + dist

                    # If there are user specified y bins...
                    if self.newYbins != False:
                        temp_hist = header.copyHistWithNewYbins(dict_hists[process][cat][dist],self.newYbins,temp_histname)#,self.oldYwidth)
                    else:
                        temp_hist = dict_hists[process][cat][dist]

                    # If there are user specified x bins...
                    for c in ['FULL','LOW','SIG','HIGH']: 
                        # Get new names
                        histname = process + '_' + cat+'_'+c+'_'+self.name
                        if dist != 'nominal':                           # if not nominal dist
                            histname = histname + '_' + dist
                        print 'Making '+histname
                        if c != 'FULL': finalhist = header.copyHistWithNewXbins(temp_hist,self.newXbins[c],histname)#,self.oldYwidth)
                        else: finalhist = header.copyHistWithNewXbins(temp_hist,self.fullXbins,histname)

                        # Test if histogram is non-zero
                        if finalhist.Integral() <= 0:
                            print 'WARNING: '+process+', '+cat+'_'+c+', '+dist+' has zero or negative events - ' + str(finalhist.Integral())
                            self.integralZero.append([process,cat+'_'+c])
                            # If it is, zero the bins except one to avoid Integral=0 errors in combine
                            for b in range(1,finalhist.GetNbinsX()*finalhist.GetNbinsY()+1):
                                finalhist.SetBinContent(b,1e-10)

                        finalhist.Write()
                        self.organizedDict[process][cat+'_'+c][dist] = finalhist.GetName()#header.copyHistWithNewXbins(temp_hist,self.newXbins[c],histname)


    def _makeCard(self):
        # Recreate file
        card_new = open(self.projPath + 'card_'+self.name+'.txt','w')

        column_width = 11+len(self.name)

        #######################################################
        # imax (bins), jmax (backgrounds), kmax (systematics) #
        #######################################################
        imax = '6'                      # pass, fail for each 'X' axis category
        channels = []
        for r in ['pass', 'fail']:
            for c in ['LOW','SIG','HIGH']:
                channels.append(r+'_'+c+'_'+self.name)                

        # Get the length of the list of all process that have CODE 2 (and ignore "HELP" key) and add 1 for qcd (which won't be in the inputConfig)
        jmax = str(len([proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and self.inputConfig['PROCESS'][proc]['CODE'] == 2]) + 1)
        # Get the length of the lsit of all systematics (and ignore "HELP" key)
        n_uncorr_systs = len([syst for syst in self.inputConfig['SYSTEMATIC'].keys() if syst != 'HELP' and 'UNCORRELATED' in self.inputConfig['SYSTEMATIC'][syst] and self.inputConfig['SYSTEMATIC'][syst]['UNCORRELATED']])
        kmax = str(len([syst for syst in self.inputConfig['SYSTEMATIC'].keys() if syst != 'HELP' and syst != 'KDEbandwidth'])+n_uncorr_systs)
        if self.rpfRatioVariations != False: kmax = str(int(kmax)+1)

        card_new.write('imax '+imax+'\n')      
        card_new.write('jmax '+jmax+'\n')
        card_new.write('kmax '+kmax+'\n')
        card_new.write('-'*120+'\n')

        ##########
        # Shapes #
        ##########
        procs_with_systs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and len(self.inputConfig['PROCESS'][proc]['SYSTEMATICS']) != 0]
        procs_without_systs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and len(self.inputConfig['PROCESS'][proc]['SYSTEMATICS']) == 0]
        if self.rpfRatioVariations == False: procs_without_systs.append('qcd')   # Again, qcd not in the input JSON but needs to be in the combine card!
        else: procs_with_systs.append('qcd')

        for proc in procs_without_systs:
            card_new.write(header.colliMate('shapes  '+proc+' * base_'+self.name+'.root w_'+self.name+':'+proc+'_$CHANNEL\n'))
        for proc in procs_with_systs:
            card_new.write(header.colliMate('shapes  '+proc+' * base_'+self.name+'.root w_'+self.name+':'+proc+'_$CHANNEL w_'+self.name+':'+proc+'_$CHANNEL_$SYSTEMATIC\n'))

        card_new.write('-'*120+'\n')

        ####################################
        # Set bin observation values to -1 #
        ####################################
        tempString = 'bin  '
        for chan in channels:
            tempString += (chan+' ')
        tempString += '\n'
        card_new.write(header.colliMate(tempString,column_width))

        tempString = 'observation  '
        for ichan in range(int(imax)):
            tempString += '-1 '
        tempString += '\n'
        card_new.write(header.colliMate(tempString,column_width))

        card_new.write('-'*120+'\n')

        ######################################################
        # Tie processes to bins and rates and simultaneously #
        # create the systematic uncertainty rows             #
        ######################################################
        bin_line = 'bin  '
        processName_line = 'process  '
        processCode_line = 'process  '
        rate_line = 'rate  '
        syst_lines = {}

        # Fill syst_lines with keys to initialized strings
        for syst in [systematic for systematic in self.inputConfig['SYSTEMATIC'].keys() if systematic != 'HELP']:
            # Get type
            if self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 0 or self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 1:        # lnN
                syst_type = 'lnN'
            elif self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 2 or self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 3:      # shape
                syst_type = 'shape'
            else:
                print 'Systematic ' + syst + ' does not have one of the four allowed codes (0,1,2,3). Quitting.'
                quit()
            
            # Build line
            if 'UNCORRELATED' in self.inputConfig['SYSTEMATIC'][syst].keys() and self.inputConfig['SYSTEMATIC'][syst]['UNCORRELATED']:
                syst_lines[syst+'_pass'] = syst + '_pass ' + syst_type + ' '
                syst_lines[syst+'_fail'] = syst + '_fail ' + syst_type + ' '
            elif syst == 'KDEbandwidth':
                continue
            else:
                syst_lines[syst] = syst + ' ' + syst_type + ' '

        if self.rpfRatioVariations != False:
            syst_lines['KDEbandwidth'] = 'KDEbandwidth shape '

        signal_procs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs' and self.inputConfig['PROCESS'][proc]['CODE'] == 0]
        MC_bkg_procs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs' and (self.inputConfig['PROCESS'][proc]['CODE'] == 2 or self.inputConfig['PROCESS'][proc]['CODE'] == 3)]

        all_procs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs']
        all_procs.append('qcd')

        for chan in channels:
            for proc in all_procs:
                # Start lines
                bin_line += (chan+' ')
                processName_line += (proc+' ')

                # If signal
                if proc in signal_procs:
                    processCode_line += (str(0-signal_procs.index(proc))+' ')
                    rate_line += ('-1 ')

                # If bkg
                elif proc in MC_bkg_procs:
                    processCode_line += (str(MC_bkg_procs.index(proc)+1)+' ')
                    if self.inputConfig['PROCESS'][proc]['CODE'] == 2:       # No floating normalization
                        rate_line += '-1 '                                            

                # If qcd
                elif proc == 'qcd':
                    processCode_line += (str(len(MC_bkg_procs)+2)+' ')
                    rate_line += '1 '

                # Fill systematic lines
                for syst_line_key in syst_lines.keys():
                    # Check for case when pass and fail are uncorrelated
                    if syst_line_key.split('_')[-1] in ['pass','fail']:
                        chan_specific = syst_line_key.split('_')[-1] 
                    else:
                        chan_specific = False

                    # If the systematic is applicable to the process
                    if proc != 'qcd':
                        base_syst_line_key = syst_line_key.replace('_pass','').replace('_fail','')
                        if base_syst_line_key in self.inputConfig['PROCESS'][proc]['SYSTEMATICS']:
                            # If we have the pass(fail) specific systematic and this is a fail(pass) region, go to next and skip the rest below
                            if chan_specific != False and chan_specific not in chan: 
                                thisVal = '-'

                            # If symmetric lnN...
                            elif self.inputConfig['SYSTEMATIC'][base_syst_line_key]['CODE'] == 0:
                                thisVal = str(self.inputConfig['SYSTEMATIC'][base_syst_line_key]['VAL'])
                            # If asymmetric lnN...
                            elif self.inputConfig['SYSTEMATIC'][base_syst_line_key]['CODE'] == 1:
                                thisVal = str(self.inputConfig['SYSTEMATIC'][base_syst_line_key]['VALDOWN']) + '/' + str(self.inputConfig['SYSTEMATIC'][base_syst_line_key]['VALUP'])
                            # If shape...
                            else:
                                if 'SCALE' not in self.inputConfig['SYSTEMATIC'][base_syst_line_key]:
                                    thisVal = '1.0'
                                else:
                                    thisVal = str(self.inputConfig['SYSTEMATIC'][base_syst_line_key]['SCALE'])
                        # Otherwise place a '-'
                        else:
                            thisVal = '-'  

                    elif syst_line_key == 'KDEbandwidth' and proc == 'qcd' and 'pass' in chan:
                        if 'SCALE' not in self.inputConfig['SYSTEMATIC'][syst_line_key]:
                            thisVal = '1.0'
                        else:
                            thisVal = self.inputConfig['SYSTEMATIC'][syst_line_key]['SCALE']
                    else:
                        thisVal = '-'

                    syst_lines[syst_line_key] += (thisVal+' ')

        card_new.write(header.colliMate(bin_line+'\n',column_width))
        card_new.write(header.colliMate(processName_line+'\n',column_width))
        card_new.write(header.colliMate(processCode_line+'\n',column_width))
        card_new.write(header.colliMate(rate_line+'\n',column_width))
        card_new.write('-'*120+'\n')

        ############################
        # Systematic uncertainties #
        ############################
        for line_key in syst_lines.keys():
            card_new.write(header.colliMate(syst_lines[line_key]+'\n',column_width))

        ######################################################
        # Mark floating values as flatParams                 # 
        # We float just the rpf params and the failing bins. #
        ######################################################
        for p in self.rpf.getFuncVarNames():
            card_new.write(header.colliMate(self.rpf.funcVars[p].GetName()+' flatParam\n',22))

        if self.fail_func != False:
            for p in self.fail_func.getFuncVarNames():
                card_new.write(header.colliMate(self.fail_func.funcVars[p].GetName()+' flatParam\n',22))

        for b in self.floatingBins:
            card_new.write(header.colliMate(b+' flatParam\n',22))
           
        card_new.close() 


def config_loop_replace(config,old,new):
    '''Self-calling function loop to find-replace one pair (old,new)
    in a nested dictionary or list (config). If old, new, and the config entry
    examined are all strings, all instances of <old> will be replaced by <new> in the
    string in the config. If this condition is not satisified, <old> must match the config
    entry in its entirety (ie. <old> == <config dict value>).

    @param config (dict,list): Nested dictionary or list where keys and values will have the
        find-replace algorithm applied.
    @param old (non-nested obj): Object to replace (string, int, float, etc). No lists or dictionaries.
    @param new (non-nested obj): Object replacement (string, int, float, etc). No lists or dictionaries.

    Raises:
        TypeError: If input is not a dict or list.

    Returns:
        dict,list: Modified dict/list.
    '''
    if isinstance(config,dict):
        for k,v in config.items():
            if old in k:
                config[k.replace(old,new)] = config.pop(k)
                k = k.replace(old,new)
            if isinstance(v,dict) or isinstance(v,list):
                config[k] = config_loop_replace(v)
            elif isinstance(v,str) and isinstance(old,str) and isinstance(new,str):
                if old in v:
                    config[k] = v.replace(old,new)
            else:
                if old == v:
                    config[k] = new
    elif isinstance(config,list):
        for i,v in enumerate(config):
            if isinstance(v,dict) or isinstance(v,list):
                config[i] = config_loop_replace(v)
            elif isinstance(v,str):
                if old in v:
                    config[i] = v.replace(old,new)
            else:
                if old == v:
                    config[i] = new
    else:
        raise TypeError('Type "%s" not accepted in config_loop_replace.')

    return config
