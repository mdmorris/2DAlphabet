import ROOT, argparse, json, pickle, os
from TwoDAlphabet.src.helpers import is_filled_list, nested_dict, open_json, parse_arg_dict
from TwoDAlphabet.src.binning import Binning, copy_hist_with_new_bins, get_bins_from_hist

_protected_keys = ["PROCESS","SYSTEMATICS","SYSTEMATIC","BINNING","OPTIONS","GLOBAL","SCALE","COLOR","TYPE","X","Y","NAME","TITLE","BINS","NBINS","LOW","HIGH"]
class Config:
    '''Class to handle the reading and manipulation of data provided 
    in 2DAlphabet JSON configuration files.
    '''
    def __init__(self,json,projPath,findreplace={},externalOptions={},loadPrevious=False):
        '''Initialize a Config object for a given JSON file and perform
        all initial checks and manipulations.

        Args:
            json (str): File name and path.
            projPath (str): Project path.
            findreplace (dict, optional): Find-replace pairs. Defaults to {}.
            externalOptions (dict, optional): Extra key-value pairs to add to the OPTIONS
                section of the JSON file. Defaults to {}.
            loadPrevious (bool, optional): Load the previous histograms made instead of remaking. Defaults to False.
        '''
        self.config = open_json(json)
        self._addFindReplace(findreplace)
        self._varReplacement()
        self.name = self.config['NAME']
        self.projPath = projPath+'/'+self.name+'/'
        self.options = self.GetOptions(externalOptions)
        self.loadPrevious = loadPrevious

    def Construct(self):
        '''Uses the config as instructions to setup binning, manipulate histograms,
        save outputs, etc. Separated from the constructor so information from the config
        can be extracted for basic sanity checks before it is processed.
        '''
        self.binning = Binning(self.config['BINNING'])
        self.processes = self.Section('PROCESSES')
        self.systematics = self.Section('SYSTEMATICS')
        if self.loadPrevious: self.organized_hists = self.ReadIn()
        else:
            self.organized_hists = organize_inputs(self)
            self.organized_hists.CreateSubRegions()

    def Section(self,key):
        '''Derive the dictionary for a given section of the configuration file
        with HELP keys removed.

        Args:
            key (str): Section name (all capital letters) of the configuration file.

        Returns:
            dict: Section of config.
        '''
        return {k:v for k,v in self.config[key].items() if k != 'HELP'}

    def _addFindReplace(self,findreplace):
        '''Add external find-replace pairs to the "GLOBAL" entry of config.

        Args:
            findreplace (dict): Find-replace pairs in non-nested dictionary.

        Raises:
            ValueError: If a "find" already exists in self.config['GLOBAL'].
        '''
        for s in findreplace.keys():
            if s in self.config['GLOBAL'].keys():
                raise ValueError('A command line string replacement (%s) conflicts with one already in the configuration file. Quitting...' %(s))
            self.config['GLOBAL'][s] = findreplace[s]

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

        parser.add_argument('blindedPlots', default=True, type=bool, nargs='?',
            help='Blind plots to SIG region. Does not blind fit. Defaults to True.')
        parser.add_argument('blindedFit', default=True, type=bool, nargs='?',
            help='Blind fit to SIG region. Does not blind plots. Defaults to True.')
        parser.add_argument('haddSignals', default=True, type=bool, nargs='?',
            help='Combine signals into one histogram for the sake of plotting. Still treated as separate in fit. Defaults to True.')

        parser.add_argument('plotTitles', default=False, type=bool, nargs='?',
            help='Include titles in plots. Defaults to False.')
        parser.add_argument('freezeFail', default=False, type=bool, nargs='?',
            help='Freeze the QCD failing bins to their pre-fit values. Defaults to False.')
        parser.add_argument('plotUncerts', default=False, type=bool, nargs='?',
            help='Plot comparison of pre-fit uncertainty shape templates in 1D projections. Defaults to False.')
        parser.add_argument('plotPrefitSigInFitB', default=False, type=bool, nargs='?',
            help='In the b-only post-fit plots, plot the signal normalized to its pre-fit value. Defaults to False.')
        parser.add_argument('plotEvtsPerUnit', default=False, type=bool, nargs='?',
            help='Post-fit bins are plotted as events per unit rather than events per bin. Defaults to False.')
        
        parser.add_argument('ySlices', default=[], type=list, nargs='?',
            help='Manually define the slices in the y-axis for the sake of plotting. Only needed if the automated algorithm is not working as intended. Defaults to empty list.')
        parser.add_argument('year', default=1, type=int, nargs='?',
            help='Year information used for the sake of plotting text. Defaults to 1 which indicates that the full Run 2 is being analyzed.')
        parser.add_argument('recycle', default=[], type=list, nargs='?',
            help='List of items to recycle. Not currently working.')
        
        return parse_arg_dict(parser,self.config['OPTIONS'])

    def SaveOut(self): # pragma: no cover
        '''Save three objects to the <self.projPath> directory:
        - a copy of the manipulated config (runConfig.json)
        - the pickled histogram map (hist_map.p)
        - the RooWorkspace (base_<self.name>.root)
        '''
        file_out = open(self.projPath+'runConfig.json', 'w')
        json.dump(self.inputConfig,file_out,indent=2,sort_keys=True)
        file_out.close()
        pickle.dump(self.organized_hists, open(self.projPath+'hist_map.p','wb'))
        self.workspace.writeToFile(self.projPath+'base_'+self.name+'.root',True)  

    def ReadIn(self): # pragma: no cover
        '''Read the histogram map from the pickled file.
        '''
        self.organized_hists = pickle.load(open(self.projPath+'hist_map.p','rb'))

    def MakeCard(self):
        # Recreate file
        card_new = open(self.projPath+'card_'+self.name+'.txt','w')

        # column_width = 11+len(self.name)

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
                raise NameError('Systematic ' + syst + ' does not have one of the four allowed codes (0,1,2,3). Quitting.')
            
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

class OrganizedHists():
    '''Class to store histograms in a consistent data structure and with accompanying
    methods to store, manipulate, and access the histograms.

    Attributes:
        name (str): Name, taken from input configObj.
        filename (str): Path to `organized_hists.root`.
        hists (dict): Three-level nested dictionary organized as [process][region][systematic variation].
        binning (Binning): Binning object, taken from configObj.
        rebinned (bool): Flag to denote if a rebinning has already occured.
        openOption (str): ROOT TFile::Open option set based on whether the file exists and if the configObj requests an overwrite.
        file (ROOT.TFile): TFile to store histograms on disk.
    '''
    def __init__(self,configObj):
        '''Constructor based on a Config object.

        Args:
            configObj (Config): Config object.
        '''
        ##
        self.name = configObj.name
        self.filename = configObj.projPath + 'organized_hists.root'
        self.hists = nested_dict(3,None) # proc, reg, syst
        self.binning = configObj.binning
        self.rebinned = False

        if os.path.exists(self.filename) and not configObj.options.overwrite:
            self.openOption = "OPEN"
        else:
            self.openOption = "RECREATE"
        self.file = ROOT.TFile.Open(self.filename,self.openOption)

    def Add(self,info):
        '''Add histogram and save information on it. Input is a HistInfo object
        used as a package of relevant info.

        Args:
            info (HistInfo): HistInfo object to add.

        Raises:
            RuntimeError: If openOption == "OPEN" in which case no new histograms can be added since this object is in read-only mode.
        '''
        if self.openOption == "OPEN":
            raise RuntimeError('Cannot add a histogram to OrganizedHists if an existing TFile was opened. Set option "overwrite" to True if you wish to delete the existing file.')
        self.hists[info.process][info.region][info.systematic] = info
        self.file.WriteObject(info.Fetch(self.binning),info.histname)

    def Get(self,histname='',process='',region='',systematic=''):
        '''Get histogram from the opened TFile. Specify the histogram
        you want via `histname` or by the combination of `process`, `region`,
        and `systematic` options. The `histname` option will take priority.

        Args:
            histname (str, optional): Name of histogram to get. Overrides other three options if specified. Defaults to ''.
            process (str, optional): Name of process to search for. Must be used in conjunction with `region` and `systematic` options. Overridden by `histname`. Defaults to ''.
            region (str, optional): Name of region to search for. Must be used in conjunction with `process` and `systematic` options. Overridden by `histname`. Defaults to ''.
            systematic (str, optional): Name of systematic to search for. Must be used in conjunction with `process` and `region` options. Overridden by `histname`. Defaults to ''.

        Returns:
            TH2F: Histogram from file.
        '''
        return self.file.Get(histname if histname != '' else '_'.join(process,region,systematic))

    @property
    def _allHists(self):
        '''List of all stored/tracked histograms.

        Returns:
            list(): [description]
        '''
        return [h.hist for p in self.hists for r in self.hists[p] for h in self.hists[p][r]]

    def CreateSubRegions(self):
        '''Sub-divide all histograms along the X axis into the regions specified in the config
        and add the new histograms to the object for tracking. Sets the attribute `rebinned` to True.

        Raises:
            RuntimeError: If `rebinned` is already True.
        '''
        if self.rebinned: raise RuntimeError('Already rebinned this OrganizedHists object.')
        self.rebinned = True

        for p in self.hists:
            for r in self.hists[p]:
                for hinfo in self.hists[p][r].values():
                    h = hinfo.hist

                    for sub in ['LOW','SIG','HIGH']:
                        this_sub_hinfo = hinfo.Clone(region=hinfo.region+'_'+sub)
                        this_sub_hinfo.hist = copy_hist_with_new_bins(hinfo.histname+'_'+sub,'X',h,self.binning.xbinByCat[sub])
                        self.Add(this_sub_hinfo)

class HistInfo():
    def __init__(self,proc,region,syst,color,scale):
        self.process = proc
        self.region = region
        self.systematic = syst
        self.new_name = '%s_%s_%s'%(proc,region,syst) if syst != 'nominal' else '%s_%s'%(proc,region)
        self.color = color
        self.scale = scale
        self.hist = None

    def Clone(self,process=None,region=None,systematic=None):
        return HistInfo(process if process!=None else self.process,
                        region if region!=None else self.region,
                        systematic if systematic!=None else self.systematic,
                        self.color,self.scale)

    def Fetch(self):
        return self.hist

class InputHistInfo(HistInfo):
    def __init__(self, hname, fname, proc, region, syst, color, scale):
        super().__init__(proc, region, syst, color, scale)
        self.histname = hname
        self.filename = fname

    def Fetch(self,binning=None):
        f = ROOT.TFile.Open(self.filename,'OPEN')
        h = f.Get(self.input_hname)
        h.SetDirectory(0)
        if isinstance(self.scale,float):
            h.Scale(self.scale)
        elif isinstance(self.scale,dict) and (self.region in self.scale):
            if 'FILE' in self.scale:
                scale_file = ROOT.TFile.Open(self.scale['FILE'])
            else:
                scale_file = f
            scale_hist = scale_file.Get(self.scale[self.region])
            h.Multiply(scale_hist)
        f.Close()

        if binning != None:
            if get_bins_from_hist("Y", h) != binning.ybinList:
                h = self.Get(self.histname)
                h = copy_hist_with_new_bins(self.histname+'_rebinY','Y',h,binning.ybinList)
            if get_bins_from_hist("X") != binning.xbinList:
                h = copy_hist_with_new_bins(self.histname+'_FULL','X',h,binning.xbinList)
            else:
                h.SetName(self.histname+'_FULL')

        if h.Integral() <= 0:
            print ('WARNING: %s, %s, %s has zero or negative events - %s'%(self.process,self.region,self.syst,h.Integral()))
            for b in range(1,h.GetNbinsX()*h.GetNbinsY()+1):
                h.SetBinContent(b,1e-10)

        return h

def organize_inputs(configObj):
    hists = OrganizedHists(configObj)

    for process,proc_info in configObj.processes.items():
        if not is_filled_list(proc_info,'SYSTEMATICS'):
            print ('---- No systematics for process ' + process)
            this_proc_shape_systs = []
        else:
            this_proc_shape_systs = [s for s in proc_info['SYSTEMATICS'] if len(set(configObj.systematics[s]) & {'VAL','VALUP','VALDOWN'}) == 0]
        
        for filepath,region_pairs in _parse_file_entries(proc_info):
            # Do Nominal
            for region_name,hist_name in region_pairs:
                hists.Add(HistInfo(hist_name,filepath,process,region_name,'nominal',proc_info['COLOR'],proc_info['SCALE']))
            # Loop over systematics
            for syst in this_proc_shape_systs:
                if syst not in configObj.systematics:
                    raise IndexError('No entry named "%s" exists in the SYSTEMATIC section of the input config.'%syst)
                syst_info = configObj.systematics[syst]
                # Get file and hist names for this syst
                for syst_filepath,syst_region_pairs in _parse_file_entries(syst_info):
                    if syst_filepath == None: syst_filepath = filepath
                    syst_filepath = syst_filepath.replace('*',process)
                    for syst_region_name,syst_hist_name in syst_region_pairs:
                        syst_hist_name = syst_hist_name.replace('*',process)
                        region_name    = syst_region_name.split('_')[:-1]
                        variation      = syst_region_name.split('_')[-1]
                        syst_name      = syst if region_name not in syst_info['UNCORRELATED'] else syst+'_'+region_name

                        if variation == 'UP':      syst_name+='Up'
                        elif variation == 'DOWN':  syst_name+='Down'
                        else:                      raise NameError('Variation for %s, %s is %s but can only be "UP" or "DOWN".'%(process,syst_region_name,variation))

                        hists.Add(InputHistInfo(syst_hist_name,syst_filepath,process,region_name,syst_name+variation,proc_info['COLOR'],proc_info['SCALE']))

    return hists

def _parse_file_entries(d):
    '''Logic to parse sub-dictionaries of the JSON config containing the file and region information. 
    A given entry in the dictionary can have two forms:
    
    - The file name/path is the key and the corresponding value is a new dictionary
    with region names as keys and the sub-values are histogram names inside
    the file specified at the first level [1].
    
    - There is no file name/path which can only happen in the case of looking
    up systematic uncertainty shape variations (in SYSTEMATIC section of config).
    This form assumes one is already nested inside the file level described above.
    Region name as the key and histogram name as the value still applies [2].

    [1] -
    {"path/to/myfile1.root":
        {
            "regA":"histA", "regB":"histB"
        },
     "path/to/myfile2.root":
        {
            "regC":"histC", "regD":"histD"
        }
    }
    [2] - 
    {"regA":"histA",
     "regB":"histB",
     "regC":"histC",
     "regD":"histD"
    }

    @param d (dict): Input dictionary to parse.

    Raises:
        FileExistsError: If file specified by key does not exist.
        ValueError: If there are duplicate regions are specified. Only one histogram per-region can be provided.

    Yields:
        tuple( str, list(tuple(str,str)) ): The first return is the file path (`None` if `d` is of form #2). The second is a list
            of the (region name, histogram name) pairs.
    '''
    file_keys = [f for f in d if (f not in _protected_keys and f.endswith('.root'))]
    if len(file_keys) > 0:
        regions = [k for f in file_keys for k in d[f]]
        if len(regions) != len(set(regions)):
            raise ValueError("There are duplicated regions in this process. Printing config process entry for debug.\n\n%s\n"%{f:r for f,r in d.items() if f in file_keys})

        for f in file_keys:
            if not os.path.exists(f):
                raise FileExistsError("Cannot access file %s"%f)
            
            yield f,[(rname,d[f][rname]) for rname in d[f]]
    else:
        region_pairs = [(rname,d[rname]) for rname in d if (rname not in _protected_keys)]
        yield None,region_pairs

def config_loop_replace(config,old,new):
    '''Self-calling function loop to find-replace one pair (old,new)
    in a nested dictionary or list (config). If old, new, and the config entry
    examined are all strings, all instances of <old> will be replaced by <new> in the
    string in the config. If this condition is not satisified, <old> must match the config
    entry in its entirety (ie. <old> == <config dict value>).

    @param config (dict,list): Nested dictionary or list where keys and values will have the
        find-replace algorithm applied.
    @param old (non-nested obj): Object to replace (of type string, int, float, etc - no lists or dictionaries).
    @param new (non-nested obj): Object replacement (of type string, int, float, etc) - no lists or dictionaries).

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
                config[k] = config_loop_replace(v,old,new)
            elif isinstance(v,str) and isinstance(old,str) and isinstance(new,str):
                if old in v:
                    config[k] = v.replace(old,new)
            else:
                if old == v:
                    config[k] = new
    elif isinstance(config,list):
        for i,v in enumerate(config):
            if isinstance(v,dict) or isinstance(v,list):
                config[i] = config_loop_replace(v,old,new)
            elif isinstance(v,str):
                if old in v:
                    config[i] = v.replace(old,new)
            else:
                if old == v:
                    config[i] = new
    else:
        raise TypeError('Type "%s" not accepted in config_loop_replace.')

    return config
