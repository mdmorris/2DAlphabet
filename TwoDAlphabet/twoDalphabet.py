import argparse, os, itertools, pandas, glob, pickle, sys, re, json
from collections import OrderedDict
from TwoDAlphabet.config import Config, OrganizedHists
from TwoDAlphabet.binning import Binning
from TwoDAlphabet.helpers import CondorRunner, execute_cmd, parse_arg_dict, unpack_to_line, make_RDH, cd
from TwoDAlphabet.alphawrap import Generic2D
from TwoDAlphabet import plot
import ROOT

class TwoDAlphabet:
    '''Class to injest and organize inputs.
    '''
    # mkdirs + rebin + organize_hists
    # track naming and fit info internally ("region information")
    def __init__(self, tag, inJSON='', findreplace={}, externalOpts={}, loadPrevious=False):
        '''Construct TwoDAlphabet object which takes as input an identification tag used to name
        the directory with stored results, JSON configuration files, find-replace pairs,
        and optional external arguments.

        Args:
            tag (str): Tag used to identify the outputs and create a project directory of the same name.
            jsons (str, optional): JSON configuration file path. Defaults to '' in which case
                the tag is assumed to already exist and the existing runConfig.json is grabbed.
            findreplace (dict, optional):  Non-nested dictionary with key-value pairs to find-replace
                in the internal class configuration dict. Defaults to {}.
            externalOpts (dict, optional): Option-value pairs. Defaults to {}.
        '''
        self.tag = tag
        if inJSON == '' and os.path.isdir(self.tag+'/'):
            print ('Attempting to grab existing runConfig ...')
            inJSON = glob.glob(self.tag+'/runConfig.json')
        if inJSON == '':
            raise RuntimeError('No JSONs were input and no existing ones could be found.')

        config = Config(inJSON, findreplace)
        optdict = config._section('OPTIONS')
        optdict.update(externalOpts)
        self.options = self.LoadOptions(optdict)
        self.df = config.FullTable()
        self.subtagTracker = {}
        self.iterWorkspaceObjs = config.iterWorkspaceObjs
        self._binningMap = {r:config._section('REGIONS')[r]['BINNING'] for r in config._section('REGIONS').keys()}
        self.ledger = Ledger(self.df)

        if not loadPrevious:
            self._setupProjDir()
            template_file = ROOT.TFile.Open(self.df.iloc[0].source_filename)
            template = template_file.Get(self.df.iloc[0].source_histname)
            template.SetDirectory(0)

            self.binnings = {}
            for kbinning in config._section('BINNING').keys():
                self.binnings[kbinning] = Binning(kbinning, config._section('BINNING')[kbinning], template)

            self.organizedHists = OrganizedHists(
                self.tag+'/', self.binnings,
                self.GetHistMap(), readOnly=False
            )

        else:
            self.binnings = pickle.load(open(self.tag+'/binnings.p','rb'))
            self.organizedHists = OrganizedHists(
                self.tag+'/', self.binnings,
                self.GetHistMap(), readOnly=True
            )
            # Does not contain the RooFit objects - just meta info
            self.ledger.alphaObjs = pandas.read_csv(self.tag+'/alphaObjs.csv')
            self.ledger.alphaParams = pandas.read_csv(self.tag+'/alphaParams.csv')
        
        if self.options.debugDraw is False:
            ROOT.gROOT.SetBatch(True)

        config.SaveOut(self.tag+'/')

    def _setupProjDir(self):
        '''Create the directory structure where results will be stored.
        '''
        if not os.path.isdir(self.tag+'/'):
            if self.options.overwrite: 
                execute_cmd('rm -rf '+self.tag)
            print ('Making dir '+self.tag+'/')
            os.mkdir(self.tag+'/')

        if self.options.plotTemplateComparisons and not os.path.isdir(self.tag+'/UncertPlots/'): 
            os.mkdir(self.tag+'/UncertPlots/')

    def LoadOptions(self, nonDefaultOpts={}):
        '''Optional arguments passed to the project.
        Options can be specified in the JSON config file (from 'OPTIONS' section)
        or via the externalOpts argument dictionary. The options provided by externalOpts
        override those from the config and modify the config in-place so that the
        version later saved reflects the conditions under which the config was used.

        @param externalOpts (dict): Option-value pairs.

        Returns:
            ArgumentParser.Namespace
        '''
        parser = argparse.ArgumentParser()
        # General
        parser.add_argument('verbosity', default=0, type=int, nargs='?',
            help="Save settings to file in JSON format. Ignored in JSON file")
        parser.add_argument('overwrite', default=False, type=bool, nargs='?',
            help="Delete project directory if it exists. Defaults to False.")
        parser.add_argument('debugDraw', default=False, type=bool, nargs='?',
            help="Draw all canvases while running for the sake of debugging. Useful for developers only. Defaults to False.")
        # Blinding
        parser.add_argument('blindedPlots', default=[], type=str, nargs='*',
            help='List of regions in which to blind plots of x-axis SIG. Does not blind fit.')
        parser.add_argument('blindedFit', default=[], type=str, nargs='*',
            help='List of regions in which to blind fit of x-axis SIG. Does not blind plots.')
        # Plotting
        parser.add_argument('haddSignals', default=True, type=bool, nargs='?',
            help='Combine signals into one histogram for the sake of plotting. Still treated as separate in fit. Defaults to True.')
        parser.add_argument('plotTitles', default=False, type=bool, nargs='?',
            help='Include titles in plots. Defaults to False.')
        parser.add_argument('plotTemplateComparisons', default=False, type=bool, nargs='?',
            help='Plot comparison of pre-fit uncertainty shape templates in 1D projections. Defaults to False.')
        parser.add_argument('plotPrefitSigInFitB', default=False, type=bool, nargs='?',
            help='In the b-only post-fit plots, plot the signal normalized to its pre-fit value. Defaults to False.')
        parser.add_argument('plotEvtsPerUnit', default=False, type=bool, nargs='?',
            help='Post-fit bins are plotted as events per unit rather than events per bin. Defaults to False.')
        parser.add_argument('year', default=1, type=int, nargs='?',
            help='Year information used for the sake of plotting text. Defaults to 1 which indicates that the full Run 2 is being analyzed.')

        if nonDefaultOpts != {}:
            out = parse_arg_dict(parser,nonDefaultOpts)
        else:
            out = parser.parse_args([])
        return out

    def _saveOut(self):
        '''Save to project directory:
        - the full model table in csv (and markdown if py3)
        - the binnings dictionary (with objects)
        - the alphaObjs and alphaParams dictionaries (without objects)
        '''
        if 'index' in self.df.columns:
               df = self.df.reset_index(drop=True).drop('index',axis=1)
        else:  df = self.df

        df.to_csv(self.tag+'/hist_ledger.csv')
        if sys.version_info.major == 3:
            df.to_markdown(self.tag+'/hist_ledger.md')

        pickle.dump(self.binnings,open(self.tag+'/binnings.p','wb'))

        alphaObjs = self.ledger.alphaObjs.drop(columns=['obj','norm'])
        alphaObjs.to_csv(self.tag+'/alphaObjs.csv')

        alphaParams = self.ledger.alphaParams.drop(columns=['obj'])
        alphaParams.to_csv(self.tag+'/alphaParams.csv')

        if self.options.plotTemplateComparisons:
            plot.make_systematic_plots(self)

# --------------AlphaObj INTERFACE ------ #
    def AddAlphaObj(self, process, region, obj, ptype='BKG', color=ROOT.kYellow):
        '''Start

        Args:
            process ([type]): [description]
            region ([type]): [description]
            obj ([type]): [description]
            ptype ([str]): 'BKG' or 'SIGNAL'.
        '''
        if not isinstance(obj, Generic2D):
            raise RuntimeError('Can only tack objects of type Generic2D.')
        if ptype not in ['BKG','SIGNAL']:
            raise RuntimeError('Process type (ptype) can only be BKG or SIGNAL.')
        self.ledger._checkAgainstConfig(process, region)

        # for cat in ['LOW','SIG','HIGH']:
        rph,norm = obj.RooParametricHist()
        model_obj_row = {
            "process": process,
            "region": region,
            "obj": rph,
            "norm": norm,
            "process_type": ptype,
            "color": color,
            'title': process
        }

        self.ledger.alphaObjs = self.ledger.alphaObjs.append(model_obj_row, ignore_index=True)

        nuis_obj_cols = ['name', 'obj', 'constraint']
        for n in obj.nuisances:
            self.ledger.alphaParams = self.ledger.alphaParams.append({c:n[c] for c in nuis_obj_cols}, ignore_index=True)

# --------------- GETTERS --------------- #
    def InitQCDHists(self):
        '''Loop over all regions and for a given region's data histogram, subtract the list of background histograms,
        and return data-bkgList.

        Returns:
            dict(region,TH2): Dictionary with regions as keys and values as histograms of data-bkgList.
        '''
        out = {}
        for region,group in self.df.groupby('region'):
            data_hist = self.organizedHists.Get(process='data_obs',region=region,systematic='')
            qcd = data_hist.Clone(data_hist.GetName().replace('data_obs','qcd'))
            qcd.SetDirectory(0)

            bkg_sources = group.loc[group.process_type.eq('BKG') & group.variation.eq('nominal')]['process']
            for process_name in bkg_sources.to_list():
                bkg_hist = self.organizedHists.Get(process=process_name,region=region,systematic='')
                qcd.Add(bkg_hist,-1)

            out[region] = qcd
            
        return out

    def GetHistMap(self, df=None):
        '''Collect information on the histograms to extract, manipulate, and save
        into organized_hists.root and store it inside a `dict` where the key is the
        filename and the value is a DataFrame with columns `source_histname`, `out_histname`,
        `scale`, `color`, and `binning`. Only accounts for "FULL" category and does not 
        contain information on subspaces.

        Args:
            df (pandas.DataFrame, optional): pandas.DataFrame to groupby. Defaults to None in which case self.df is used.

        Returns:
            dict(str:pandas.DataFrame): Map of file name to DataFrame with information on histogram name, scale factor, fill color, and output name.
        '''
        def _get_out_name(row):
            '''Per-row processing to create the output histogram name.

            Args:
                row (pandas.Series): Row to process.

            Returns:
                str: New name.
            '''
            if row.variation == 'nominal':
                return row.process+'_'+row.region+'_FULL'
            else:
                return row.process+'_'+row.region+'_FULL_'+row.variation+row.direction
        
        if not isinstance(df,pandas.DataFrame):
            df = self.df

        hists = {}
        for g, group_df in df.groupby(['source_filename']):
            out_df = group_df.copy(True)
            out_df = out_df[out_df['variation'].eq('nominal') | out_df["syst_type"].eq("shapes")]
            out_df['out_histname'] = out_df.apply(_get_out_name, axis=1)
            out_df['binning'] = out_df.apply(lambda row: self._binningMap[row.region], axis=1)
            hists[g] = out_df[['source_histname','out_histname','scale','color','binning']]
        return hists

    def GetBinningFor(self, region):
        for r,b in self._binningMap.items():
            if r == region:
                return self.binnings[b], b
        
        raise RuntimeError('Cannot find region (%s) in config:\n\t%s'%(region,self._binningMap))

    def _getCatNameRobust(self, hname):
        if hname.split('_')[-1] in ['FULL','SIG','HIGH','LOW']: # simplest case
            out =  hname.split('_')[-1]
        else: # systematic variation so need to be careful
            this_rname = False
            for rname in self.ledger.GetRegions():
                if rname in hname:
                    this_rname = rname
                    break

            if not this_rname: raise RuntimeError('Could not find a region name from %s in "%s"'%(self.ledger.GetRegions(),hname))

            start_idx = hname.index(this_rname)+len(this_rname)+1 #+1 for the trailing underscore
            end_idx = hname.index('_',start_idx)
            out = hname[start_idx:end_idx]
        
        return out

# ---------- FIRST STEP CONSTRUCTION ------ #
    def _makeWorkspace(self):
        var_lists = {}
        for binningName in self.binnings.keys():
            var_lists[binningName] = {
                c:ROOT.RooArgList(self.binnings[binningName].xVars[c], self.binnings[binningName].yVar) for c in ['LOW','SIG','HIGH']
            }

        print ("Making workspace...")
        workspace = ROOT.RooWorkspace("w")
        for hname in self.organizedHists.GetHistNames():
            cat = self._getCatNameRobust(hname)
            if cat == 'FULL':
                continue
            
            binningName = self.organizedHists.BinningLookup(hname.replace(cat,'FULL'))

            print ('Making RooDataHist... %s'%hname)
            rdh = make_RDH(self.organizedHists.Get(hname), var_lists[binningName][cat])
            getattr(workspace,'import')(rdh)

        for rph in self.ledger.alphaObjs.obj:
            for rph_cat in rph.values():
                print ('Adding RooParametricHist... %s'%rph_cat.GetName())
                getattr(workspace,'import')(rph_cat,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        for norm in self.ledger.alphaObjs.norm:
            for norm_cat in norm.values():
                print ('Adding RooParametricHist norm... %s'%norm_cat.GetName())
                getattr(workspace,'import')(norm_cat,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())

        out = ROOT.TFile.Open('base.root','RECREATE')
        out.cd()
        workspace.Write()
        out.Close()

    def Construct(self):
        with cd(self.tag):
            self._makeWorkspace()
        self._saveOut()

    def MakeCard(self, subledger, subtag, workspaceDir='../', toyData=None):
        with cd(self.tag):
            _runDirSetup(subtag)
            MakeCard(subledger, subtag, workspaceDir, toyData)

# -------- STAT METHODS ------------------ #
    def MLfit(self, subtag, rMin=-1, rMax=10, setParams={}, verbosity=0, usePreviousFit=False, toyData=None):
        with cd(self.tag+'/'+subtag):
            _runMLfit(
                subtag='',
                blinding=self.options.blindedFit,
                verbosity=verbosity, 
                rMin=rMin, rMax=rMax,
                setParams=setParams,
                usePreviousFit=usePreviousFit)
            make_postfit_workspace('')
            # systematic_analyzer_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/systematicsAnalyzer.py card.txt --all -f html > systematics_table.html'
            # execute_cmd(systematic_analyzer_cmd)    

    def StdPlots(self, subtag, ledger=None):
        run_dir = self.tag+'/'+subtag
        with cd(run_dir):
            if ledger == None:
                ledger = pickle.load(open('ledger.p','rb'))
            plot.nuis_pulls()
            plot.save_post_fit_parametric_vals()
            plot.make_correlation_matrix( # Ignore nuisance parameters that are bins
                varsToIgnore=self.ledger.alphaParams.name.str.contains('_bin_\d+-\d+').to_list()
            )
            plot.gen_post_fit_shapes()
            plot.gen_projections(ledger, self, 'b')
            plot.gen_projections(ledger, self, 's')
            
    def GetParamsOnMatch(self, regex='', subtag='', b_or_s='b'):
        out = {}

        f = ROOT.TFile.Open(self.tag+'/'+subtag+'/fitDiagnosticsTest.root')
        fr = f.Get('fit_'+b_or_s)
        final_pars = ROOT.RooArgList(fr.floatParsFinal())

        for i in range(final_pars.getSize()):
            par = final_pars[i]
            if re.search(regex, par.GetName()):
                out[par.GetName()] = {'val': par.getValV(), 'error': par.getError()}

        f.Close()

        return out

    def GenerateToys(self, name, subtag='', card=None, workspace=None, ntoys=1, seed=123456, expectSignal=0, setParams={}, freezeParams=[], run=True):
        if card == None and workspace == None:
            raise IOError('Either card or workspace must be provided relative to the directory where the generation will be run.')
        elif card != None and workspace != None:
            raise IOError('Only one of card or workspace can be provided.')

        if card:
            if isinstance(card, str):
                card_name = card
            elif isinstance(card, bool) and card == True:
                card_name = 'card.txt'
            workspace_file = '%s_gen_workspace.root'%(name if subtag=='' else subtag+'/'+name)
            execute_cmd(
                'text2workspace.py -b {0}/{1} -o {0}/{2} --channel-masks --X-no-jmax'.format(self.tag,card_name, workspace_file)
            )
            
            input_opt = '-d %s'%workspace_file
        
        elif workspace:
            if isinstance(workspace, str):
                workspace_file = workspace.split(':')
                input_opt = '-d %s --snapshotName %s'%(workspace.split(':'))
            elif isinstance(workspace, bool) and workspace == True:
                workspace_file = 'initialFitWorkspace.root'
                input_opt = '-d %s --snapshotName initialFit'%workspace_file

        masks_off = ['%s=0'%mask for mask in self._getMasks(self.tag+'/'+workspace_file)]

        run_dir = self.tag+'/'+subtag
        with cd(run_dir):
            param_vals = ['%s=%s'%(k,v) for k,v in setParams.items()] + masks_off
            param_opt = ''
            if len(param_vals) > 0:
                param_opt = '--setParameters '+','.join(param_vals),

            freeze_opt = ''
            if len(freezeParams) > 0:
                freeze_opt = '--freezeParameters '+','.join(freezeParams)

            gen_command_pieces = [
                'combine -M GenerateOnly',
                 input_opt, param_opt,
                '--toysFrequentist --bypassFrequentistFit',
                '-t %s'%ntoys,
                '--saveToys -s %s'%seed,
                '--expectSignal %s'%expectSignal,
                '-n _%s'%name,
                freeze_opt
            ]
            if run:
                execute_cmd(' '.join(gen_command_pieces))
            else:
                print ('Not running:\n\t'+' '.join(gen_command_pieces))

        return '%shiggsCombine_%s.GenerateOnly.mH120.%s.root'%(subtag+'/' if subtag!='' else '',name,seed)

    def _getMasks(self, filename):
        masked_regions = []
        f = ROOT.TFile.Open(filename)
        w = f.Get('w')
        allVars = ROOT.RooArgList(w.allVars())

        for i in range(allVars.getSize()):
            var = allVars[i]
            if 'mask_' in var.GetName() and '_SIG_' in var.GetName():
                if var.getValV() == 1:
                    masked_regions.append(var.GetName())
        f.Close()

        return masked_regions

    def GoodnessOfFit(self):
        pass

    def SignalInjection(self, r):
        pass

    def Limit(self, subtag, loadFitDir, blindData=True, verbosity=0, setParams={}, location='', eosRootfiles=''):
        if subtag == '': 
            raise RuntimeError('The subtag for limits must be non-empty so that the limit will be run in a nested directory.')
        if location == '':
            print ('Automatically determining whether location is "local" or "condor".')
            if 'condor' in os.getcwd():
                location = 'condor'
            else:
                location = 'local'
        elif location not in ['local','condor']:
            raise RuntimeError('Limit location can only be "", "local", or "condor". If empty string, will determine if condor based on absolute path.')

        run_dir = self.tag+'/'+subtag
        _runDirSetup(run_dir)
        workspace_dir = '' if subtag=='' else '../'
        with cd(run_dir):
            subledger = pickle.load(open('ledger.p','rb'))
            self.MakeCard(subledger, subtag, workspace_dir)
            limit_cmd = _runLimit(loadFitDir, blindData, verbosity, setParams, location) # runs on this line if location == 'local'
            
            if location == 'condor':
                condor = CondorRunner(
                    name=self.tag+'_'+subtag,
                    primaryCmds=[limit_cmd],
                    toPkg=run_dir,
                    toGrab=run_dir+'/higgsCombineTest.AsymptoticLimits.mH120.root',
                    eosRootfileTarball=eosRootfiles
                )
                condor.submit()
                
    def Impacts(self):
        pass

class Ledger():
    def __init__(self, df):
        self.df = df
        self.alphaObjs = pandas.DataFrame(columns=['process','region','obj','norm','process_type','color','combine_idx','title'])
        self.alphaParams = pandas.DataFrame(columns=['name','obj','constraint'])

    def append(self, toAppend):
        self.df.append(toAppend, ignore_index=True if isinstance(toAppend, dict) else False)

    def select(self,f,*args):
        eval_lamba = lambda row: f(row,args)
        df = self.df.loc[self.df.apply(eval_lamba, axis=1)]
        new_ledger = Ledger(df)
        new_ledger.alphaObjs = self.alphaObjs.loc[self.alphaObjs.apply(eval_lamba, axis=1)]
        return new_ledger

    def GetRegions(self):
        return self.df.region.unique()

    def GetProcesses(self, ptype='', includeNonConfig=True, includeConfig=True):
        if ptype not in ['','SIGNAL','BKG','DATA']:
            raise ValueError('Process type "%s" not accepted. Must be empty string or one of "SIGNAL","BKG","DATA".'%ptype)

        proc_list = []
        if includeConfig:
            if ptype == '':
                to_add = self.df.process.unique()
            else:
                to_add = self.df[self.df.process_type.eq(ptype)].process.unique()
            
            proc_list.extend(list(to_add))

        if includeNonConfig and self.alphaObjs.process.unique().size > 0:
            if ptype == '':
                to_add = self.alphaObjs.process.unique()
            else:
                to_add = self.alphaObjs[self.alphaObjs.process_type.eq(ptype)].process.unique()

            proc_list.extend(list(to_add))
        return proc_list

    def GetProcRegPairs(self):
        return [g[0] for g in self.df.groupby(['process','region'])]+[g[0] for g in self.alphaObjs.groupby(['process','region'])]

    def GetShapeSystematics(self):
        return self.df.variation.unique()

    def GetAlphaSystematics(self):
        return self.alphaParams.name.unique()

    def GetAllSystematics(self):
        return self.GetShapeSystematics()+self.GetAlphaSystematics()

    def _getProcessAttrBase(self, procName, attrName):
        if procName in self.df.process.unique():
            df = self.df
        elif procName in self.alphaObjs.process.unique():
            df = self.alphaObjs
        else:
            raise NameError('Process "%s" does not exist.'%procName)
        return get_process_attr(df, procName, attrName)

    def GetProcessColor(self, procName):
        return self._getProcessAttrBase(procName,'color')

    def GetProcessType(self, procName):
        return self._getProcessAttrBase(procName,'process_type')

    def GetProcessTitle(self, procName):
        return self._getProcessAttrBase(procName,'title')

    def IsSignal(self, procName):
        return self.GetProcessType(procName) == 'SIGNAL'

    def IsBackground(self, procName):
        return self.GetProcessType(procName) == 'BKG'

    def IsData(self, procName):
        return self.GetProcessType(procName) == 'DATA'

    @property
    def nsignals(self):
        return self.df[self.df.process_type.eq('SIGNAL')].process.nunique() + self.alphaObjs[self.alphaObjs.process_type.eq('SIGNAL')].process.nunique()

    @property
    def nbkgs(self):
        return self.df[self.df.process_type.eq('BKG')].process.nunique() + self.alphaObjs[self.alphaObjs.process_type.eq('BKG')].process.nunique()

    def _checkAgainstConfig(self, process, region):
        if (process,region) in self.GetProcRegPairs():
            raise RuntimeError('Attempting to track an object for process-region pair (%s,%s) that already exists among those defined in the config:\n\t%s'%(process,region,self.GetProcesses()))
        if region not in self.GetRegions():
            raise RuntimeError('Attempting to track an object for region "%s" but that region does not exist among those defined in the config:\n\t%s'%(region,self.GetRegions()))

    def _getCombineIdxMap(self):
        all_signals = self.df[self.df.process_type.eq('SIGNAL')].process.unique().tolist() + self.alphaObjs[self.alphaObjs.process_type.eq('SIGNAL')].process.unique().tolist()
        all_bkgs    = self.df[self.df.process_type.eq('BKG')].process.unique().tolist()    + self.alphaObjs[self.alphaObjs.process_type.eq('BKG')].process.unique().tolist()

        signal_map = pandas.DataFrame({'process': all_signals, 'combine_idx': [-1*i for i in range(0,len(all_signals))] })
        bkg_map    = pandas.DataFrame({'process': all_bkgs,    'combine_idx': [i for i in range(1,len(all_bkgs)+1)] })

        out = pandas.concat([signal_map, bkg_map])
        return out

def _runDirSetup(runDir):
    dirs_to_make = [
        runDir+'/',
        runDir+'/plots_fit_b/',
        runDir+'/plots_fit_s/',
        runDir+'/plots_fit_b/base_figs/',
        runDir+'/plots_fit_s/base_figs/',
    ]
    for d in dirs_to_make:
        if not os.path.isdir(d):
            os.mkdir(d)
    
    return runDir

def MakeCard(ledger, subtag, workspaceDir, toyData=None):
    combine_idx_map = ledger._getCombineIdxMap()

    card_new = open('%s/card.txt'%subtag,'w')
    # imax (bins), jmax (backgrounds+signals), kmax (systematics) 
    imax = 3*len(ledger.GetRegions()) # pass, fail for each 'X' axis category    
    jmax = ledger.nbkgs + ledger.nsignals -1
    kmax = len(ledger.GetShapeSystematics())-1 # -1 for nominal, does not include alphaParams
    channels = ['_'.join(r) for r in itertools.product(ledger.GetRegions(),['LOW','SIG','HIGH'])]
    
    card_new.write('imax %s\n'%imax)      
    card_new.write('jmax %s\n'%jmax)
    card_new.write('kmax %s\n'%kmax)
    card_new.write('-'*120+'\n')

    # Shapes
    shape_line = 'shapes  {0:20} * {1} w:{0}_$CHANNEL w:{0}_$CHANNEL_$SYSTEMATIC\n'
    for proc in ledger.GetProcesses(includeNonConfig=False):
        if proc == 'data_obs': continue
        card_new.write(shape_line.format(proc, workspaceDir+'base.root'))

    shape_line_nosyst = shape_line.replace(' w:{0}_$CHANNEL_$SYSTEMATIC','')
    for proc in list(ledger.alphaObjs.process.unique()):
        card_new.write(shape_line_nosyst.format(proc, workspaceDir+'base.root'))

    # Do data - importing toy if needed
    if toyData == None:
        card_new.write(shape_line_nosyst.format('data_obs', workspaceDir+'base.root'))
    else:
        card_new.write(shape_line_nosyst.replace(' w:{0}_$CHANNEL','').format('data_obs', toyData))

    card_new.write('-'*120+'\n')

    # Set bin observation values to -1
    card_new.write('bin                 %s\n'%(unpack_to_line(channels)))
    card_new.write('observation %s\n'%unpack_to_line([-1 for i in range(imax)]))

    card_new.write('-'*120+'\n')

    ######################################################
    # Tie processes to bins and rates and simultaneously #
    # create the systematic uncertainty rows             #
    ######################################################
    bin_line         = '{0:20} {1:20}'.format('bin','')
    processName_line = '{0:20} {1:20}'.format('process','')
    processCode_line = '{0:20} {1:20}'.format('process','')
    rate_line        = '{0:20} {1:20}'.format('rate','')
    syst_lines = OrderedDict()

    # Fill syst_lines with keys to initialized strings
    for syst,syst_group in ledger.df.groupby(by='variation',sort=True):
        if syst == 'nominal': continue
        syst_type = syst_group.iloc[0].syst_type
        syst_lines[syst] = '{0:20} {1:20} '.format(syst, syst_type)

    # Work with template bkgs first
    for pair, group in ledger.df.groupby(['process','region']):
        proc, region = pair
        if proc == 'data_obs': continue
        combine_idx = combine_idx_map[combine_idx_map.process.eq(proc)].combine_idx.iloc[0]
        for cat in ['LOW','SIG','HIGH']:
            chan = '%s_%s'%(region,cat)

            bin_line += '{0:20} '.format(chan)
            processName_line += '{0:20} '.format(proc)
            processCode_line += '{0:20} '.format(combine_idx)
            rate_line += '{0:20} '.format('-1')

            for syst in syst_lines.keys():
                if syst in group.variation.unique():
                    syst_effect = group.loc[group.variation.eq(syst)].apply(lambda row: row[row.syst_type],axis=1).iloc[0]
                else:
                    syst_effect = '-'

                syst_lines[syst] += '{0:20} '.format(syst_effect)

    # Now work with alpha objects
    # NOTE: duplicated code but no good way to combine without making things confusing
    for pair,group in ledger.alphaObjs.groupby(['process', 'region']):
        proc,region = pair
        combine_idx = combine_idx_map[combine_idx_map.process.eq(proc)].combine_idx.iloc[0]
        for cat in ['LOW','SIG','HIGH']:
            chan = '%s_%s'%(region, cat)

            bin_line += '{0:20} '.format(chan)
            processName_line += '{0:20} '.format(proc)
            processCode_line += '{0:20} '.format(combine_idx)
            rate_line += '{0:20} '.format('1')

            for syst in syst_lines.keys():
                syst_lines[syst] += '{0:20} '.format('-')

    card_new.write(bin_line+'\n')
    card_new.write(processName_line+'\n')
    card_new.write(processCode_line+'\n')
    card_new.write(rate_line+'\n')
    card_new.write('-'*120+'\n')
    for line_key in syst_lines.keys():
        card_new.write(syst_lines[line_key]+'\n')

    ######################################################
    # Mark floating values as flatParams                 #
    # We float just the rpf params and the failing bins. #
    ######################################################
    for param in ledger.alphaParams.itertuples():
        card_new.write('{0:40} {1}\n'.format(param.name, param.constraint))
    
    card_new.close()
    pickle.dump(ledger, open('%s/ledger.p'%subtag,'wb'))

def _runMLfit(subtag, blinding, verbosity, rMin, rMax, setParams, usePreviousFit=False):
    if usePreviousFit: param_options = ''
    else:              param_options = '--text2workspace "--channel-masks" '
    params_to_set = ','.join(['mask_%s_SIG=1'%r for r in blinding]+['%s=%s'%(p,v) for p,v in setParams.items()]+['r=1'])
    param_options += '--setParameters '+params_to_set

    fit_cmd = 'combine -M FitDiagnostics {card_or_w} {param_options} --saveWorkspace --cminDefaultMinimizerStrategy 0 --rMin {rmin} --rMax {rmax} -v {verbosity}'
    fit_cmd = fit_cmd.format(
        card_or_w='initifalFitWorkspace.root --snapshotName initialFit' if usePreviousFit else 'card%s.txt'%subtag,
        param_options=param_options,
        rmin=rMin,
        rmax=rMax,
        verbosity=verbosity
    )

    with open('FitDiagnostics_command.txt','w') as out:
        out.write(fit_cmd)

    if os.path.isfile('fitDiagnosticsTest.root'):
        execute_cmd('rm fitDiagnosticsTest.root')

    execute_cmd(fit_cmd)

def _runLimit(loadFitDir, blindData, verbosity, setParams, location, condor_name=''):
    fr_path = '../'+loadFitDir+'/fitDiagnosticsTest.root'
    if not os.path.exists(fr_path):
        raise OSError('The file %s does not exist. Please check that MLfit() finished correctly.'%fr_path)

    import_fitresult('card.txt', fr_path)

    param_options = ''
    if len(setParams) > 0:
        param_options = '--setParameters '+','.join('%s=%s'%(k,v) for k,v in setParams.items())

    limit_cmd = 'combine -M AsymptoticLimits -d morphedWorkspace.root --snapshotName morphedModel --saveWorkspace --cminDefaultMinimizerStrategy 0 {param_opt} {blind_opt} -v {verb}' 
    limit_cmd = limit_cmd.format(
        blind_opt='--run=blind' if blindData else '',
        param_opt=param_options,
        verb=verbosity
    )

    # Run combine if not on condor
    if location == 'local':   
        with open('Limit_command.txt','w') as out:
            out.write(limit_cmd) 
        execute_cmd(limit_cmd)

    return limit_cmd

def get_process_attr(df, procName, attrName):
    return df.loc[df.process.eq(procName)][attrName].iloc[0]

def make_postfit_workspace(d=''):
    w_f = ROOT.TFile.Open(d+'higgsCombineTest.FitDiagnostics.mH120.root')
    w = w_f.Get('w')
    fr_f = ROOT.TFile.Open(d+'fitDiagnosticsTest.root')
    fr = fr_f.Get('fit_b')
    myargs = ROOT.RooArgSet(fr.floatParsFinal())
    w.saveSnapshot('initialFit',myargs,True)
    fout = ROOT.TFile('initialFitWorkspace.root', "recreate")
    fout.WriteTObject(w, 'w')
    fout.Close()

# TODO: Add ability to freeze parameters via var.setConstant() while looping over floatParsFinal()
def import_fitresult(inCard, fitResult, toDrop=[]):
    # First convert the card and open workspace
    execute_cmd('text2workspace.py -b %s -o morphedWorkspace.root --channel-masks --X-no-jmax'%inCard)
    w_f = ROOT.TFile.Open('morphedWorkspace.root', 'UPDATE')
    w = w_f.Get('w')
    # Open fit result we want to import
    print ('Importing %s...'%fitResult)
    fr_f = ROOT.TFile.Open(fitResult)
    fr = fr_f.Get('fit_b') # b-only fit result (fit_s for s+b)
    myargs = fr.floatParsFinal()
    outargs = ROOT.RooArgSet()

    _cache = [] # to avoid seg faults
    for i in range(myargs.getSize()):
        var = myargs.at(i)
        if var.GetName() in toDrop: continue

        if not w.allVars().contains(var):
            print ('WARNING: Could not find %s'%var.GetName())
            continue
                
        outargs.add(var)
        _cache.append(var)
    
    w.saveSnapshot('morphedModel',outargs,True)
    w_f.WriteTObject(w,'w',"Overwrite")
    w_f.Close()