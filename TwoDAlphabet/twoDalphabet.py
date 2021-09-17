import argparse, os, itertools, pandas, pickle
from collections import OrderedDict
from TwoDAlphabet.config import Config
from TwoDAlphabet.helpers import execute_cmd, get_config_dirs, parse_arg_dict, unpack_to_line, make_RDH
from TwoDAlphabet.alphawrap import Generic2D
import ROOT

class TwoDAlphabet:
    '''Class to injest and organize inputs.
    '''
    # mkdirs + rebin + organize_hists
    # track naming and fit info internally ("region information")
    def __init__(self,tag,jsons=[],findreplace={},externalOpts={}):
        '''Construct TwoDAlphabet object which takes as input an identification tag used to name
        the directory with stored results, JSON configuration files, find-replace pairs,
        and optional external arguments.

        Args:
            tag (str): Tag used to identify the outputs and create a project directory of the same name.
            jsons (list(str), optional): List of JSON configuration file paths. Defaults to [] in which case
                the tag is assumed to already exist and the existing runConfig.json is grabbed.
            findreplace (dict, optional):  Non-nested dictionary with key-value pairs to find-replace
                in the internal class configuration dict. Defaults to {}.
            externalOpts (dict, optional): Option-value pairs. Defaults to {}.
        '''
        self.tag = tag
        self.config = None
        self.options = self.GetOptions(externalOpts)
        if jsons == [] and os.path.isdir(self.tag+'/'):
            print ('Attempting to grab existing runConfig ...')
            jsons = [d+'/runConfig.json' for d in get_config_dirs(self.tag+'/')]
        if jsons == []:
            raise RuntimeError('No jsons were input and no existing ones could be found.')
        for j in jsons:
            self.AddConfig(j,findreplace)

        self.config.Construct()
        self.nbkgs, self.nsignals = self.config.nbkgs, self.config.nsignals
        self.SetupProjDir()
        if self.options.draw == False:
            ROOT.gROOT.SetBatch(True)

        self.alphaObjs = pandas.DataFrame(columns=['process','region','nuisances','obj'],index=['process','region'])
        self.alphaParams = pandas.DataFrame(columns=['name','obj','constraint'],index='name') # "name":"constraint"

    def GetOptions(self,externalOpts):
        '''General arguments passed to the project. Options specified in 
        Config.GetOptions() can also be provided to set an option globally
        to all Config objects being tracked.

        @param externalOpts (dict): Option-value pairs.

        Returns:
            ArgumentParser.Namespace
        '''
        parser = argparse.ArgumentParser()
        parser.add_argument('verbosity', default=0, type=int, nargs='?',
            help="Save settings to file in json format. Ignored in json file")
        parser.add_argument('overwrite', default=False, type=bool, nargs='?',
            help="Delete project directory if it exists. Defaults to False.")
        parser.add_argument('debugDraw', default=False, type=bool, nargs='?',
            help="Draw all canvases while running for the sake of debugging. Useful for developers only. Defaults to False.")
        return parse_arg_dict(parser,externalOpts)

    def AddConfig(self,jsonFileName,findreplace,onlyOn=['process','region']):
        '''Add a json configuration to process and track.

        Args:
            jsonFileName (str): JSON configuration file name/path.
            findreplace (dict): Non-nested dictionary with key-value pairs to find-replace
                in the internal class configuration dict.
            onlyOn (list(str),str): Column name or list of columns to match for determining duplicates.
                Options are either 'process' or 'region'. Defaults is ['process','region'] and both are considered.
        '''
        inputConfig = Config(jsonFileName,findreplace,externalOptions=vars(self.options))
        if self.configs == None:
            self.config == inputConfig
        else:
            self.config.Add(inputConfig,onlyOn)

    def AddAlphaObj(self,process,region,obj,ptype='BKG'):
        '''Start

        Args:
            process ([type]): [description]
            region ([type]): [description]
            obj ([type]): [description]
            ptype ([str]): 'BKG' or 'SIGNAL'.
        '''
        if not isinstance(Generic2D):
            raise RuntimeError('Can only tack objects of type Generic2D.')
        if ptype not in ['BKG','SIGNAL']:
            raise RuntimeError('Process type (ptype) can only be BKG or SIGNAL.')
        self._checkAgainstConfig(process,region)

        for cat in ['LOW','SIG','HIGH']:
            rph = obj.RooParametricHist()
            model_obj_row = {
                "process": process,
                "region": region+'_'+cat,
                "obj": rph[cat],
                "process_type": ptype,
                "combine_idx": self.nbkgs+1 if ptype == 'BKG' else -1*self.nsignals
            }
            self.alphaObjs.append(model_obj_row,verify_integrity=True)
        
        if ptype == 'BKG': self.nbkgs+=1
        elif ptype == 'SIGNAL': self.nsignals+=1

        nuis_obj_cols = ['name','obj','constraint']
        for n in obj.nuisances:
            self.alphaParams.append({c:n[c] for c in nuis_obj_cols},verify_integrity=True)

    def _checkAgainstConfig(self,process,region):
        all_pairs = [g[0] for g in self.config.df.groupby(['p','r'])]
        if (process,region) not in all_pairs:
            raise RuntimeError('Attempting to track an object for process "%s" and region "%s" but that pair does not exist among those defined in the config:\n\t%s'%(process,region,all_pairs))

    def SetupProjDir(self):
        '''Create the directory structure where results will be stored.
        '''
        if not os.path.isdir(self.tag+'/'):
            if self.options.overwrite: 
                execute_cmd('rm -rf '+self.tag)
            print ('Making dir '+self.tag+'/')
            os.mkdir(self.tag+'/')

        for c in self.configs.values():
            dirs_to_make = [
                self.tag+'/'+c.name,
                self.tag+'/'+c.name+'/plots/',
                self.tag+'/'+c.name+'/plots/fit_b/',
                self.tag+'/'+c.name+'/plots/fit_s/',
            ]
            if c.options.plotUncerts and not os.path.isdir(self.tag+'/'+c.name+'/UncertPlots/'): 
                dirs_to_make.append(self.tag+'/'+c.name+'/UncertPlots/')

            for d in dirs_to_make:
                if not os.path.isdir(d):
                    os.mkdir(d)

    def _saveOut(self):
        '''Save individual configs to project directory.
        '''
        self.config.SaveOut()
        pickle.dump(self,self.projDir+'twoDobj_%s.p'%self.tag)

    def _makeCard(self):
        card_new = open(self.config.projPath+'card_'+self.tag+'.txt','w')
        # imax (bins), jmax (backgrounds+signals), kmax (systematics) 
        imax = str(3*self.GetNregions()) # pass, fail for each 'X' axis category           
        jmax = self.nbkgs + self.nsignals
        kmax = self.config.GetNsystematics() # does not include alphaParams
        channels = ['_'.join(r)+'_'+self.name for r in itertools.product(self.df.regions.unique(),['LOW','SIG','HIGH'])]
        
        card_new.write('imax '+imax+'\n')      
        card_new.write('jmax '+jmax+'\n')
        card_new.write('kmax '+kmax+'\n')
        card_new.write('-'*120+'\n')

        # Shapes
        shape_line = 'shapes  {0:20} * base_{1}.root w_{1}:{0}_$CHANNEL w_{1}:{0}_$CHANNEL_$SYSTEMATIC\n'
        for proc in self.config.df.process.unique():
            if proc == 'data_obs': continue
            card_new.write(shape_line.format(proc,self.name))

        shape_line_nosyst = shape_line.replace('_$SYSTEMATIC','')
        for proc in ['data_obs']+self.alphaObjs.process.unique():
            card_new.write(shape_line_nosyst.format(proc,self.name))

        card_new.write('-'*120+'\n')

        # Set bin observation values to -1
        card_new.write('bin %s'%(unpack_to_line(channels)))
        card_new.write('observation %s'%unpack_to_line([-1 for i in range(imax)]))

        card_new.write('-'*120+'\n')

        ######################################################
        # Tie processes to bins and rates and simultaneously #
        # create the systematic uncertainty rows             #
        ######################################################
        bin_line         = '{0:20} '.format('bin')
        processName_line = '{0:20} '.format('process')
        processCode_line = '{0:20} '.format('process')
        rate_line        = '{0:20} '.format('rate')
        syst_lines = OrderedDict()

        # Fill syst_lines with keys to initialized strings
        for syst,syst_group in self.config.df.groupby(by='variation',sort=True):
            if syst == 'nominal': continue
            syst_type = syst_group.loc[0,'syst_type']
            syst_lines[syst] = '{0:20} {1:20} '%(syst, syst_type)

        # Work with template bkgs first
        for pair,group in self.config.df.groupby(['process','region']):
            proc,region = pair
            combine_idx = group.combine_idx[0]
            for cat in ['LOW','SIG','HIGH']:
                chan = '%s_%s_%s'%(region,cat,self.name)

                bin_line += '{0:20} '.format(chan)
                processName_line += '{0:20} '.format(proc)
                processCode_line += '{0:20} '.format(combine_idx)
                rate_line += '{0:20} '.format('1')

                for syst in syst_lines.keys():
                    if syst in group.systematics.unique():
                        syst_effect = group.apply(lambda row: row[row.syst_type])[0]
                    else:
                        syst_effect = '-'

                    syst_lines[syst] += '{0:20} '.format(syst_effect)

        # Now work with alpha objects
        # NOTE: duplicated code but no good way to combine without making things confusing
        for pair,group in self.alphaObjs.groupby(['process','region']):
            proc,region = pair
            combine_idx = group.combine_idx[0]
            for cat in ['LOW','SIG','HIGH']:
                chan = '%s_%s_%s'%(region,cat,self.name)

                bin_line += '{0:20} '.format(chan)
                processName_line += '{0:20} '.format(proc)
                processCode_line += '{0:20} '.format(combine_idx)
                rate_line += '{0:20} '.format('-1')

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
        for param in self.alphaParams.itertuples():
            card_new.write('{0:40} {1}'.format(param.name,param.constraint))
        
        card_new.close() 

    def _makeWorkspace(self):
        var_lists = {c:ROOT.RooArgList(self.binning.xVars[c],self.binning.yVar) for c in ['LOW','SIG','HIGH']}

        print ("Making workspace...")
        workspace = ROOT.RooWorkspace("w_"+self.name)
        for hname in self.config.organizedHists.GetHistNames():
            if hname.endswith('_FULL'):
                continue
            
            print ('Making RooDataHist... %s'%hname)
            cat = hname.split('_')[-1]
            rdh = make_RDH(self.config.organizedHists.Get(hname), var_lists[cat])
            getattr(workspace,'import')(rdh,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())

        for rph in self.alphaObjs.obj:
            print ('Adding RooParametricHist... %s'%rph.GetName())
            getattr(workspace,'import')(rph,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())

    def Construct(self):
        self._makeCard()
        self._makeWorkspace()
        self._saveOut()

    def plot(self):
        # plotter.methodA()
        # plotter.methodB()
        # ...
        pass

    def MLfit(self):
        pass

    def GoodnessOfFit(self):
        pass

    def SignalInjection(self,r):
        pass

    def Limits(self):
        pass

    def Impacts(self):
        pass

def SetSnapshot(d=''):
    w_f = ROOT.TFile.Open(d+'higgsCombineTest.FitDiagnostics.mH120.root')
    w = w_f.Get('w')
    fr_f = ROOT.TFile.Open(d+'fitDiagnosticsTest.root')
    fr = fr_f.Get('fit_b')
    myargs = ROOT.RooArgSet(fr.floatParsFinal())
    w.saveSnapshot('initialFit',myargs,True)
    fout = ROOT.TFile('initialFitWorkspace.root',"recreate")
    fout.WriteTObject(w,'w')
    fout.Close()