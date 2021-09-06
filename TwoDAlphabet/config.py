import ROOT, argparse, json, pickle, os, pandas, re
from numpy import nan
import pprint
pp = pprint.PrettyPrinter(indent=4)
from TwoDAlphabet.helpers import copy_update_dict, is_filled_list, nested_dict, open_json, parse_arg_dict, replace_multi
from TwoDAlphabet.binning import Binning, copy_hist_with_new_bins, get_bins_from_hist

_protected_keys = ["PROCESSES","SYSTEMATICS","REGIONS","BINNING","OPTIONS","GLOBAL","SCALE","COLOR","TYPE","X","Y","NAME","TITLE","BINS","NBINS","LOW","HIGH"]
_syst_col_defaults = {
    # 'variation': nan,
    'lnN': nan,
    'lnN_asym': nan,
    'shape_sigma': nan,
    'syst_type': nan,
    'source_filename': nan,
    'source_histname': nan,
    'direction': nan
}
class Config:
    '''Class to handle the reading and manipulation of data provided 
    in 2DAlphabet JSON configuration files. Constructor initializes
    a Config object for a given set of JSON files and performs
    all initial checks and manipulations.

    Args:
        json (str): File name and path.
        projPath (str): Project path.
        findreplace (dict, optional): Find-replace pairs. Defaults to {}.
        externalOptions (dict, optional): Extra key-value pairs to add to the OPTIONS
            section of the JSON file. Defaults to {}.
        loadPrevious (bool, optional): Load the previous histograms made instead of remaking. Defaults to False.

    Attributes:
        config (dict): JSON config opened as a dict.
        name (str): Name, inherited from the input JSON config.
        projPath (str): The project path + self.name.
        options (Namespace): Final set of options to use.
        loadPrevious (bool): Equal to the constructor arg of the same name.
    '''
    def __init__(self,json,projPath,findreplace={},externalOptions={},loadPrevious=False):
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
        self.df = self.FullTable()
        template_file = ROOT.TFile.Open(self.df.iloc[0].source_filename)
        template = template_file.Get(self.df.iloc[0].source_histname)
        template.SetDirectory(0)
        self.binning = Binning(self.name, self.Section('BINNING'), template)
        self.processes = self.Section('PROCESSES')
        self.systematics = self.Section('SYSTEMATICS')
        if self.loadPrevious: 
            self.organized_hists = OrganizedHists(self,readOnly=True)
        else:
            self.organized_hists = OrganizedHists(self,readOnly=False)

        return self

    def Section(self,key):
        '''Derive the dictionary for a given section of the configuration file
        with HELP keys removed.

        Args:
            key (str): Section name (all capital letters) of the configuration file.

        Returns:
            dict: Section of config.
        '''
        if not isinstance(self.config[key],dict):
            raise TypeError('Section access is only for sub-dictionaries. Try accessing directly reading config[key].')
        return {k:v for k,v in self.config[key].items() if k != 'HELP'}

    def _addFindReplace(self,findreplace):
        '''Add external find-replace pairs to the "GLOBAL" entry of config.

        Args:
            findreplace (dict): Find-replace pairs in non-nested dictionary.

        Raises:
            ValueError: If a "find" already exists in self.config['GLOBAL'].
        '''
        for s in findreplace:
            if s in self.Section('GLOBAL'):
                raise ValueError('A command line string replacement (%s) conflicts with one already in the configuration file. Quitting...' %(s))
            self.Section('GLOBAL')[s] = findreplace[s]

    def _varReplacement(self):
        '''Do string substitution for config entries based on the dictionary of find
        and replace values found in self.config['GLOBAL']. Call self._addFindReplace()
        before running this function to add in external find-replace pairs.

        Args:
            findreplace (dict): Non-nested dictionary with key-value pairs to find-replace
                in the internal class configuration dict.

        Returns:
            None.
        '''
        print ("Doing GLOBAL variable replacement in input json...")
        for old_string in self.Section('GLOBAL'):
            if old_string == "HELP":
                print ('WARNING: The HELP entry is deprecated and checking for it will be removed in the future. Please remove it from your config.')
                continue
            new_string = self.Section('GLOBAL')[old_string]
            self.config = config_loop_replace(self.config, old_string, new_string)

    def GetOptions(self,externalOpts={}):
        '''Loads options specific to this config file (from 'OPTIONS' section).
        External options (externalOpts) can be provided but will overwrite those in
        the config (and modify the config in-place so that the version later saved
        reflects the conditions under which the config was used).

        Args:
            externalOpts (dict, optional): Option-value pairs. Defaults to {}.

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
        parser.add_argument('overwrite', default=False, type=bool, nargs='?',
            help='Overwrite any existing files.')
        
        parser.add_argument('ySlices', default=[], type=list, nargs='?',
            help='Manually define the slices in the y-axis for the sake of plotting. Only needed if the automated algorithm is not working as intended. Defaults to empty list.')
        parser.add_argument('year', default=1, type=int, nargs='?',
            help='Year information used for the sake of plotting text. Defaults to 1 which indicates that the full Run 2 is being analyzed.')
        parser.add_argument('recycle', default=[], type=list, nargs='?',
            help='List of items to recycle. Not currently working.')
        
        return parse_arg_dict(parser,self.Section('OPTIONS'))

    def SaveOut(self): # pragma: no cover
        '''Save two objects to the <self.projPath> directory:
        - a copy of the manipulated config (runConfig.json)
        - the pickled histogram map (hist_map.p)
        '''
        file_out = open(self.projPath+'runConfig.json', 'w')
        json.dump(self.config,file_out,indent=2,sort_keys=True)
        file_out.close()
        # pickle.dump(self.organized_hists, open(self.projPath+'hist_map.p','wb'))
        # self.workspace.writeToFile(self.projPath+'base_'+self.name+'.root',True)  

    def FullTable(self):
        '''Generate full table of processes, regions, and systematic variations
        to account for, including relevant information for each. The table is
        returned as a pandas DataFrame for convenient manipulation.

        Returns:
            pandas.DataFrame: Table
        '''
        regions = self._regionTable()
        processes = self._processTable()
        systematics = self._systematicsTable()

        proc_syst = processes.merge(systematics,right_index=True,left_on='variation',how='left',suffixes=['','_syst'])
        proc_syst = _df_condense_nameinfo(proc_syst,'source_histname')
        proc_syst = _df_condense_nameinfo(proc_syst,'source_filename')

        final = regions.merge(proc_syst,right_index=True,left_on='process',how='left')
        final = _keyword_replace(final, ['source_filename', 'source_histname'])
        _df_sanity_checks(final)
        return final

    def _regionTable(self):
        '''Generate the table of region information based on the JSON config.
        Columns are `process` and `region`.

        Returns:
            pandas.DataFrame
        '''
        def _data_not_included(region):
            '''Check if the list of processes associated with a region
            includes the one marked as `type` `DATA` in the `PROCESSES`
            section of the config. If it doesn't included it but it exists,
            return the name so it can be added to the list of processes
            for the region.

            Args:
                region (str): Name of region to check.

            Raises:
                RuntimeError: No process of TYPE == "DATA" provided.
                RuntimeError: Multiple processes of TYPE == "DATA" provided in PROCESSES section.
                RuntimeError: Multiple processes of TYPE == "DATA" provided in REGIONS subsection.

            Returns:
                str: The name of the data key if it's not already included in the list of processes
                for the region.
            '''
            region_df = pandas.DataFrame({'process':self.Section('REGIONS')[region]})
            process_df = pandas.DataFrame(self.Section('PROCESSES')).T[['TYPE']]
            region_df = region_df.merge(process_df,
                                        left_on='process',
                                        right_index=True,
                                        how='left')

            # Check DATA type is even provided in the PROCESSES
            if (process_df['TYPE'] == 'DATA').sum() == 0:
                raise RuntimeError('No process of TYPE == "DATA" provided.')
            elif (process_df['TYPE'] == 'DATA').sum() > 1:
                raise RuntimeError('Multiple processes of TYPE == "DATA" provided in PROCESSES section.')
            else:
                data_key = process_df[process_df['TYPE'] == 'DATA'].index[0]
            # Check if it's included in the regions
            if (region_df['TYPE'] == 'DATA').sum() == 0:
                out = data_key
            elif ((region_df['TYPE'] == 'DATA').sum() == 1):
                out = False
            else:
                raise RuntimeError('Multiple processes of TYPE == "DATA" provided in REGIONS subsection.')

            return out

        out_df = pandas.DataFrame(columns=['process','region'])
        for r in self.Section('REGIONS'):
            data_key = _data_not_included(r)
            if data_key:
                out_df = out_df.append(pandas.Series({'process':data_key,'region':r}),ignore_index=True)

            for p in self.Section('REGIONS')[r]:
                out_df = out_df.append(pandas.Series({'process':p,'region':r}),ignore_index=True)
            
        return out_df

    def _processTable(self):
        '''Generate the table of process information based on the JSON config.
        Columns are `color`, `process_type`, `scale`, `variation`, `source_filename`,
        and `source_histname`.

        Returns:
            pandas.DataFrame
        '''
        out_df = pandas.DataFrame(columns=['color','process_type','scale','variation','source_filename','source_histname'])
        for p in self.Section('PROCESSES'):
            this_proc_info = self.Section('PROCESSES')[p]
            for s in this_proc_info['SYSTEMATICS']+['nominal']:
                out_df = out_df.append(pandas.Series(
                                {'color': nan if 'COLOR' not in this_proc_info else this_proc_info['COLOR'],
                                'process_type': this_proc_info['TYPE'],
                                'scale': 1.0 if 'SCALE' not in this_proc_info else this_proc_info['SCALE'],
                                'source_filename': this_proc_info['LOC'].split(':')[0],
                                'source_histname': this_proc_info['LOC'].split(':')[1],
                                'variation': s},
                                name=p)
                            )
        return out_df

    def _systematicsTable(self):
        '''Generate the table of process information based on the JSON config.
        Columns are  'lnN', 'lnN_asym', 'shape_sigma', 'syst_type', 'source_filename',
        'source_histname', and 'direction' (ie. NaN, 'Up', or 'Down').

        Returns:
            pandas.DataFrame
        '''
        out_df = pandas.DataFrame(columns=_syst_col_defaults.keys())
        for s in self.Section('SYSTEMATICS'):
            for syst in _get_syst_attrs(s,self.Section('SYSTEMATICS')[s]):
                out_df = out_df.append(syst)
        return out_df

    def get_hist_map(self,df=None):
        '''Collect information on the histograms to extract, manipulate, and save
        into organized_hists.root and store it inside a `dict` where the key is the
        filename and the value is a DataFrame with columns `source_histname`, `out_histname`,
        `scale`, and `color`.

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
                return row.process+'_'+row.region+'_'+row.variation+"_FULL"
            else:
                return row.process+'_'+row.region+'_'+row.variation+row.direction+"_FULL"

        if not isinstance(df,pandas.DataFrame):
            df = self.df

        hists = {}
        for g in df.groupby(['source_filename']):
            group_df = g[1].copy()
            group_df = group_df[group_df['variation'].eq('nominal') | group_df["syst_type"].eq("shapes")]
            group_df['out_histname'] = group_df.apply(_get_out_name,axis=1)
            hists[g[0]] = group_df[['source_histname','out_histname','scale','color']]
        return hists

class OrganizedHists():
    '''Class to store histograms in a consistent data structure and with accompanying
    methods to access the histograms.

    Attributes:
        name (str): Name, taken from input configObj.
        filename (str): Path to `organized_hists.root`.
        hists (dict): Three-level nested dictionary organized as [process][region][systematic variation].
        binning (Binning): Binning object, taken from configObj.
        rebinned (bool): Flag to denote if a rebinning has already occured.
        file (ROOT.TFile): TFile to store histograms on disk.

    Args:
        configObj (Config): Config object.
    '''
    def __init__(self,configObj,readOnly=False):
        self.name = configObj.name
        self.filename = configObj.projPath + 'organized_hists.root'
        self.binning = configObj.binning
        self.hist_map = configObj.get_hist_map()

        if os.path.exists(self.filename) and readOnly:
            self.file = ROOT.TFile.Open(self.filename,"OPEN")
        else:
            self.file = ROOT.TFile.Open(self.filename,"RECREATE")
            self.Add()

    def Add(self):
        '''Manipulate all histograms in self.hist_map and save them to organized_hists.root.

        Returns:
            None
        '''
        for infilename,histdf in self.hist_map.items():
            infile = ROOT.TFile.Open(infilename)
            for row in histdf.itertuples():
                h = infile.Get(row.source_histname)
                h.SetDirectory(0)
                h.Scale(row.scale)

                if get_bins_from_hist("Y", h) != self.binning.ybinList:
                    h = copy_hist_with_new_bins(row.out_histname+'_rebinY','Y',h,self.binning.ybinList)
                if get_bins_from_hist("X", h) != self.binning.xbinList:
                    h = copy_hist_with_new_bins(row.out_histname,'X',h,self.binning.xbinList)
                else:
                    h.SetName(self.info['new_name']+'_FULL')

                h.SetTitle(row.out_histname)
                h.SetFillColor(row.color)

                if h.Integral() <= 0:
                    print ('WARNING: %s has zero or negative events - %s'%(row.out_histname, h.Integral()))
                    for b in range(1,h.GetNbinsX()*h.GetNbinsY()+1):
                        h.SetBinContent(b,1e-10)

                self.file.WriteObject(h, row.out_histname)

                self.CreateSubRegions(h)

            infile.Close()

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

    def CreateSubRegions(self,h):
        '''Sub-divide input histogram along the X axis into the regions specified in the config
        and write the new histogram to organized_hists.root.

        Returns:
            None
        '''
        for sub in self.binning.xbinByCat.keys():
            hsub = h.Clone()
            hsub = copy_hist_with_new_bins(h.GetName().replace('_FULL','_'+sub),'X',h,self.binning.xbinByCat[sub])
            hsub.SetTitle(hsub.GetName())
            self.file.WriteObject(hsub, hsub.GetName())

def _keyword_replace(df,col_strs):
    '''Given a DataFrame and list of column names,
    find and replace the three keywords ("$process", "$region$", "$syst") with their
    respective values in the row for the DataFrame.

    Args:
        df (pandas.DataFrame): DataFrame in which to do the find-replace and to find the per-row keyword matches.
        col_strs (list(str)): List of column names to consider for the find-replace.

    Returns:
        pandas.DataFrame: The manipulated DataFrame copy.
    '''
    def _batch_replace(row,s=None):
        if pandas.isna(row[s]):
            return nan
        else:
            return replace_multi(
                row[col_str],
                {'$process': row.process,
                 '$region':  row.region,
                 '$syst':    row.variation}
            )

    for col_str in col_strs:
        df[col_str] = df.apply(_batch_replace, axis='columns', s=col_str)
    return df

def _get_syst_attrs(name,syst_dict):
    '''Parse an entry in the `"SYSTEMATICS"` section of the JSON config and
    generate the row(s) to append to the systematics DataFrame based on
    that information.

    Args:
        name (str): Name of the systematic variation.
        syst_dict (dict): Dictionary of the config["SYSTEMATICS"][name] section of the JSON config.

    Raises:
        RuntimeError: Systematic variation type could not be determined.

    Returns:
        list(pands.Series): List of new rows to append to the main systematics DataFrame.
    '''
    if 'VAL' in syst_dict:
        out = [{
            'lnN':syst_dict['VAL'],
            'syst_type': 'lnN'
        }]
    elif 'VALUP' in syst_dict and 'VALDOWN' in syst_dict:
        out = [{
            'lnN_asym':[syst_dict['VALUP'], syst_dict['VALDOWN']],
            'syst_type': 'lnN_asym'
        }]
    elif 'UP' in syst_dict and 'DOWN' in syst_dict:
        out = [
            {
                'shape_sigma':syst_dict['SIGMA'],
                'syst_type': 'shapes',
                'source_filename': syst_dict['UP'].split(':')[0],
                'source_histname': syst_dict['UP'].split(':')[1],
                'direction': 'Up'
            }, {
                'shape_sigma':syst_dict['SIGMA'],
                'syst_type': 'shapes',
                'source_filename': syst_dict['DOWN'].split(':')[0],
                'source_histname': syst_dict['DOWN'].split(':')[1],
                'direction': 'Down'
            }
        ]
    else:
        raise RuntimeError('Systematic variation type could not be determined for "%s".'%name)

    out = [pandas.Series(copy_update_dict(_syst_col_defaults, d), name=name) for d in out]
    return out

def _df_condense_nameinfo(df,baseColName):
    '''Condense information after the left-join of the process and systematics DataFrames which creates duplicates
    of the `source_*` columns. Manipulate `df` and return it so that `baseColName+"_syst"` replaces `baseColName`
    and is then dropped.

    Args:
        df (pandas.DataFrame): Input DataFrame to condense.
        baseColName (str): Name of column to condense into.

    Returns:
        pandas.DataFrame: Condensed DataFrame.
    '''
    df[baseColName] = df.apply(lambda row: row[baseColName+'_syst'] if pandas.notna(row[baseColName+'_syst']) else row[baseColName], axis='columns')
    df.drop(baseColName+'_syst',axis='columns',inplace=True)
    return df

def _df_sanity_checks(df):
    '''Check for duplicated (process,region,variation,source_filename,source_histname).
    Prints any duplicate rows if they exist (and raises RuntimeError).

    Args:
        df (pandas.DataFrame): DataFrame to check.

    Raises:
        RuntimeError: Duplicates exist (duplicated rows are printed).
    '''
    # check for duplicate process+region+variation
    dupes = df[df.duplicated(subset=['process','region','variation','source_filename','source_histname'],keep=False)]
    if dupes.shape[0] > 0:
        raise RuntimeError('Duplicates exist. Printing them...\n%s'%dupes)
      
def config_loop_replace(config,old,new,inGLOBAL=False):
    '''Self-calling function loop to find-replace one pair (old,new)
    in a nested dictionary or list (config). If old, new, and the config entry
    examined are all strings, all instances of <old> will be replaced by <new> in the
    string in the config. If this condition is not satisified, <old> must match the config
    entry in its entirety (ie. <old> == <config dict value>).

    Args:
        config (dict,list): Nested dictionary or list where keys and values will have the
            find-replace algorithm applied.
        old (non-nested obj): Object to replace (of type string, int, float, etc - no lists or dictionaries).
        new (non-nested obj): Object replacement (of type string, int, float, etc) - no lists or dictionaries).

    Raises:
        TypeError: If input is not a dict or list.

    Returns:
        dict,list: Modified dict/list.
    '''
    next_is_global = False
    if isinstance(config,dict):
        for k,v in config.items():
            if k == 'GLOBAL':
                next_is_global = True
            if old in k and not inGLOBAL:
                config[re.sub(r'\b%s\b'%old,new,k)] = config.pop(k)
                k = re.sub(r'\b%s\b'%old,new,k)
            if isinstance(v,dict) or isinstance(v,list):
                config[k] = config_loop_replace(v,old,new,inGLOBAL=next_is_global)
            elif isinstance(v,str) and isinstance(old,str) and isinstance(new,str):
                if old in v:
                    config[k] = re.sub(r'\b%s\b'%old,new,v)
            else:
                if old == v:
                    config[k] = new
    elif isinstance(config,list):
        for i,v in enumerate(config):
            if isinstance(v,dict) or isinstance(v,list):
                config[i] = config_loop_replace(v,old,new)
            elif isinstance(v,str):
                if old in v:
                    config[i] = re.sub(r'\b%s\b'%old,new,v)
            else:
                if old == v:
                    config[i] = new
    else:
        raise TypeError('Type "%s" not accepted in config_loop_replace.')

    return config
