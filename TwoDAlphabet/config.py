import ROOT, argparse, json, os, pandas, re, warnings, sys
from numpy import nan
import pprint
pp = pprint.PrettyPrinter(indent=4)
from TwoDAlphabet.helpers import copy_update_dict, open_json, parse_arg_dict, replace_multi
from TwoDAlphabet.binning import Binning, copy_hist_with_new_bins, get_bins_from_hist

_protected_keys = ["PROCESSES","SYSTEMATICS","REGIONS","BINNING","OPTIONS","GLOBAL","SCALE","COLOR","TYPE","X","Y","NAME","TITLE","BINS","NBINS","LOW","HIGH"]
_syst_col_defaults = {
    # 'variation': nan,
    'lnN': nan,
    'shapes': nan, # shape sigma
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
        jsonPath (str): File name and path.
        projPath (str): Project path.
        findreplace (dict, optional): Find-replace pairs. Defaults to {}.
        externalOptions (dict, optional): Extra key-value pairs to add to the OPTIONS
            section of the JSON file. Defaults to {}.
        loadPrevious (bool, optional): Load the previous histograms made instead of remaking. Defaults to False.

    Attributes:
        config (dict): JSON config opened as a dict.
        name (str): Name, inherited from the input JSON config.
        projPath (str): The project path (`tag` attribute in `TwoDAlphabet()`).
        options (Namespace): Final set of options to use.
        loadPrevious (bool): Equal to the constructor arg of the same name.
        nsignals (int): Number of signal processes. Zero before running Construct().
        nbkgs (int): Number of signal processes. Zero before running Construct().
        constructed (bool): Whether Construct() has been called yet or not.
        binning (Binning): Binning object. Only exists after Construct().
        organizedHists (OrganizedHists): Object for storing, manipulating, and accessing histograms for fit.
    '''
    def __init__(self,jsonPath,projPath,findreplace={},externalOptions={},loadPrevious=False):
        self.config = open_json(jsonPath)
        self._addFindReplace(findreplace)
        if 'GLOBAL' in self.config.keys(): self._varReplacement()
        self.name = self.config['NAME']
        self.projPath = projPath+'/'
        self.options = self.GetOptions(externalOptions)
        self.loadPrevious = loadPrevious
        self.nsignals, self.nbkgs = 0, 0
        self.df = self.FullTable()
        self.constructed = False
        self._addedConfigs = []

    def Construct(self):
        '''Uses the config as instructions to setup binning, manipulate histograms,
        save outputs, etc. Separated from the constructor so information from the config
        can be extracted for basic sanity checks before it is processed.
        '''
        self.constructed = True
        template_file = ROOT.TFile.Open(self.df.iloc[0].source_filename)
        template = template_file.Get(self.df.iloc[0].source_histname)
        template.SetDirectory(0)
        self.binnings = {}
        for kbinning in self._section('BINNING').keys():
            self.binnings[kbinning] = Binning(kbinning, self._section('BINNING')[kbinning], template)
        if self.loadPrevious: 
            self.organizedHists = OrganizedHists(self,readOnly=True)
        else:
            self.organizedHists = OrganizedHists(self,readOnly=False)

        return self

    def _section(self,key):
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
            if s in self._section('GLOBAL'):
                raise ValueError('A command line string replacement (%s) conflicts with one already in the configuration file.' %(s))
            self.config['GLOBAL'][s] = findreplace[s]

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
        for old_string in self._section('GLOBAL'):
            if old_string == "HELP":
                print ('WARNING: The HELP entry is deprecated and checking for it will be removed in the future. Please remove it from your config.')
                continue
            new_string = self._section('GLOBAL')[old_string]
            self.config = config_loop_replace(self.config, old_string, new_string)

    def GetOptions(self,externalOptions={}):
        '''Loads options specific to this config file (from 'OPTIONS' section).
        External options (externalOptions) can be provided but will overwrite those in
        the config (and modify the config in-place so that the version later saved
        reflects the conditions under which the config was used).

        Args:
            externalOptions (dict, optional): Option-value pairs. Defaults to {}.

        Returns:
            ArgumentParser.Namespace
        '''
        if 'OPTIONS' in self.config.keys():
            self.config['OPTIONS'].update(externalOptions)

        parser = argparse.ArgumentParser()
        # Blinding
        parser.add_argument('blindedPlots', type=list, nargs='?', required=True,
            help='List of regions in which to blind plots of x-axis SIG. Does not blind fit.')
        parser.add_argument('blindedFit', type=list, nargs='?', required=True,
            help='List of regions in which to blind fit of x-axis SIG. Does not blind plots.')
        # Plotting
        parser.add_argument('haddSignals', default=True, type=bool, nargs='?',
            help='Combine signals into one histogram for the sake of plotting. Still treated as separate in fit. Defaults to True.')
        parser.add_argument('plotTitles', default=False, type=bool, nargs='?',
            help='Include titles in plots. Defaults to False.')
        parser.add_argument('plotUncerts', default=False, type=bool, nargs='?',
            help='Plot comparison of pre-fit uncertainty shape templates in 1D projections. Defaults to False.')
        parser.add_argument('plotPrefitSigInFitB', default=False, type=bool, nargs='?',
            help='In the b-only post-fit plots, plot the signal normalized to its pre-fit value. Defaults to False.')
        parser.add_argument('plotEvtsPerUnit', default=False, type=bool, nargs='?',
            help='Post-fit bins are plotted as events per unit rather than events per bin. Defaults to False.')
        parser.add_argument('year', default=1, type=int, nargs='?',
            help='Year information used for the sake of plotting text. Defaults to 1 which indicates that the full Run 2 is being analyzed.')
        
        if 'OPTIONS' in self.config.keys():
            out = parse_arg_dict(parser,self._section('OPTIONS'))
        else:
            out = parser.parse_args([])

        return out

    def SaveOut(self): # pragma: no cover
        '''Save the histogram table to the <self.projPath> directory in csv
        and markdown formats. No copy of the input JSONs is saved since multiple
        could be provided with the only final combination making it to the histogram table.
        '''
        for c in self._addedConfigs+[self]:
            file_out = open(c.projPath+'runConfig_%s.json'%c.name, 'w')
            json.dump(c.config,file_out,indent=2,sort_keys=True)
            file_out.close()

        if 'index' in self.df.columns:
            self.df = self.df.reset_index(drop=True).drop('index',axis=1)

        self.df.to_csv(self.projPath+'hist_table.csv')
        if sys.version_info.major == 3: self.df.to_markdown(self.projPath+'hist_table.md')

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
        final = _keyword_replace(final, ['source_filename', 'source_histname']).reset_index(drop=True)
        _df_sanity_checks(final)
        return final

    def _regionTable(self):
        '''Generate the table of region information based on the JSON config.
        Columns are `process`, `region`, and `binning`.

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
            region_df = pandas.DataFrame({'process':self._section('REGIONS')[region]['PROCESSES']})
            process_df = pandas.DataFrame(self._section('PROCESSES')).T[['TYPE']]
            region_df = region_df.merge(process_df,
                                        left_on='process',
                                        right_index=True,
                                        how='left')

            # Check DATA type is even provided in the PROCESSES
            if (process_df['TYPE'] == 'DATA').sum() == 0:
                warnings.warn('No process of TYPE == "DATA" provided. Ignoring...', RuntimeWarning)
                data_key = False
            elif (process_df['TYPE'] == 'DATA').sum() > 1:
                raise RuntimeError('Multiple processes of TYPE == "DATA" provided in PROCESSES section.')
            else:
                data_key = 'data_obs'
            # Check if it's included in the regions
            if (region_df['TYPE'] == 'DATA').sum() == 0:
                out = data_key
            elif (region_df['TYPE'] == 'DATA').sum() == 1:
                out = False
            else:
                raise RuntimeError('Multiple processes of TYPE == "DATA" provided in REGIONS subsection.')

            return out

        out_df = pandas.DataFrame(columns=['process','region','binning'])
        for r in self._section('REGIONS'):
            data_key = _data_not_included(r)
            if data_key:
                out_df = out_df.append(pandas.Series({'process':data_key,'region':r, 'binning':self._section('REGIONS')[r]['BINNING']}),ignore_index=True)

            for p in self._section('REGIONS')[r]['PROCESSES']:
                if p not in self._section('PROCESSES'):
                    raise RuntimeError('Process "%s" listed for region "%s" not defined in PROCESSES section.'%(p,r))
                out_df = out_df.append(pandas.Series({'process':p,'region':r,'binning':self._section('REGIONS')[r]['BINNING']}),ignore_index=True)
            
        return out_df

    def _processTable(self):
        '''Generate the table of process information based on the JSON config.
        Columns are `color`, `process_type`, `scale`, `variation`, `source_filename`,
        and `source_histname`.

        Returns:
            pandas.DataFrame
        '''
        out_df = pandas.DataFrame(columns=['color','process_type','scale','variation','source_filename','source_histname','alias','combine_idx'])
        for p in self._section('PROCESSES'):
            this_proc_info = self._section('PROCESSES')[p]
            combine_idx = self._getCombineIdx(this_proc_info)
            if this_proc_info['TYPE'] == 'DATA' and p != 'data_obs':
                raise RuntimeError('Any process of type DATA must have section key "data_obs".')
            for s in this_proc_info['SYSTEMATICS']+['nominal']:
                out_df = out_df.append(pandas.Series(
                                {'color': nan if 'COLOR' not in this_proc_info else this_proc_info['COLOR'],
                                'process_type': this_proc_info['TYPE'],
                                'scale': 1.0 if 'SCALE' not in this_proc_info else this_proc_info['SCALE'],
                                'source_filename': this_proc_info['LOC'].split(':')[0],
                                'source_histname': this_proc_info['LOC'].split(':')[1],
                                'alias': p if 'ALIAS' not in this_proc_info.keys() else this_proc_info['ALIAS'], #in file name
                                'title': p if 'TITLE' not in this_proc_info.keys() else this_proc_info['TITLE'], #in legend entry
                                'variation': s,
                                'combine_idx':combine_idx},
                                name=p)
                            )
        return out_df

    def _systematicsTable(self):
        '''Generate the table of process information based on the JSON config.
        Columns are  'lnN', 'shapes', 'syst_type', 'source_filename',
        'source_histname', and 'direction' (ie. NaN, 'Up', or 'Down').

        Note that 'shapes' is short for 'shape sigma'.

        Returns:
            pandas.DataFrame
        '''
        out_df = pandas.DataFrame(columns=_syst_col_defaults.keys())
        for s in self._section('SYSTEMATICS'):
            for syst in _get_syst_attrs(s,self._section('SYSTEMATICS')[s]):
                out_df = out_df.append(syst)
        return out_df

    def _getCombineIdx(self,procdict):
        if procdict['TYPE'] == 'SIGNAL':
            combine_idx = '-%s'%self.nsignals # first signal idxed at 0 so set it *before* incrementing
            self.nsignals += 1
        else:
            self.nbkgs += 1
            combine_idx = '%s'%self.nbkgs # first signal idxed at 0 so set it *before* incrementing
        return combine_idx

    def GetHistMap(self,df=None):
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
        for g in df.groupby(['source_filename']):
            group_df = g[1].copy()
            group_df = group_df[group_df['variation'].eq('nominal') | group_df["syst_type"].eq("shapes")]
            group_df['out_histname'] = group_df.apply(_get_out_name,axis=1)
            hists[g[0]] = group_df[['source_histname','out_histname','scale','color','binning']]
        return hists

    def Add(self,cNew,onlyOn=['process','region']):
        def _drop(row,dupes_list):
            drop = False
            for d in dupes_list:
                to_compare = []
                for i,c in enumerate(onlyOn): # build bool from onlyOn cols
                    to_compare.append(row[c] == d[i])
                drop = pandas.Series(to_compare).all()

                if drop == True: # if we found a match to a
                    break
            return drop

        if self.constructed == True:
            raise RuntimeError('This config has already been constructed so no additions can be made.')
        if isinstance(onlyOn,str):
            if onlyOn not in ['process','region']:
                raise RuntimeError('Can only add configs together on the "process" or "region" information.')
            onlyOn = [onlyOn]
        elif onlyOn != ['process','region']:
            raise RuntimeError('Can only add configs together on the "process" or "region" information.')
        
        df_modified_base         = self.df.append(cNew.df).reset_index(drop=True)
        df_modified_nominal_only = df_modified_base[df_modified_base.variation.eq('nominal')]
        df_modified_dupes        = df_modified_nominal_only[ df_modified_nominal_only.duplicated(subset=onlyOn,keep='first') ]

        dupes_list = set(zip(*(df_modified_dupes[k] for k in onlyOn)))
        if len(dupes_list) > 0: # if duplicates, replace old with new
            print ('Found duplicates in attempting to modify base Config. Replacing...')
            for d in dupes_list:
                print('\t(%s)'%(','.join(d)))
            df_final = self.df.loc[
                            ~self.df.apply(_drop,args=[dupes_list],axis='columns')
                        ].append(cNew.df).reset_index(drop=True)

        else: # if no duplicates, just use the appended df
            df_final = df_modified_base

        self._addedConfigs.append(cNew)
        self.df = df_final

    def GetNregions(self):
        return self.df.region.nunique()

    def GetNsystematics(self):
        return self.df.variation.nunique()-1

    def InitQCDHists(self):
        '''Loop over all regions and for a given region's data histogram, subtract the list of background histograms,
        and return data-bkgList.

        Returns:
            dict(region,TH2): Dictionary with regions as keys and values as histograms of data-bkgList.
        '''
        out = {}
        for region,group in self.df.groupby('region'):
            data_sources = group.loc[group.process_type.eq('DATA')][['source_filename','source_histname']]
            if data_sources.shape[0] > 1:
                raise RuntimeError('More than one data source found in\n%s'%data_sources)

            data_file = ROOT.TFile.Open(data_sources.iloc[0]['source_filename'])
            data_hist = data_file.Get(data_sources.iloc[0]['source_histname'])
            qcd = data_hist.Clone(data_hist.GetName().replace('data_obs','qcd'))
            qcd.SetDirectory(0)

            bkg_sources = group.loc[group.process_type.eq('BKG') & group.variation.eq('nominal')][['source_filename','source_histname']]
            for bkg_file_name,bkg_hist_names in bkg_sources.groupby('source_filename'):
                if bkg_file_name == data_file.GetName():
                    bkg_file = data_file
                else:
                    bkg_file = ROOT.TFile.Open(bkg_file_name)

                if bkg_hist_names.shape[0] > 1:
                    raise RuntimeError('More than one bkg histogram found in %s\n%s'%(bkg_file_name,bkg_hist_names))

                bkg_hist = bkg_file.Get(bkg_hist_names.iloc[0].source_histname)
                qcd.Add(bkg_hist,-1)

            out[region] = qcd

        return out

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
        self.filename = configObj.projPath + 'organized_hists.root'
        self.binnings = configObj.binnings
        self.hist_map = configObj.GetHistMap()

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
                binning = self.binnings[row.binning]

                if get_bins_from_hist("Y", h) != binning.ybinList:
                    h = copy_hist_with_new_bins(row.out_histname+'_rebinY','Y',h,binning.ybinList)
                if get_bins_from_hist("X", h) != binning.xbinList:
                    h = copy_hist_with_new_bins(row.out_histname,'X',h,binning.xbinList)
                else:
                    h.SetName(row.out_histname)

                h.SetTitle(row.out_histname)
                h.SetFillColor(row.color)

                if h.Integral() <= 0:
                    print ('WARNING: %s has zero or negative events - %s'%(row.out_histname, h.Integral()))
                    for b in range(1,h.GetNbinsX()*h.GetNbinsY()+1):
                        h.SetBinContent(b,1e-10)

                self.file.WriteObject(h, row.out_histname)
                self.CreateSubRegions(h, binning)

            infile.Close()

    def Get(self,histname='',process='',region='',systematic='',subspace='FULL'):
        '''Get histogram from the opened TFile. Specify the histogram
        you want via `histname` or by the combination of `process`, `region`,
        and `systematic` options. The `histname` option will take priority.

        Args:
            histname (str, optional): Name of histogram to get. Overrides other three options if specified. Defaults to ''.
            process (str, optional): Name of process to search for. Must be used in conjunction with `region` and `systematic` options. Overridden by `histname`. Defaults to ''.
            region (str, optional): Name of region to search for. Must be used in conjunction with `process` and `systematic` options. Overridden by `histname`. Defaults to ''.
            systematic (str, optional): Name of systematic to search for. Must be used in conjunction with `process` and `region` options. Overridden by `histname`. Defaults to ''.
            subspace (str, optional): Name of subspace. Default is 'FULL' with other options being 'LOW', 'SIG', and 'HIGH'.

        Raises:
            NameError: If subspace option is not 'FULL','LOW','SIG', or 'HIGH'.

        Returns:
            TH2F: Histogram from file.
        '''
        if subspace not in ['FULL','LOW','SIG','HIGH']:
            raise NameError("Subspace '%s' not accepted. Options are 'FULL','LOW','SIG','HIGH'.")
        return self.file.Get(histname if histname != '' else '_'.join([process,region,subspace,systematic]))

    def GetHistNames(self):
        return [hkey.GetName() for hkey in self.file.GetListOfKeys()]

    def BinningLookup(self,histname):
        all_hists = pandas.concat([v[['out_histname','binning']] for v in self.hist_map.values()])
        return all_hists.loc[all_hists.out_histname.eq(histname)].iloc[0].binning

    def CreateSubRegions(self,h,binning):
        '''Sub-divide input histogram along the X axis into the regions specified in the config
        and write the new histogram to organized_hists.root.

        Returns:
            None
        '''
        for sub in binning.xbinByCat.keys():
            hsub = h.Clone()
            hsub = copy_hist_with_new_bins(h.GetName().replace('_FULL','_'+sub),'X',h,binning.xbinByCat[sub])
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
                {'$process': row.alias,
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
            'lnN':str(syst_dict['VAL']),
            'syst_type': 'lnN'
        }]
    elif 'VALUP' in syst_dict and 'VALDOWN' in syst_dict:
        out = [{
            'lnN':'%s/%s'%(syst_dict['VALDOWN'], syst_dict['VALUP']),
            'syst_type': 'lnN'
        }]
    elif 'UP' in syst_dict and 'DOWN' in syst_dict:
        out = [
            {
                'shapes':syst_dict['SIGMA'],
                'syst_type': 'shapes',
                'source_filename': syst_dict['UP'].split(':')[0],
                'source_histname': syst_dict['UP'].split(':')[1],
                'direction': 'Up'
            }, {
                'shapes':syst_dict['SIGMA'],
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
            next_is_global = False # never consider next k to be GLOBAL
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
