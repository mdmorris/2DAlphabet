import argparse, os, warnings, itertools
from TwoDAlphabet.config import Config
from TwoDAlphabet.helpers import execute_cmd, get_config_dirs, parse_arg_dict
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
        self.SetupProjDir()
        if self.options.draw == False:
            ROOT.gROOT.SetBatch(True)

        # self.tf = 

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

    def SaveOut(self):
        '''Save individual configs to project directory.
        '''
        for c in self.configs.values():
            c.SaveOut()

    def MakeCard(self):
        card_new = open(self.projPath+'card_'+self.name+'.txt','w')

        # imax (bins), jmax (backgrounds), kmax (systematics) 
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

    def MakeWorkspace(self):
        self.floatingBins = [] # This holds the names of all of the variables that we want to float.
                           # These are typically bins in the RPH2D 

        ################################
        # Establish our axis variables #
        ################################
        x_vars,y_var = self._getRRVs()  # x_vars is a dict with different RRVs for LOW,SIG,HIGH (keys)
        self.allVars.extend([x_vars,y_var])
        var_lists = {}
        for c in x_vars.keys():
            var_lists[c] = RooArgList(x_vars[c],y_var)

        #########################
        #   Make RooDataHists   #
        #########################
        # It may have seemed crazy to keep this dictionary of TH2s around but it has two great things
        # 1 - structure, 2 - the TH2s we need to make into RDHs
        # However, we will do one thing for convenience - copy it and replace the TH2s in the copy with RDHs
        # if the process has CODE 0,1,2 and a PDF with a normalization if the CODE is 3

        Roo_dict = header.dictCopy(self.organizedDict)

        # For procees, cat, dict...
        for process in self.organizedDict.keys():
            if process == 'qcdmc': continue
            for cat in ['pass','fail']:
                for c in ['LOW','SIG','HIGH']:
                    for dist in self.organizedDict[process][cat+'_'+c].keys():
                        # For each category
                        Roo_dict[process][cat+'_'+c][dist] = {}
                        var_list = var_lists[c]
                        Roo_dict[process][cat+'_'+c][dist] = {}
                        print ('Making RDH '+self.organizedDict[process][cat+'_'+c][dist])
                        Roo_dict[process][cat+'_'+c][dist]['RDH'] = header.makeRDH(self.orgFile.Get(self.organizedDict[process][cat+'_'+c][dist]),var_list)


        #############################################################################################
        # Everything from here on is only dealing with the QCD estimate - everything else is done   #
        #############################################################################################
                    
        ######################################
        # Build the RooParametricHist2D bins #
        ######################################
        Roo_dict['qcd'] = {}
        for r in ['pass','fail']:
            for c in ['LOW','SIG','HIGH']:
                Roo_dict['qcd'][r+'_'+c] = {}
        
        if self.rpfRatio != False:
            TH2_qcdmc_ratios = {}
            TH2_qcdmc_fail = self.orgFile.Get(self.organizedDict['qcdmc']['fail_FULL']['nominal'])
            TH2_qcdmc_pass = self.orgFile.Get(self.organizedDict['qcdmc']['pass_FULL']['nominal'])

            TH2_qcdmc_ratios['FULL'] = TH2_qcdmc_pass.Clone('qcdmc_rpf_full')
            TH2_qcdmc_ratios['FULL'].Divide(TH2_qcdmc_fail)
            for c in ['LOW','SIG','HIGH']:
                TH2_qcdmc_ratios[c] = header.copyHistWithNewXbins(TH2_qcdmc_ratios['FULL'],self.newXbins[c],'qcdmc_rpf_'+c+'_smooth')
        
        TH2_data_toy_ratios = {}
        TH2_data_pass_toys = {}
        TH2_data_fail_toys = {}
        # Need to build for each category
        for c in ['LOW','SIG','HIGH']:
            bin_list_fail = RooArgList()
            bin_list_pass = RooArgList()

            TH2_data_fail = self.orgFile.Get(self.organizedDict['data_obs']['fail_'+c]['nominal'])
            TH2_data_pass = self.orgFile.Get(self.organizedDict['data_obs']['pass_'+c]['nominal'])
            
            TH2_data_fail_toy = TH2_data_fail.Clone()
            TH2_data_pass_toy = TH2_data_pass.Clone()

            for process in self.organizedDict.keys():
                if process == 'qcdmc': continue
                elif self.inputConfig['PROCESS'][process]['CODE'] == 2: 
                    to_subtract_fail = self.orgFile.Get(self.organizedDict[process]['fail_'+c]['nominal'])
                    to_subtract_pass = self.orgFile.Get(self.organizedDict[process]['pass_'+c]['nominal'])
                    
                    TH2_data_fail_toy.Add(to_subtract_fail,-1)
                    TH2_data_pass_toy.Add(to_subtract_pass,-1)

            
            if self.rpfRatio != False:
                TH2_data_toy_ratios[c] = TH2_data_pass_toy.Clone()
                TH2_data_toy_ratios[c].Divide(TH2_data_fail_toy)
                
            else:
                TH2_data_fail_toys[c] = TH2_data_fail_toy
                TH2_data_pass_toys[c] = TH2_data_pass_toy

            # Get each bin
            for ybin in range(1,len(self.newYbins)):
                for xbin in range(1,len(self.newXbins[c])):
                    this_full_xbin = self._getFullXbin(xbin,c)
                    # Now that we're in a specific bin, we need to process it

                    # Make a name for the bin RRV
                    fail_bin_name = 'Fail_bin_'+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name
                    pass_bin_name = 'Pass_bin_'+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name

                    # Initialize contents
                     # First check if we want a parametric fail
                    if self.parametricFail != False and (self.newYbins[ybin-1] > self.parametricFail['START']):
                        binRRV = self.fail_func.Eval(this_full_xbin,ybin,fail_bin_name)
                        # Store the bin
                        bin_list_fail.add(binRRV)
                        self.allVars.append(binRRV)

                        # And now get the Rpf function value for this bin 
                        this_rpf = self.rpf.Eval(this_full_xbin,ybin)

                        if self.rpfRatio == False:
                            formula_arg_list = RooArgList(binRRV,this_rpf)
                            this_bin_pass = RooFormulaVar(pass_bin_name, pass_bin_name, "@0*@1",formula_arg_list)
                            
                        else:
                            mc_ratio_var = RooConstVar("mc_ratio_x_"+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name, "mc_ratio_x_"+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name, TH2_qcdmc_ratios[c].GetBinContent(xbin,ybin))
                            formula_arg_list = RooArgList(binRRV,this_rpf,mc_ratio_var)
                            this_bin_pass = RooFormulaVar(pass_bin_name, pass_bin_name, "@0*@1*@2",formula_arg_list)
                            self.allVars.append(mc_ratio_var)

                        bin_list_pass.add(this_bin_pass)
                        self.allVars.append(formula_arg_list)
                        self.allVars.append(this_bin_pass)
                        self.allVars.append(this_rpf)

                    else:
                        bin_content    = TH2_data_fail.GetBinContent(xbin,ybin)
                        bin_range_up   = bin_content*3 
                        bin_range_down = 1e-9
                        bin_err_up     = TH2_data_fail.GetBinErrorUp(xbin,ybin)
                        bin_err_down   = TH2_data_fail.GetBinErrorLow(xbin,ybin)

                        # Now subtract away the MC
                        for process in self.organizedDict.keys():
                            this_TH2 = self.orgFile.Get(self.organizedDict[process]['fail_'+c]['nominal'])

                            # Check the code and change bin content and errors accordingly
                            if process == 'qcdmc': continue
                            elif self.inputConfig['PROCESS'][process]['CODE'] == 0: continue # signal
                            elif self.inputConfig['PROCESS'][process]['CODE'] == 1: continue # data
                            elif self.inputConfig['PROCESS'][process]['CODE'] == 2: # MC
                                bin_content    = bin_content     - this_TH2.GetBinContent(xbin,ybin)
                                bin_err_up     = bin_err_up      + this_TH2.GetBinErrorUp(xbin,ybin) #- this_TH2.GetBinContent(xbin,ybin)             # Just propagate errors normally
                                bin_err_down   = bin_err_down    - this_TH2.GetBinErrorLow(xbin,ybin) #- this_TH2.GetBinContent(xbin,ybin)

                        # If fail bin content is <= 0, treat this bin as a RooConstVar at value close to 0
                        if (bin_content <= 0):# or (this_pass_bin_zero == True):
                            binRRV = RooConstVar(fail_bin_name, fail_bin_name, max(1e-9,bin_content))
                            bin_list_fail.add(binRRV)
                            self.allVars.append(binRRV)

                            # Now get the Rpf function value for this bin 
                            self.rpf.Eval(this_full_xbin,ybin) # store rpf for this bin but dont need return

                            this_bin_pass = RooConstVar(pass_bin_name, pass_bin_name, 1e-9)
                            bin_list_pass.add(this_bin_pass)
                            self.allVars.append(this_bin_pass)

                        else:
                            # Create the fail bin
                            if self.freezeFail:
                                binRRV = RooConstVar(fail_bin_name, fail_bin_name, bin_content)

                            else:
                                if bin_content < 1: # Give larger floating to range to bins with fewer events
                                    binRRV = RooRealVar(fail_bin_name, fail_bin_name, max(bin_content,0.1), 1e-9, 10)
                                    print (fail_bin_name + ' < 1')

                                elif bin_content < 10: # Give larger floating to range to bins with fewer events
                                    binRRV = RooRealVar(fail_bin_name, fail_bin_name, max(bin_content,1), 1e-9, 50)
                                    print (fail_bin_name + ' < 10')
                                else:
                                    binRRV = RooRealVar(fail_bin_name, fail_bin_name, bin_content, max(1e-9,bin_range_down), bin_range_up)

                                if bin_content - bin_err_down < 0.0001:
                                    bin_err_down = bin_content - 0.0001 # For the edge case when bin error is larger than the content
                                
                                binRRV.setAsymError(bin_err_down,bin_err_up)
                                self.floatingBins.append(fail_bin_name)

                            # Store the bin
                            bin_list_fail.add(binRRV)
                            self.allVars.append(binRRV)

                            # And now get the Rpf function value for this bin 
                            this_rpf = self.rpf.Eval(this_full_xbin,ybin)

                            if self.rpfRatio == False:
                                formula_arg_list = RooArgList(binRRV,this_rpf)
                                this_bin_pass = RooFormulaVar(pass_bin_name, pass_bin_name, "@0*@1",formula_arg_list)
                                
                            else:
                                mc_ratio_var = RooConstVar("mc_ratio_x_"+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name, "mc_ratio_x_"+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name, TH2_qcdmc_ratios[c].GetBinContent(xbin,ybin))
                                formula_arg_list = RooArgList(binRRV,this_rpf,mc_ratio_var)
                                this_bin_pass = RooFormulaVar(pass_bin_name, pass_bin_name, "@0*@1*@2",formula_arg_list)
                                self.allVars.append(mc_ratio_var)

                            bin_list_pass.add(this_bin_pass)
                            self.allVars.append(formula_arg_list)
                            self.allVars.append(this_bin_pass)
                            self.allVars.append(this_rpf)


            print ("Making RPH2Ds")
            Roo_dict['qcd']['fail_'+c] = {}
            Roo_dict['qcd']['pass_'+c] = {}

            Roo_dict['qcd']['fail_'+c]['RPH2D'] = RooParametricHist2D('qcd_fail_'+c+'_'+self.name,'qcd_fail_'+c+'_'+self.name,x_vars[c], y_var, bin_list_fail, TH2_data_fail)
            Roo_dict['qcd']['fail_'+c]['norm']  = RooAddition('qcd_fail_'+c+'_'+self.name+'_norm','qcd_fail_'+c+'_'+self.name+'_norm',bin_list_fail)
            Roo_dict['qcd']['pass_'+c]['RPH2D'] = RooParametricHist2D('qcd_pass_'+c+'_'+self.name,'qcd_pass_'+c+'_'+self.name,x_vars[c], y_var, bin_list_pass, TH2_data_fail)
            Roo_dict['qcd']['pass_'+c]['norm']  = RooAddition('qcd_pass_'+c+'_'+self.name+'_norm','qcd_pass_'+c+'_'+self.name+'_norm',bin_list_pass)

        if self.rpfRatio != False:
            mc_rpf = TH2_qcdmc_ratios['FULL']#header.stitchHistsInX('mc_ratio',self.fullXbins,self.newYbins,[TH2_qcdmc_ratios['LOW'],TH2_qcdmc_ratios['SIG'],TH2_qcdmc_ratios['HIGH']])
            data_rpf = header.stitchHistsInX('data_ratio',self.fullXbins,self.newYbins,[TH2_data_toy_ratios['LOW'],TH2_data_toy_ratios['SIG'],TH2_data_toy_ratios['HIGH']],blinded=[1] if self.blindedPlots else [])
            rpf_ratio = data_rpf.Clone()
            rpf_ratio.Divide(mc_rpf)
            rpf_ratio.SetMaximum(2.5)
            # rpf_ratio.GetZaxis().SetLabelSize(0.04)

            header.makeCan('rpf_ratio',self.projPath,[data_rpf,mc_rpf,rpf_ratio],
                titles=["Data Ratio;%s;%s;R_{P/F}"%(self.xVarTitle,self.yVarTitle),"MC Ratio;%s;%s;R_{P/F}"%(self.xVarTitle,self.yVarTitle),"Ratio of ratios;%s;%s;R_{Ratio}"%(self.xVarTitle,self.yVarTitle)],
                year=self.year)
        else: 
            data_fail = header.stitchHistsInX('data_fail',self.fullXbins,self.newYbins,[TH2_data_fail_toys['LOW'],TH2_data_fail_toys['SIG'],TH2_data_fail_toys['HIGH']],blinded=[1] if self.blindedPlots else [])
            data_pass = header.stitchHistsInX('data_fail',self.fullXbins,self.newYbins,[TH2_data_pass_toys['LOW'],TH2_data_pass_toys['SIG'],TH2_data_pass_toys['HIGH']],blinded=[1] if self.blindedPlots else [])
            data_rpf = data_pass.Clone()
            data_rpf.Divide(data_fail)
            header.makeCan('data_ratio',self.projPath,[data_pass,data_fail,data_rpf],titles=["Data Pass","Data Fail","R_{P/F}"],year=self.year)
            header.makeCan('data_ratio_lego',self.projPath,[data_pass,data_fail,data_rpf],titles=["Data Pass","Data Fail","R_{P/F}"],year=self.year,datastyle='lego')

        # Determine whether we need rpfRatio templates
        if self.rpfRatio != False and self.rpfRatioVariations != False:
            for v in self.rpfRatioVariations:
                TH2_qcdmc_fail = self.orgFile.Get(self.organizedDict['qcdmc']['fail_FULL']['KDEbandwidth'+v.capitalize()])
                TH2_qcdmc_pass = self.orgFile.Get(self.organizedDict['qcdmc']['pass_FULL']['KDEbandwidth'+v.capitalize()])

                TH2_qcdmc_ratios['FULL_'+v] = TH2_qcdmc_pass.Clone('qcdmc_rpf_full_'+v)
                TH2_qcdmc_ratios['FULL_'+v].Divide(TH2_qcdmc_fail)

                for c in ['LOW','SIG','HIGH']:
                    TH2_qcdmc_ratios[c+'_'+v] = header.copyHistWithNewXbins(TH2_qcdmc_ratios['FULL_'+v],self.newXbins[c],'qcdmc_rpf_'+c+'_'+v)
                    bin_list_fail = Roo_dict['qcd']['fail_'+c]['RPH2D'].getPars()
                    bin_list_pass = RooArgList()
                    for ybin in range(1,len(self.newYbins)):
                        for xbin in range(1,len(self.newXbins[c])):
                            this_full_xbin = self._getFullXbin(xbin,c)
                            # Now that we're in a specific bin, we need to process it
                            fail_bin_name = 'Fail_bin_'+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name
                            pass_bin_name = 'Pass_bin_'+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name+'_'+v

                            binRRV = bin_list_fail.find(fail_bin_name)
                            bin_content = binRRV.getValV()
                            # If fail bin content is <= 0, treat this bin as a RooConstVar at value close to 0
                            if isinstance(binRRV,RooConstVar):# or (this_pass_bin_zero == True):
                                this_bin_pass = RooConstVar(pass_bin_name, pass_bin_name, 1e-9)
                                bin_list_pass.add(this_bin_pass)
                                self.allVars.append(this_bin_pass)

                            else:
                                # And now get the Rpf function value for this bin 
                                this_rpf = self.rpf.getFuncBinRRV(c,this_full_xbin,ybin)
                                mc_ratio_var = RooConstVar("mc_ratio_"+v+"_x_"+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name, 
                                                           "mc_ratio_"+v+"_x_"+str(this_full_xbin)+'-'+str(ybin)+'_'+self.name, 
                                                           TH2_qcdmc_ratios[c+'_'+v].GetBinContent(xbin,ybin))
                                formula_arg_list = RooArgList(binRRV,this_rpf,mc_ratio_var)
                                this_bin_pass = RooFormulaVar(pass_bin_name, pass_bin_name, "@0*@1*@2",formula_arg_list)
                                
                                bin_list_pass.add(this_bin_pass)
                                self.allVars.append(formula_arg_list)
                                self.allVars.append(this_bin_pass)
                                self.allVars.append(mc_ratio_var)

                    Roo_dict['qcd']['pass_'+c+'_'+v] = {}
                    Roo_dict['qcd']['pass_'+c+'_'+v]['RPH2D'] = RooParametricHist2D('qcd_pass_'+c+'_'+self.name+'_KDEbandwidth'+v.capitalize(),
                                                                                    'qcd_pass_'+c+'_'+self.name+'_KDEbandwidth'+v.capitalize(),
                                                                                    x_vars[c], y_var, bin_list_pass, TH2_qcdmc_ratios[c+'_'+v])
                    Roo_dict['qcd']['pass_'+c+'_'+v]['norm']  = RooAddition('qcd_pass_'+c+'_'+self.name+'_KDEbandwidth'+v.capitalize()+'_norm',
                                                                            'qcd_pass_'+c+'_'+self.name+'_KDEbandwidth'+v.capitalize()+'_norm',
                                                                            bin_list_pass)

        print ("Making workspace...")
        # Make workspace to save in
        self.workspace = RooWorkspace("w_"+self.name)
        for process in Roo_dict.keys():
            for cat in [k for k in Roo_dict[process].keys() if 'file' not in k and 'FULL' not in k]:
                if process == 'qcd':
                    rooObj = Roo_dict[process][cat]
                    for itemkey in rooObj.keys():
                        print ("Importing " + rooObj[itemkey].GetName() + ' from ' + process + ', ' + cat + ', ' + itemkey)
                        getattr(self.workspace,'import')(rooObj[itemkey],RooFit.RecycleConflictNodes(),RooFit.Silence())
                elif process == 'qcdmc':
                    continue
                else:
                    for dist in Roo_dict[process][cat].keys():
                        rooObj = Roo_dict[process][cat][dist]
                        for itemkey in rooObj.keys():
                            print ("Importing " + rooObj[itemkey].GetName() + ' from ' + process + ', ' + cat  +', ' +dist+ ', ' + itemkey)
                            getattr(self.workspace,'import')(rooObj[itemkey],RooFit.RecycleConflictNodes(),RooFit.Silence())

    def _makeCard(self):
        # card_locs = []
        # for c in self.all_configs:
        #     card_locs.append(c.makeCard())
        # combine_cards_cmd = 'combineCards.py ... %s'%(' '.join(card_locs))
        # execute_cmd(combine_cards_cmd)
        pass

    def plot(self):
        # plotter.methodA()
        # plotter.methodB()
        # ...
        pass

def proj_info_lookup(projDir,card_tag):
    # Check if there was more than one 2DAlphabet object run over
    more_than_one = False
    run_card = open(projDir+'/card_'+card_tag+'.txt','r')
    firstline = run_card.readline()
    if 'Combination of ' in firstline:
        more_than_one = True
    
    proj_info = {}
    if more_than_one:
        # Look up all categories run in the most recent fit
        twoD_names = getTwoDAlphaNames(firstline)
        for n in twoD_names:
            proj_info[n] = pickle.load(open(projDir+'/'+n+'/saveOut.p','r'))
            proj_info[n]['rpfVarNames'] = proj_info[n]['rpf'].getFuncVarNames()

    elif not more_than_one:
        proj_info[card_tag] = pickle.load(open(projDir+'/saveOut.p','r'))
        proj_info['rpfVarNames'] = proj_info['rpf'].getFuncVarNames()

    return proj_info

def GetTwoDAlphaNames(line):
    card_locs = [loc for loc in line.split(' ') if (loc != ' ' and loc != '' and 'txt' in loc)]
    proj_names = [n.split('/')[0] for n in card_locs]
    return proj_names

def SetSnapshot(d=''):
    w_f = TFile.Open(d+'higgsCombineTest.FitDiagnostics.mH120.root')
    w = w_f.Get('w')
    fr_f = TFile.Open(d+'fitDiagnosticsTest.root')
    fr = fr_f.Get('fit_b')
    myargs = RooArgSet(fr.floatParsFinal())
    w.saveSnapshot('initialFit',myargs,True)
    fout = TFile('initialFitWorkspace.root',"recreate")
    fout.WriteTObject(w,'w')
    fout.Close()