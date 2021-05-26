#####################################################################################################
# 2DAlphabet.py - written by Lucas Corcodilos, 1/1/19                                               #
# ---------------------------------------------------                                               #
# This is the rewrite of the 2DAlphabet.py wrapper which has changed the workflow to a class.       #
# The input is still a properly formatted JSON file. Options have been removed for the sake of      #
# simplicity and are replaced by options in the JSON file.                                          #
#####################################################################################################

#########################################################
#                       Imports                         #
#########################################################
from optparse import OptionParser
import subprocess, CMS_lumi
import cPickle as pickle
import os, sys, array, json
import pprint
pp = pprint.PrettyPrinter(indent = 2)

import header
from header import ConvertToEvtsPerUnit 
import RpfHandler

import ROOT
from ROOT import *

gStyle.SetOptStat(0)

class TwoDAlphabet:
    # If you just want to do everything yourself
    def __init__(self):
        pass

    # def __del__(self):
    #     del self.workspace

    # Initialization setup to just build workspace. All other steps must be called with methods
    def __init__ (self,jsonFileName,quicktag='',recycleAll=False,recycle=False,CLoptions=[],stringSwaps={}): # jsonFileNames is a list
        self.allVars = []    # This is a list of all RooFit objects made. It never gets used for anything but if the
                        # objects never get saved here, then the python memory management will throw them out
                        # because of conflicts with the RooFit memory management. It's a hack.

        # Setup config
        self.jsonFileName = jsonFileName
        self.stringSwaps = stringSwaps
        self.inputConfig = header.openJSON(self.jsonFileName)
        self.CLoptions = CLoptions

        # Setup name
        # NAME - unique to config
        # TAG - used to tie configs together
        if 'name' in self.inputConfig['OPTIONS'].keys():
            self.name = self.inputConfig['OPTIONS']['name']
        else:
            self.name = jsonFileName.split('.json')[0].split('input_')[1]
        if quicktag != '':
            self.tag = quicktag
        elif 'tag' in self.inputConfig['OPTIONS'].keys():
            self.tag = self.inputConfig['OPTIONS']['tag']
        else:
            self.tag = ''

        self.xVarName = self.inputConfig['BINNING']['X']['NAME']
        self.yVarName = self.inputConfig['BINNING']['Y']['NAME']
        self.xVarTitle = self.inputConfig['BINNING']['X']['TITLE']
        self.yVarTitle = self.inputConfig['BINNING']['Y']['TITLE']
        self.sigStart = self.inputConfig['BINNING']['X']['SIGSTART']
        self.sigEnd = self.inputConfig['BINNING']['X']['SIGEND']

        # Setup options
        self.freezeFail = self._getOption('freezeFail')
        self.blindedPlots = self._getOption('blindedPlots')
        self.blindedFit = self._getOption('blindedFit')
        self.plotUncerts = self._getOption('plotUncerts')
        self.draw = self._getOption('draw')
        self.verbosity = str(self._getOption('verbosity')) # will eventually need to be fed as a string so just convert now
        self.fitGuesses = self._getOption('fitGuesses')
        self.prerun = self._getOption('prerun')
        self.overwrite = self._getOption('overwrite')
        self.recycle = recycle if recycle != False else self._getOption('recycle')
        self.rpfRatio = self._getOption('rpfRatio')
        self.parametricFail = self._getOption('parametricFail')
        self.year = self._getOption('year')
        self.plotPrefitSigInFitB = self._getOption('plotPrefitSigInFitB')
        self.plotTitles = self._getOption('plotTitles')
        self.addSignals = self._getOption('addSignals')
        self.plotEvtsPerUnit = self._getOption('plotEvtsPerUnit')
        self.ySlices = self._getOption('ySlices')
        self.recycleAll = recycleAll

        # Setup a directory to save
        self.projPath = self._projPath()

        # Determine whether we need rpfRatio templates
        if self.rpfRatio != False:
            if 'SYSTEMATICS' in self.inputConfig['OPTIONS']['rpfRatio'].keys() and len(self.inputConfig['OPTIONS']['rpfRatio']['SYSTEMATICS']) > 0:
                if len(self.inputConfig['OPTIONS']['rpfRatio']['SYSTEMATICS']) == 1 and self.inputConfig['OPTIONS']['rpfRatio']['SYSTEMATICS'][0] == 'KDEbandwidth':
                    self.rpfRatioVariations = ['up','down']
                else:
                    raise ValueError('Only "KDEbandwidth" is accepted as a systematic uncertainty for the qcd mc Rpf ratio.')
            else:
                self.rpfRatioVariations = False
        else:
            self.rpfRatioVariations = False

        # Check if doing an external import
        importOrgDict = False
        importExternalTemplates = False
        print 'Checking for external imports...'
        for i in self.recycle:
            if 'organizedDict' in i:
                importOrgDict = True
                if '{' in i:
                    importExternalTemplates = i[i.find('{')+1:i.find('}')]
                    print 'Will import organized files from %s'%importExternalTemplates
                else:
                    importExternalTemplates = False

                break

        # Pickle reading if recycling
        if (self.recycle != False and self.recycle != []) or self.recycleAll:
            if not (len(self.recycle) == 1 and importExternalTemplates != False):
                self.pickleFile = pickle.load(open(self.projPath+'saveOut.p','rb'))
        
        # Dict to pickle at the end
        self.pickleDict = {}
 
        # Draw if desired
        if self.draw == False:
            gROOT.SetBatch(kTRUE)

        # Replace signal with command specified

        # Global var replacement
        if not self.recycleAll or 'runConfig' not in self.recycle:
            self._configGlobalVarReplacement()
        else:
            self.inputConfig = self._readIn('runConfig')

        # Get binning for three categories
        if not ('newXbins' in self.recycle and 'newYbins' in self.recycle) and not recycleAll:
            self.newXbins, self.newYbins, self.oldXwidth, self.oldYwidth = self._getBinning(self.inputConfig['BINNING']) # sets self.new*bins (list of user specified bins or False if none specified)
        else:
            self.newXbins = self._readIn('newXbins')
            self.newYbins = self._readIn('newYbins')

        # Get one list of the x bins over the full range
        self.fullXbins = self._getFullXbins(self.newXbins)

        print '\n'
        print self.fullXbins
        print self.newYbins

        # Run pseudo2D for fit guesses and make the config to actually run on
        if ("runConfig" not in self.recycle and not self.recycleAll):
            if self.fitGuesses: self._makeFitGuesses()

        # Initialize rpf class       
        if 'rpf' not in self.recycle and not self.recycleAll:
            self.rpf = RpfHandler.RpfHandler(self.inputConfig['FIT'],self.name,self._dummyTH2(),self.tag)
        else:
            self.rpf = self._readIn('rpf')

        # Initialize parametericFail 
        if self.parametricFail != False:
            self.fail_func = RpfHandler.BinnedParametricFunc(self.inputConfig['OPTIONS']['parametricFail'],self.name,self._dummyTH2(),'failShape')
        else:
            self.fail_func = False

        # Organize everything for workspace building
        if importOrgDict or self.recycleAll:
            if importExternalTemplates != False:
                self.organizedDict = self._readIn('organizedDict',importExternalTemplates)
                header.executeCmd('cp '+importExternalTemplates+'/'+self.name+'/organized_hists.root '+self.projPath+'organized_hists.root')
                # self.orgFile = TFile.Open(importExternalTemplates+'/'+self.name+'/organized_hists.root') # have to save out the histograms to keep them persistent past this function
            else:
                self.organizedDict = self._readIn('organizedDict')
            self.orgFile = TFile.Open(self.projPath+'organized_hists.root') # have to save out the histograms to keep them persistent past this function
        else:
            self.orgFile = TFile(self.projPath+'organized_hists.root','RECREATE') # have to save out the histograms to keep them persistent past this function
            self.organizedDict = {}
            self._inputOrganizer()

        # Make systematic uncertainty plots
        if self.plotUncerts and not self.recycleAll:
            self._makeSystematicPlots()

        # Build the workspace
        if 'workspace' in self.recycle or self.recycleAll:
            self.workspace = self._readIn('workspace')
            self.floatingBins = self._readIn('floatingBins')
        else:
            self._buildFitWorkspace()

        # Make the card
        if 'card' not in self.recycle and not self.recycleAll:
            self._makeCard()

        # Do a prerun where we fit just this pass-fail pair and set the rpf to result
        if self.prerun and not self.recycleAll:
            print 'Pre-running '+self.tag+' '+self.name+' to get a better estimate of the transfer function'
            self.workspace.writeToFile(self.projPath+'base_'+self.name+'.root',True)  
            runMLFit([self],'0','5','',skipPlots=True,prerun=True)    
            prerun_file = TFile.Open(self.projPath+'/fitDiagnosticsTest.root')
            if prerun_file:
                if prerun_file.GetListOfKeys().Contains('fit_b'):
                    prerun_result = prerun_file.Get('fit_b').floatParsFinal()
                elif prerun_file.GetListOfKeys().Contains('fit_s'):
                    prerun_result = prerun_file.Get('fit_s').floatParsFinal()
                else:
                    prerun_result = False
                if prerun_result != False:
                    for v in self.rpf.funcVars.keys():
                        prerun_coeff = prerun_result.find(v)
                        self.rpf.funcVars[v].setVal(prerun_coeff.getValV())
                        self.rpf.funcVars[v].setError(prerun_coeff.getValV()*0.5)
                        self.workspace.var(v).setVal(prerun_coeff.getValV())
                        self.workspace.var(v).setError(prerun_coeff.getValV()*0.5)

            else:
                raw_input('WARNING: Pre-run for '+self.tag+' '+self.name+'failed. Using original Rp/f parameters. Press any key to continue.')

        # Save out at the end
        if not self.recycleAll:
            self._saveOut()
            pickle.dump(self.pickleDict, open(self.projPath+'saveOut.p','wb'))

        # Very last thing - get a seg fault otherwise
        del self.workspace
        self.orgFile.Close()
        del self.pickleDict
        # del self.allVars

    # FUNCTIONS USED IN INITIALIZATION
    def _dummyTH2(self): # stores binning of space
        dummyTH2 = TH2F('dummyTH2','dummyTH2',len(self.fullXbins)-1,array.array('d',self.fullXbins),len(self.newYbins)-1,array.array('d',self.newYbins))
        return dummyTH2

# WRAPPER FUNCTIONS
def runMLFit(twoDs,rMin,rMax,systsToSet,skipPlots=False,prerun=False):
    # Set verbosity - chosen from first of configs
    verbose = ''
    if twoDs[0].verbosity != False:
        verbose = ' -v '+twoDs[0].verbosity
    
    # Set signal strength range - chosen from first of configs
    sig_option = ' --rMin '+rMin+' --rMax '+rMax

    # Set blinding (mask pass for each twoD that requires it)
    # For channel masking, need to send text2workspace arg to text2workspace.py via `--text2workspace "--channel-masks"`
    blind_option = ''
    blindedFit = False
    for twoD in twoDs:
        if twoD.blindedFit == True:
            blindedFit = True
            if blind_option == '':
                blind_option += ' mask_pass_SIG_'+twoD.name+'=1'
            else:
                blind_option += ',mask_pass_SIG_'+twoD.name+'=1'

    if blindedFit:
        blind_option = '--text2workspace "--channel-masks" --setParameters' + blind_option

    # Set card name and project directory
    card_name = 'card_'+twoDs[0].tag+'.txt' if not prerun else 'card_'+twoDs[0].name+'.txt'
    projDir = twoDs[0].tag if not prerun else twoDs[0].projPath

    # Determine if any nuisance/sysetmatic parameters should be set before fitting
    if systsToSet != '':
        if blind_option != '': blind_option = blind_option+','+systsToSet
        else: blind_option = '--setParameters '+systsToSet

    # Always set r to start at 1
    if blind_option != '': blind_option = blind_option+',r=1'
    else: blind_option = '--setParameters r=1'

    # Run Combine
    FitDiagnostics_command = 'combine -M FitDiagnostics -d '+card_name+' '+blind_option+' --saveWorkspace --cminDefaultMinimizerStrategy 0 ' + sig_option +verbose 

    with header.cd(projDir):
        command_saveout = open('FitDiagnostics_command.txt','w')
        command_saveout.write(FitDiagnostics_command)
        command_saveout.close()

        if os.path.isfile('fitDiagnosticsTest.root'):
            header.executeCmd('rm fitDiagnosticsTest.root')

        header.executeCmd(FitDiagnostics_command)

        if not os.path.isfile('fitDiagnosticsTest.root'):
            print "Combine failed and never made fitDiagnosticsTest.root. Quitting..."
            for i in twoDs:
                del i
            quit()

        diffnuis_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnosticsTest.root --abs -g nuisance_pulls.root'
        header.executeCmd(diffnuis_cmd)

        systematic_analyzer_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/systematicsAnalyzer.py '+card_name+' --all -f html > systematics_table.html'
        header.executeCmd(systematic_analyzer_cmd)

        # Make a PDF of the nuisance_pulls.root
        if os.path.exists('nuisance_pulls.root'):
            nuis_file = TFile.Open('nuisance_pulls.root')
            nuis_can = nuis_file.Get('nuisances')
            nuis_can.Print('nuisance_pulls.pdf','pdf')
            nuis_file.Close()

    # Save out Rp/f to a text file and make a re-run config
    for twoD in twoDs:
        for fittag in ['b','s']:
            param_out = open(twoD.projPath+'rpf_params_'+twoD.name+'fit'+fittag+'.txt','w')
            rerun_config = header.dictCopy(twoD.inputConfig)

            try:
                coeffs_final = TFile.Open(projDir+'/fitDiagnosticsTest.root').Get('fit_'+fittag).floatParsFinal()
                coeffIter_final = coeffs_final.createIterator()
                coeff_final = coeffIter_final.Next()
                while coeff_final:
                    if coeff_final.GetName() in twoD.rpf.funcVars.keys():
                        # Text file
                        param_out.write(coeff_final.GetName()+': ' + str(coeff_final.getValV()) + ' +/- ' + str(coeff_final.getError())+'\n')
                        # Re run config
                        for k in rerun_config['FIT'].keys():
                            if 'generic'+k in coeff_final.GetName():
                                rerun_config['FIT'][k]['ERROR'] = coeff_final.getError()
                                rerun_config['FIT'][k]['NOMINAL'] = coeff_final.getValV()
                                rerun_config['FIT'][k]['MIN'] = coeff_final.getValV()-3*coeff_final.getError()
                                rerun_config['FIT'][k]['MAX'] = coeff_final.getValV()+3*coeff_final.getError()

                    # Next
                    coeff_final = coeffIter_final.Next()
            except:
                print 'Unable to write fit_'+fittag+ ' parameters to text file'

            # Close text file
            param_out.close()
            # Write out dictionary
            rerun_out = open(twoD.projPath+'rerunConfig_fit'+fittag+'.json', 'w')
            json.dump(rerun_config,rerun_out,indent=2, sort_keys=True)
            rerun_out.close()

    if not skipPlots:
        with header.cd(projDir):
            bshapes_cmd = 'PostFit2DShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -o postfitshapes_b.root -f fitDiagnosticsTest.root:fit_b --postfit --sampling --samples 100 --print 2> PostFitShapes2D_stderr_b.txt'
            header.executeCmd(bshapes_cmd)
            sshapes_cmd = 'PostFit2DShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -o postfitshapes_s.root -f fitDiagnosticsTest.root:fit_s --postfit --sampling --samples 100 --print 2> PostFitShapes2D_stderr_s.txt'
            header.executeCmd(sshapes_cmd)

            covMtrx_File = TFile.Open('fitDiagnosticsTest.root')
            fit_result = covMtrx_File.Get("fit_b")
            if hasattr(fit_result,'correlationMatrix'):
                corrMtrx = header.reducedCorrMatrixHist(fit_result)
                corrMtrxCan = TCanvas('c','c',1400,1000)
                corrMtrxCan.cd()
                corrMtrxCan.SetBottomMargin(0.22)
                corrMtrxCan.SetLeftMargin(0.17)
                corrMtrxCan.SetTopMargin(0.06)

                corrMtrx.Draw('colz text')
                corrMtrxCan.Print('correlation_matrix.png','png')
            else:
                print 'WARNING: Not able to produce correlation matrix.'


