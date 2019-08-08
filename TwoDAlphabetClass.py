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
import subprocess
import cPickle as pickle
import os, sys, array, json
import pprint
pp = pprint.PrettyPrinter(indent = 2)

import header
import RpfHandler

import ROOT
from ROOT import *

gStyle.SetOptStat(0)

class TwoDAlphabet:
    # If you just want to do everything yourself
    def __init__(self):
        pass

    def __del__(self):
        del self.workspace

    # Initialization setup to just build workspace. All other steps must be called with methods
    def __init__ (self,jsonFileName,quicktag='',recycleAll=False,stringSwaps={}): # jsonFileNames is a list
        self.allVars = []    # This is a list of all RooFit objects made. It never gets used for anything but if the
                        # objects never get saved here, then the python memory management will throw them out
                        # because of conflicts with the RooFit memory management. It's a hack.

        # Setup config
        self.jsonFileName = jsonFileName
        self.stringSwaps = stringSwaps
        self.inputConfig = header.openJSON(self.jsonFileName)

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

        # Setup bool options
        self.freezeFail = self._getOption('freezeFail')
        self.blindedPlots = self._getOption('blindedPlots')
        self.blindedFit = self._getOption('blindedFit')
        self.plotUncerts = self._getOption('plotUncerts')
        self.draw = self._getOption('draw')
        self.verbosity = str(self._getOption('verbosity')) # will eventually need to be fed as a string so just convert now
        self.fitGuesses = self._getOption('fitGuesses')
        self.prerun = self._getOption('prerun')
        self.overwrite = self._getOption('overwrite')
        self.recycle = self._getOption('recycle')
        self.plotTogether = self._getOption('plotTogether')
        self.recycleAll = recycleAll

        # Setup a directory to save
        self.projPath = self._projPath()

        # Pickle reading if recycling
        if self.recycle != [] or self.recycleAll:
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
        self.fullXbins = list(self.newXbins['LOW']) # need list() to make a copy - not a reference
        for c in ['SIG','HIGH']:
            self.fullXbins.extend(self.newXbins[c][1:])

        print self.fullXbins

        # Make systematic uncertainty plots
        if self.plotUncerts and not self.recycleAll:
            self._makeSystematicPlots()

        # Run pseudo2D for fit guesses and make the config to actually run on
        if ("runConfig" not in self.recycle and not self.recycleAll):
            self._makeFitGuesses()

        # Initialize rpf class
        if 'organizedDict' not in self.recycle and not self.recycleAll:
        #     self.rpf = self._readIn('rpf')
        # else:
            self.rpf = RpfHandler.RpfHandler(self.inputConfig['FIT'],self.name,self._dummyTH2())

        # Organize everything for workspace building
        if 'organizedDict' in self.recycle or self.recycleAll:
            self.organizedDict = self._readIn('organizedDict')
            self.orgFile = TFile.Open(self.projPath+'organized_hists.root') # have to save out the histograms to keep them persistent past this function
        else:
            self.orgFile = TFile(self.projPath+'organized_hists.root','RECREATE') # have to save out the histograms to keep them persistent past this function
            self.organizedDict = {}
            self._inputOrganizer()

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
            runMLFit([self],'0','5',skipPlots=True)    
            prerun_file = TFile.Open(self.projPath+'/fitDiagnostics.root')
            if prerun_file:
                if prerun_file.GetListOfKeys().Contains('fit_b'):
                    prerun_result = prerun_file.Get('fit_b').floatParsFinal()
                elif prerun_file.GetListOfKeys().Contains('fit_s'):
                    prerun_result = prerun_file.Get('fit_s').floatParsFinal()
                else:
                    prerun_result = False
                if prerun_result != False:
                    for v in self.rpf.rpfVars.keys():
                        prerun_coeff = prerun_result.find(v)
                        self.rpf.rpfVars[v].setVal(prerun_coeff.getValV())
                        self.rpf.rpfVars[v].setError(prerun_coeff.getValV()*0.5)
                        # self.rpf.rpfVars[v].setMin(max(self.rpf.rpfVars[v].getMin(), prerun_coeff.getValV()-2*prerun_coeff.getError()))
                        # self.rpf.rpfVars[v].setMax(min(self.rpf.rpfVars[v].getMax(), prerun_coeff.getValV()+2*prerun_coeff.getError()))
                        self.workspace.var(v).setVal(prerun_coeff.getValV())
                        self.workspace.var(v).setError(prerun_coeff.getValV()*0.5)
                        # self.workspace.var(v).setMin(max(self.rpf.rpfVars[v].getMin(), prerun_coeff.getValV()-2*prerun_coeff.getError()))
                        # self.workspace.var(v).setMax(min(self.rpf.rpfVars[v].getMax(), prerun_coeff.getValV()+2*prerun_coeff.getError()))
            else:
                raw_input('WARNING: Pre-run for '+self.tag+' '+self.name+'failed. Using original Rp/f parameters. Press any key to continue.')

        # Save out at the end
        if not self.recycleAll:
            self._saveOut()
            pickle.dump(self.pickleDict, open(self.projPath+'saveOut.p','wb'))

        # Very last thing - get a seg fault otherwise
        del self.workspace

    # FUNCTIONS USED IN INITIALIZATION
    def _configGlobalVarReplacement(self):
        #####################################
        # Do GLOBAL variable substitution   #
        # --------------------------------- #
        # Relies on certain JSON structure. #
        # Anything marked as 'changeable'   #
        # needs to be checked for GLOBAL    #
        # variable.                         # 
        #                                   #
        # - HELP is unchangeable            #
        # - mainkeys are unchangeable       #
        # - subkeys are changeable in       #
        #   PROCESS, SYSTEMATIC, and FIT    #
        # - subsubkeys are changeable       #
        # - subsubkey values are changeable #
        #                                   #
        # CURRENTLY ONLY SUPPORTS STRING    #
        # REPLACEMENT                       #
        #####################################

        # Add possible string swaps to the 'GLOBAL' dict of the config
        for s in self.stringSwaps.keys():
            if s in self.inputConfig['GLOBAL'].keys():
                print 'ERROR: A command line string replacement (%s) conflicts with one already in the configuration file. Quitting...' %(s)
                quit()
            self.inputConfig['GLOBAL'][s] = self.stringSwaps[s]

        print "Doing GLOBAL variable replacement in input json...",
        for old_string in self.inputConfig['GLOBAL'].keys():
            new_string = self.inputConfig['GLOBAL'][old_string]
            
            if old_string != "HELP":                                            # For each key in GLOBAL that is not HELP
                for mainkey in self.inputConfig.keys():                         # For each main (top level) key in config that isn't GLOBAL
                    if mainkey != 'GLOBAL' and mainkey != 'OPTIONS':            # Mainkeys are all unchangeable (uppercase) so no check performed
                        for subkey in self.inputConfig[mainkey].keys():         # For each subkey of main key dictionary
                            if old_string in subkey:                            # Check subkey for old_string
                                self.inputConfig[mainkey][subkey.replace(old_string,new_string)] = self.inputConfig[mainkey].pop(subkey)  # replace it
                                subkey = subkey.replace(old_string,new_string)

                            # If the subkey value is not a dict, then check one more time
                            if type(self.inputConfig[mainkey][subkey]) != dict:
                                if old_string in self.inputConfig[mainkey][subkey]:
                                    self.inputConfig[mainkey][subkey] = self.inputConfig[mainkey][subkey].replace(old_string,new_string)
                            # If it is a dict, go deeper
                            else:
                                for subsubkey in self.inputConfig[mainkey][subkey].keys():  # so loop though subsubkeys
                                    if old_string in subsubkey:                                   # check subsubkey
                                        self.inputConfig[mainkey][subkey][subsubkey.replace(old_string,new_string)] = self.inputConfig[mainkey][subkey].pop(subsubkey)   # replace it
                                        subsubkey = subsubkey.replace(old_string,new_string)

                                    if type(self.inputConfig[mainkey][subkey][subsubkey]) == str:
                                        if old_string in self.inputConfig[mainkey][subkey][subsubkey]:                               # check subsubkey val
                                            self.inputConfig[mainkey][subkey][subsubkey] = self.inputConfig[mainkey][subkey][subsubkey].replace(old_string,new_string) # replace it

    def _dummyTH2(self): # stores binning of space
        dummyTH2 = TH2F('dummyTH2','dummyTH2',len(self.fullXbins)-1,array.array('d',self.fullXbins),len(self.newYbins)-1,array.array('d',self.newYbins))
        return dummyTH2

    def _getRRVs(self):
        xRRVs = {}
        # Y
        yname = self.yVarName+'_'+self.name
        ylow = self.newYbins[0]
        yhigh = self.newYbins[-1]
        yRRV = RooRealVar(yname,yname,ylow,yhigh)
        yBinArray = array.array('d',self.newYbins)
        yRooBinning = RooBinning(len(self.newYbins)-1,yBinArray)
        yRRV.setBinning(yRooBinning)

        # X
        for c in ['LOW','SIG','HIGH']:
            xname = self.xVarName+'_'+c+'_'+self.name
            xlow = self.newXbins[c][0]
            xhigh = self.newXbins[c][-1]
            xRRVs[c] = RooRealVar(xname,xname,xlow,xhigh)
            xBinArray = array.array('d',self.newXbins[c])
            xRooBinning = RooBinning(len(self.newXbins[c])-1,xBinArray)
            xRRVs[c].setBinning(xRooBinning)

        return xRRVs,yRRV

    def _projPath(self):
        if self.tag != '':
            if not os.path.isdir(self.tag+'/'):
                print 'Making dir '+self.tag+'/'
                os.mkdir(self.tag+'/')
            elif self.overwrite:
                subprocess.call(['rm -rf '+self.tag],shell=True)
                print 'Making dir '+self.tag+'/'
                os.mkdir(self.tag+'/')

            if not os.path.isdir(self.tag+'/'+self.name): os.mkdir(self.tag+'/'+self.name)
            if not os.path.isdir(self.tag+'/'+self.name+'/plots/'): os.mkdir(self.tag+'/'+self.name+'/plots/')
            if not os.path.isdir(self.tag+'/'+self.name+'/plots/fit_b/'): os.mkdir(self.tag+'/'+self.name+'/plots/fit_b/')
            if not os.path.isdir(self.tag+'/'+self.name+'/plots/fit_s/'): os.mkdir(self.tag+'/'+self.name+'/plots/fit_s/')
            if self.plotUncerts and not os.path.isdir(self.tag+'/'+self.name+'/UncertPlots/'): os.mkdir(self.tag+'/'+self.name+'/UncertPlots/')

            dirname = self.tag+'/'+self.name+'/'
        
        else:
            if not os.path.isdir(self.name+'/'): os.mkdir(self.name+'/')
            if not os.path.isdir(self.name+'/plots/'): os.mkdir(self.name+'/plots/')
            if not os.path.isdir(self.name+'/plots/fit_b/'): os.mkdir(self.name+'/plots/fit_b/')
            if not os.path.isdir(self.name+'/plots/fit_s/'): os.mkdir(self.name+'/plots/fit_s/')
            if self.plotUncerts and not os.path.isdir(self.name+'/UncertPlots/'): os.mkdir(self.name+'/UncertPlots/')

            elif self.overwrite:
                subprocess.call(['rm -rf '+self.name],shell=True)
                os.mkdir(self.name+'/')
                os.mkdir(self.name+'/plots/')
                os.mkdir(self.name+'/plots/fit_b/')
                os.mkdir(self.name+'/plots/fit_s/')

            dirname = self.name+'/'

        return dirname

    def _getBinning(self, binDict):
        # DOCUMENT

        # If running a blinded fit, then we want to do a combined fit over 
        # two categories: below and above the signal region. This requires
        # generating histograms in those three regions and it's useful
        # to have different binning for all of those. If the signal region
        # is not blinded then we can fit the entire region but it's convenient
        # to still do three categories for later parts of the code. So here are
        # the options.
        # 1) It may be desired or even required to bin the fit in three categories
        # each with its own binning structure (say statistics are good in 
        # region below the signal region but bad above it so you'd like to
        # use more and fewer bins, respectively). 
        # 2) Additionally, variable binning can be used for each category. 
        # 3) Single binning strategy across all three regions and only defined
        # once in the configuration.
        # The only requirement for any of this is that the bin walls of the new
        # binning match a bin wall of the input histograms (you can't make bins smaller or split bins!)

        # For config input, this breaks down into
        # Standard bins over one category - one NBINS,MIN,MAX
        # Standard bins over three categories - three NBINS,MIN,MAX (organized by dict)
        # Custom bins over one category - list of bin walls
        # Custom bins over three cateogories - three lists of bin walls (organized by dict)

        # Finally, we need to get the bin normalizations correct. So we look at one of the
        # input histograms to get the bin widths and use that as the base to normalize the 
        # rebinning to.

        temp_input_file = TFile.Open(self.inputConfig['PROCESS']['data_obs']['FILE'])
        temp_input_hist = temp_input_file.Get(self.inputConfig['PROCESS']['data_obs']['HISTFAIL'])
        oldXwidth = (temp_input_hist.GetXaxis().GetXmax() - temp_input_hist.GetXaxis().GetXmin())/temp_input_hist.GetNbinsX()
        oldYwidth = (temp_input_hist.GetYaxis().GetXmax() - temp_input_hist.GetYaxis().GetXmin())/temp_input_hist.GetNbinsY()

        for v in ['X','Y']:
            # ONE CATEGORY - VARIABLE
            if 'BINS' in binDict[v].keys():
                # If X, take one list of bins, split around signal region into three lists, feed back as dictionary
                if v == 'X':
                    new_bins = header.splitBins(binDict[v]['BINS'],self.sigStart,self.sigEnd)
                else: 
                    new_bins = binDict[v]['BINS']

            # ONE CATEGORY - CONSTANT
            elif ('MIN' in binDict[v].keys()) and ('MAX' in binDict[v].keys()) and ('NBINS' in binDict[v].keys()):
                new_min = binDict[v]['MIN']
                new_max = binDict[v]['MAX']
                new_nbins = binDict[v]['NBINS']
                new_width = float(new_max-new_min)/float(new_nbins)

                bin_walls = []
                for i in range(new_nbins):
                    b = new_min + new_width*i
                    bin_walls.append(b)
                bin_walls.append(new_max)

                if v == 'X':
                    new_bins = header.splitBins(bin_walls,self.sigStart,self.sigEnd)
                else:
                    new_bins = bin_walls

            # THREE CATEGORIES but only if in X
            elif v == 'X':
                if ('LOW' in binDict[v].keys()) and ('SIG' in binDict[v].keys()) and ('HIGH' in binDict[v].keys()):
                    new_bins = {}
                    for c in ['LOW','SIG','HIGH']:
                        # Check each category if variable or not
                        if 'BINS' in binDict[v][c].keys():
                            new_bins[c] = binDict[v][c]['BINS']
                        # or constant in each category
                        elif ('MIN' in binDict[v][c].keys()) and ('MAX' in binDict[v][c].keys()) and ('NBINS' in binDict[v][c].keys()):
                            new_min = binDict[v][c]['MIN']
                            new_max = binDict[v][c]['MAX']
                            new_nbins = binDict[v][c]['NBINS']
                            new_width = float(new_max-new_min)/float(new_nbins)

                            new_bins[c] = []
                            b = new_min
                            while (new_max - b) > new_width:
                                new_bins[c].append(b)
                                b += new_width
                            new_bins[c].append(new_max)

            else:
                print 'No user bins specified. Will use binning of input histograms (data_obs_pass).'
                temp_file = TFile.Open(self.inputConfig['PROCESS']['data_obs']['FILE'])
                temp_TH2 = temp_file.Get(self.inputConfig['PROCESS']['data_obs']['HISTPASS'])
                new_bins = []
                if v == 'X':
                    for b in range(1,temp_TH2.GetNbinsX()+1):
                        new_bins.append(temp_TH2.GetXaxis().GetBinLowEdge(b))
                    new_bins.append(temp_TH2.GetXaxis().GetXmax())

                elif v == 'Y':
                    for b in range(1,temp_TH2.GetNbinsY()+1):
                        new_bins.append(temp_TH2.GetYaxis().GetBinLowEdge(b))
                    new_bins.append(temp_TH2.GetYaxis().GetXmax())

            if v == 'X':
                newXbins = new_bins
        
            elif v == 'Y':
                newYbins = new_bins

        return newXbins,newYbins,oldXwidth,oldYwidth

    def _getFullXbin(self,xbin,c):
        # Evaluate for the bin - a bit tricky since it was built with separate categories
        # Determine the category and x bin from that
        # Ex.
        # full_x_bins = [a,b,c,d,e,f,g]; newXbins[cat] = [c,d,e]
        # newXbins[cat][xbin] = upper wall of xbin
        # full_x_bins.index(newXbins[cat][xbin]) = index of upper global wall OR index of bin win want in the the histogram (remember those bin indices start at 1 not 0)

        return self.fullXbins.index(self.newXbins[c][xbin])

    def _getOption(self,optionName):
        # If it's in the config, just set it
        if optionName in self.inputConfig['OPTIONS'].keys():
            option_return = self.inputConfig['OPTIONS'][optionName]
        # Default to true
        elif optionName in ['blindedPlots','blindedFit']:
            print 'WARNING: '+optionName+' boolean not set explicitely. Default to True.'
            option_return = True
        # Default to false
        elif optionName in ['freezeFail','fitGuesses','plotUncerts','prerun']:
            print 'WARNING: '+optionName+' boolean not set explicitely. Default to False.'
            option_return = False
        elif optionName == 'verbosity':
            option_return = 0
        else:
            if optionName == 'recycle':
                print 'WARNING: '+optionName+' boolean not set explicitely. Default to [].'
                option_return = []
            else:
                print 'WARNING: '+optionName+' boolean not set explicitely. Default to False.'
                option_return = False
            

        return option_return

    def _saveOut(self):
        # runConfig
        file_out = open(self.projPath+'runConfig.json', 'w')
        json.dump(self.inputConfig,file_out,indent=2, sort_keys=True)
        file_out.close()

        self.pickleDict['name'] = self.name
        self.pickleDict['tag'] = self.tag
        self.pickleDict['xVarName'] = self.xVarName
        self.pickleDict['yVarName'] = self.yVarName
        self.pickleDict['xVarTitle'] = self.xVarTitle
        self.pickleDict['yVarTitle'] = self.yVarTitle
        self.pickleDict['sigStart'] = self.sigStart
        self.pickleDict['sigEnd'] = self.sigEnd
        self.pickleDict['freezeFail'] = self.freezeFail
        self.pickleDict['blindedFit'] = self.blindedFit
        self.pickleDict['plotTogether'] = self.plotTogether

        # Setup a directory to save
        self.projPath = self._projPath()

        # bins
        self.pickleDict['newXbins'] = self.newXbins
        self.pickleDict['full_x_bins'] = self.fullXbins
        self.pickleDict['newYbins'] = self.newYbins

        # rpf - Don't do this - takes up +5 GB
        # self.pickleDict['rpf'] = self.rpf
        self.pickleDict['rpfVarNames'] = self.rpf.getRpfVarNames()

        # organizedDict
        self.pickleDict['organizedDict'] = self.organizedDict

        # floatingBins
        self.pickleDict['floatingBins'] = self.floatingBins

        # workspace
        self.workspace.writeToFile(self.projPath+'base_'+self.name+'.root',True)  

    def _readIn(self,attrname):
        if attrname == 'runConfig':
            return header.openJSON(self.projPath+'runConfig.json')

        elif attrname == 'newXbins': 
            return self.pickleFile['newXbins']

        elif attrname == 'newYbins':
            return self.pickleFile['newYbins']

        # elif attrname == 'rpf': 
        #     return self.pickleFile['rpf']

        elif attrname == 'organizedDict':
            return self.pickleFile['organizedDict']

        elif attrname == 'floatingBins':
            return self.pickleFile['floatingBins']

        elif attrname == 'workspace':
            return TFile.Open(self.projPath+'base_'+self.name+'.root').Get('w_'+self.name)

    def _makeSystematicPlots(self):
        for proc in self.inputConfig['PROCESS'].keys():
            if proc != 'HELP':
                for syst in self.inputConfig['PROCESS'][proc]['SYSTEMATICS']:
                    # For each systematic in each process

                    print proc + '_'+syst

                    # If code 2 or 3 (shape based), grab the up and down shapes for passing and project onto the Y axis
                    if self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 2:
                        thisFile = TFile.Open(self.inputConfig['PROCESS'][proc]['FILE'])
                        passUp2D = thisFile.Get(self.inputConfig['SYSTEMATIC'][syst]['HISTPASS_UP'])
                        passDown2D = thisFile.Get(self.inputConfig['SYSTEMATIC'][syst]['HISTPASS_DOWN'])
                        passUp = passUp2D.ProjectionY(proc + '_'+syst,21,28)
                        passDown = passDown2D.ProjectionY(proc + '_'+syst+'_down',21,28)

                    elif self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 3:
                        if 'FILE_UP_'+proc in self.inputConfig['SYSTEMATIC'][syst]:
                            thisFileUp = TFile.Open(self.inputConfig['SYSTEMATIC'][syst]['FILE_UP_'+proc])
                            thisFileDown = TFile.Open(self.inputConfig['SYSTEMATIC'][syst]['FILE_DOWN_'+proc])

                        elif 'FILE_UP_*' in self.inputConfig['SYSTEMATIC'][syst]:
                            thisFileUp = TFile.Open(self.inputConfig['SYSTEMATIC'][syst]['FILE_UP_*'].replace('*',proc))
                            thisFileDown = TFile.Open(self.inputConfig['SYSTEMATIC'][syst]['FILE_DOWN_*'].replace('*',proc))

                        else:
                            print 'Could not identify file for ' + proc +', '+syst

                        passUp = thisFileUp.Get(self.inputConfig['SYSTEMATIC'][syst]['HISTPASS']).ProjectionY(proc + '_'+syst,21,28)
                        passDown = thisFileDown.Get(self.inputConfig['SYSTEMATIC'][syst]['HISTPASS']).ProjectionY(proc + '_'+syst+'_down',21,28)

                    else:
                        continue

                    # Setup the nominal shape (consistent across all of the pltos)
                    fileNom = TFile.Open(self.inputConfig['PROCESS'][proc]['FILE'])
                    passNom = fileNom.Get(self.inputConfig['PROCESS'][proc]['HISTPASS']).ProjectionY(proc + '_'+syst+'_nom',21,28)

                    thisCan = TCanvas('canvas_'+proc+'_'+syst,'canvas_'+proc+'_'+syst,700,600)
                    thisCan.cd()
                    passNom.SetLineColor(kBlack)
                    passNom.SetFillColor(kYellow-9)
                    passUp.SetLineColor(kRed)
                    passDown.SetLineColor(kBlue)

                    passUp.SetLineStyle(9)
                    passDown.SetLineStyle(9)
                    passUp.SetLineWidth(2)
                    passDown.SetLineWidth(2)

                    histList = [passNom,passUp,passDown]

                    # Set the max of the range so we can see all three histograms on the same plot
                    yMax = histList[0].GetMaximum()
                    maxHist = histList[0]
                    for h in range(1,len(histList)):
                        if histList[h].GetMaximum() > yMax:
                            yMax = histList[h].GetMaximum()
                            maxHist = histList[h]
                    for h in histList:
                        h.SetMaximum(yMax*1.1)

                    passNom.SetXTitle(self.inputConfig['BINNING']['Y']['TITLE'])
                    # passUp.SetXTitle(inputConfig['BINNING']['Y']['TITLE'])
                    # passDown.SetXTitle(inputConfig['BINNING']['Y']['TITLE'])
                    passNom.SetTitle(proc + ' - ' + syst + ' uncertainty')
                    passNom.SetTitleOffset(1.2,"X")

                    passNom.Draw('hist')
                    passUp.Draw('same hist')
                    passDown.Draw('same hist')
                    

                    thisCan.Print(self.projPath+'/UncertPlots/Uncertainty_'+proc+'_'+syst+'.pdf','pdf')

    def _makeFitGuesses(self,nslices=6,sigma=5):
        # Grab everything
        infile = TFile.Open(self.inputConfig['PROCESS']['data_obs']['FILE'])

        # Initiliaze QCD estimate
        data_pass = infile.Get(self.inputConfig['PROCESS']['data_obs']['HISTPASS'])
        data_pass = data_pass.Clone('qcd_pass')
        data_fail = infile.Get(self.inputConfig['PROCESS']['data_obs']['HISTFAIL'])
        data_fail = data_fail.Clone('qcd_fail')

        # Subtract away non-qcd processes
        for nonqcd in [process for process in self.inputConfig['PROCESS'].keys() if process != 'HELP']:
            if self.inputConfig['PROCESS'][nonqcd]['CODE'] > 1: # remove data and signal from considerations
                print 'Subtracting ' + nonqcd
                nonqcd_file = TFile.Open(self.inputConfig['PROCESS'][nonqcd]['FILE'])
                nonqcd_pass = nonqcd_file.Get(self.inputConfig['PROCESS'][nonqcd]['HISTPASS'])
                nonqcd_fail = nonqcd_file.Get(self.inputConfig['PROCESS'][nonqcd]['HISTFAIL'])
                data_pass.Add(nonqcd_pass,-1)
                data_fail.Add(nonqcd_fail,-1)

        # Zero any negative bins
        for reghist in [data_pass,data_fail]:
            for xbin in range(1,data_pass.GetNbinsX()+1):
                for ybin in range(1,data_pass.GetNbinsY()+1):
                    if reghist.GetBinContent(xbin,ybin) < 0:
                        reghist.SetBinContent(xbin,ybin,0)
                    if reghist.GetBinContent(xbin,ybin) < 0:
                        reghist.SetBinContent(xbin,ybin,0)

        ###########
        # Binning #
        ###########
        if self.fitGuesses == True: # Use binning defined in 'BINNING' section
            xbins_pseudo = self.newXbins['LOW']
            xbins_pseudo.extend(self.newXbins['SIG'][1:])
            xbins_pseudo.extend(self.newXbins['HIGH'][1:])
            ybins_pseudo = self.newYbins

        elif self.fitGuesses == 'auto': # Use x binning from 'BINNING' section and do auto binning for y
            xbins_pseudo = self.newXbins['LOW']
            xbins_pseudo.extend(self.newXbins['SIG'][1:])
            xbins_pseudo.extend(self.newXbins['HIGH'][1:])
            ybins_pseudo = 'auto'

        elif type(self.fitGuesses) == dict: # Use separate binning scheme provided by user
            xbins_pseudo_dict, ybins_pseudo, oldXwidth_pseudo, oldYwidth_pseudo = self._getBinning(self.fitGuesses)
            xbins_pseudo = list(xbins_pseudo_dict['LOW']) # need list() to make a copy - not a reference
        
            for c in ['SIG','HIGH']:
                xbins_pseudo.extend(xbins_pseudo_dict[c][1:])

        else:
            xbins_pseudo = self.fullXbins
            ybins_pseudo = self.newYbins

        print xbins_pseudo
        print ybins_pseudo
        ##############################
        #   Auto derive new y bins   #
        ##############################
        if self.fitGuesses == 'auto' and nslices > 0:
            # Need to figure out how to slice up the y-axis based on the total statistics (pass+fail) for each 2D bin
            data_total = data_pass.Clone('data_total')     
            data_total.Add(data_fail)

            # Get an average number of events per slice
            total_events = data_total.Integral()
            events_per_slice = total_events/float(nslices)

            # Loop over ybins, get the number of events in that ybin, and check if less than events_per_slice
            ysum = 0            # ysum for after we've confirmed that the next bin does not lead to a number greater than events_per_slice (for printing only)
            ysum_temp = 0       # ysum allowed to overflow (for counting)
            new_y_bins = []
            for ybin in range(1,len(self.newYbins)):                 # For each ybin
                # If on last ybin, add the final edge to the list
                # - - - - - - - - - - - - - - - - - - - - - - - - 
                # There's an interesting note to make here about the number of slices that one can do. There are three important facts here:
                # 1) this method only reduces the number of bins
                # 2) every bin edge in the new binning was a bin edge in the previous binning
                # 3) new bins have fewer than the target events_per_slice UNLESS the new bin is only made up of one bin in which case it has more than events_per_slice
                #
                # This means two scenarios happen
                # Example 1) if old bins 1-3 have fewer events than events_per_slice but old bins 1-4 have more than events_per_slice,
                #            then new bin 1 will contain the contents of old bins 1-3 (no 4). 
                # Example 2) But if the content of bin 4 is GREATER THAN events_per_slice,
                #            then new bin 2 will contain the contents of ONLY old bin 4. 
                #
                # Here is the interesting part. If the Example 2 scenario happens enough times (as is possible with a peaked distribution), then 
                # there won't be enough events to "go around" for the tail of the distribution. So you might ask for 9 slices but because 
                # slices 2-4 have many more events in them than 3*events_per_slice, then the next slices will have to count more events in
                # the tail and you'll "use" all of the events in the tail before you get to the 9th slice. In other words, this method naturally caps itself at a certain number of slices.
                
                if ybin == data_total.GetNbinsY() or len(new_y_bins) == nslices: # If final bin or max number of slices reached
                    new_y_bins.append(self.newYbins[-1])  # Add the final edge
                    break                                                   # This final bin will most likely have more than events_per_slice in it

                # Otherwise, if we're still building the slice list...
                else:
                    for xbin in range(1,len(xbins_pseudo)):             # For each xbin
                        ysum_temp += data_total.GetBinContent(xbin,ybin)     # Add the 2D bin content to the temp sum for the ybin

                    # If less, set ysum and go onto the next bin
                    if ysum_temp < events_per_slice:
                        ysum = ysum_temp
                        
                    # Otherwise, cut off the slice at the previous ybin and restart the sum with this ybin
                    else:
                        ysum_temp = 0
                        for xbin in range(1,len(xbins_pseudo)):
                            ysum_temp += data_total.GetBinContent(xbin,ybin)
                        new_y_bins.append(data_total.GetYaxis().GetBinLowEdge(ybin))

            ybins_pseudo = new_y_bins

        print 'Will bin y-axis for fit guesses using bins ',
        print ybins_pseudo

        # Rebin x and y
        rebinned_x_pass = header.copyHistWithNewXbins(data_pass,xbins_pseudo,'rebinned_x_pass')#,self.oldXwidth)
        rebinned_x_fail = header.copyHistWithNewXbins(data_fail,xbins_pseudo,'rebinned_x_fail')#,self.oldXwidth)

        rebinned_pass = header.copyHistWithNewYbins(rebinned_x_pass,ybins_pseudo,'rebinned_pass')#,self.oldYwidth)
        rebinned_fail = header.copyHistWithNewYbins(rebinned_x_fail,ybins_pseudo,'rebinned_fail')#,self.oldYwidth)

        ######################################################
        #   Rebin the distributions according to new y bins  #
        ######################################################

        # Blind if necessary
        if self.blindedFit:
            final_pass = header.makeBlindedHist(rebinned_pass,[self.sigStart,self.sigEnd])
            final_fail = header.makeBlindedHist(rebinned_fail,[self.sigStart,self.sigEnd])

        # Otherwise just get the Rpf
        else:
            final_pass = rebinned_pass
            final_fail = rebinned_fail

        final_pass.SetName('final_pass')
        final_pass.SetTitle('final_pass')
        final_fail.SetName('final_fail')
        final_fail.SetTitle('final_fail')

        final_pass.Sumw2()
        final_fail.Sumw2()
        RpfToRemap = final_pass.Clone('RpfToRemap')
        RpfToRemap.Divide(final_fail)

        Rpf = header.remapToUnity(RpfToRemap)

        # Plot comparisons out
        header.makeCan('prefit_pass_fail',self.projPath,[final_pass,final_fail],xtitle=self.xVarName,ytitle=self.yVarName)
        header.makeCan('prefit_rpf_lego',self.projPath,[RpfToRemap,Rpf],xtitle=self.xVarName,ytitle=self.yVarName)

        ###############################################
        # Determine fit function from the inputConfig #
        ###############################################
        if self.fitGuesses != False:
            if 'XPFORM' in self.inputConfig['FIT'].keys() and 'YPFORM' in self.inputConfig['FIT'].keys():
                # Do some quick checks to make sure these are formatted correctly
                header.checkFitForm(self.inputConfig['FIT']['XPFORM'],self.inputConfig['FIT']['YPFORM'])
                # Determine number of params in each direction
                nxparams = max([int(param[1:]) for param in self.inputConfig['FIT'].keys() if param.find('X') != -1 and param != 'XPFORM'])
                nyparams = max([int(param[1:]) for param in self.inputConfig['FIT'].keys() if param.find('Y') != -1 and param != 'YPFORM'])
                # Get the strings
                xFuncString = header.RFVform2TF1(self.inputConfig['FIT']['XPFORM'],-1)
                yFuncString = header.RFVform2TF1(self.inputConfig['FIT']['YPFORM'],-1)
                yFuncString = yFuncString.replace('y','x')

                # Make the fixed form of the formula
                funcString = ''
                paramIndex = 0
                for xparam in range(nxparams):
                    for yparam in range(nyparams):
                        funcString += '['+str(paramIndex)+']*x**'+str(xparam)+'*y**'+str(yparam)+'+'
                        paramIndex += 1
                funcString = funcString[:-1]


            elif 'PFORM' in self.inputConfig['FIT'].keys():
                funcString = self.inputConfig['FIT']['PFORM']
                # Reconstruct x
                xFuncString = ''
                nxparams = 0
                for xparam in self.inputConfig['FIT'].keys():
                    if 'X' in xparam:                                           # For each X*Y0
                        if 'Y0' in xparam:
                            powerIndex = xparam[xparam.find('X')+1]
                            xFuncString += '['+str(powerIndex)+']*x**'+str(powerIndex)+'+'
                            nxparams += 1
                xFuncString = xFuncString[:-1]

                # Reconstruct y
                yFuncString = ''
                nyparams = 0
                for yparam in self.inputConfig['FIT'].keys():
                    if 'X0Y' in yparam:                                           # For each X0Y*
                        powerIndex = yparam[yparam.find('Y')+1]
                        yFuncString += '['+str(powerIndex)+']*x**'+str(powerIndex)+'+'
                        nyparams += 1
                yFuncString = yFuncString[:-1]


            elif 'FORM' in self.inputConfig['FIT'].keys():
                # Need to take a 2D function and freeze it in one dimension for each y slice
                funcString = header.RFVform2TF1(self.inputConfig['FIT']['FORM'],0)
                yFuncString = '[0]' # since all y dependence will be absorbed by the funcString, there should be no y dependence on each of the parameters and thus we should fit each with a constant
                nxparams = max([int(param) for param in self.inputConfig['FIT'].keys() if param != 'FORM' and param != 'HELP']) +1

            elif 'XFORM' in self.inputConfig['FIT'].keys() and 'YFORM' in self.inputConfig['FIT'].keys():
                xFuncString = self.inputConfig['FIT']['XFORM']
                yFuncString = self.inputConfig['FIT']['YFORM']
                funcString = '('+header.RFVform2TF1(self.inputConfig['FIT']['XFORM'],0)+')*('+header.RFVform2TF1(self.inputConfig['FIT']['YFORM'],0).replace('x','y')+')'
                nxparams = max([int(param[1:]) for param in self.inputConfig['FIT'].keys() if param.find('X') != -1 and param != 'XFORM'])
                nyparams = max([int(param[1:]) for param in self.inputConfig['FIT'].keys() if param.find('Y') != -1 and param != 'YFORM'])

            else:
                print 'Fit form not supported in get_fit_guesses.py. Quitting...'
                quit()

            pseudo2D_results = TFile.Open(self.projPath+'/pseudo2D_results.root','RECREATE')
            pseudo2D_results.cd()

            if 'FORM' in self.inputConfig['FIT'].keys():# or ('XFORM' in self.inputConfig['FIT'].keys() and 'YFORM' in self.inputConfig['FIT'].keys()):
                pseudo2D_Rpf = TF2('pseudo2D_Rpf',funcString,0,1,0,1)

                for p in range(nxparams):
                    pseudo2D_Rpf.SetParameter(p,self.inputConfig['FIT'][str(p)]['NOMINAL'])
                    pseudo2D_Rpf.SetParLimits(p,self.inputConfig['FIT'][str(p)]['MIN'],self.inputConfig['FIT'][str(p)]['MAX'])

                Rpf.Fit(pseudo2D_Rpf)
                
                print 'Resetting fit parameters in input config'
                self.inputConfig['FIT']['FORM'] = funcString
                for ix in range(nxparams):
                    # for iy in range(nyparams):
                        param = str(ix)
                        pseudo2D_Rpf.SetParameter(ix,pseudo2D_Rpf.GetParameter(ix))
                        pseudo2D_Rpf.SetParError(ix,pseudo2D_Rpf.GetParError(ix))
                        if self.fitGuesses != False:
                            # self.inputConfig['FIT'][param] = {'NOMINAL':None,'MIN':None,'MAX':None}
                            self.inputConfig['FIT'][param]['NOMINAL'] = pseudo2D_Rpf.GetParameter(ix)
                            self.inputConfig['FIT'][param]['MAX'] = min(self.inputConfig['FIT'][param]['MAX'],pseudo2D_Rpf.GetParameter(ix)+sigma*pseudo2D_Rpf.GetParError(ix))
                            self.inputConfig['FIT'][param]['MIN'] = max(self.inputConfig['FIT'][param]['MIN'],pseudo2D_Rpf.GetParameter(ix)-sigma*pseudo2D_Rpf.GetParError(ix))
                            self.inputConfig['FIT'][param]['ERROR'] = pseudo2D_Rpf.GetParError(ix)

                pp.pprint(self.inputConfig['FIT'])
                            
            else:
                ##################################
                # Now do the fit in the y slices #
                ##################################
                fitResults = {}

                # Book TGraphs to store the fit results as a function of y bins
                unitYbins = array.array('d',[Rpf.GetYaxis().GetBinLowEdge(b) for b in range(1,Rpf.GetNbinsY()+1)]+[1])
                for xparam in range(nxparams):
                    fitResults['xparam_'+str(xparam)+'_vs_y'] = TH1F('xparam_'+str(xparam)+'_vs_y','xparam_'+str(xparam)+'_vs_y',Rpf.GetNbinsY(),unitYbins)

                # Project each y-axis bin and fit along x - save out coefficients to booked tgraph
                projXs = []
                for ybin in range(1,Rpf.GetNbinsY()+1):
                    # If doing FORM, define xFuncString here with y bin center plugged in
                    if 'FORM' in self.inputConfig['FIT'].keys():
                        thisYBinCenter = Rpf.GetYaxis().GetBinCenter(ybin)
                        xFuncString = funcString.replace('y',str(thisYBinCenter))

                    fitResults['fitSlice_'+str(ybin)] = TF1('fitSlice_'+str(ybin),xFuncString,0,1)

                    if 'FORM' in self.inputConfig['FIT'].keys():
                        for p in range(nxparams):
                            fitResults['fitSlice_'+str(ybin)].SetParameter(p,self.inputConfig['FIT'][str(p)]['NOMINAL'])
                            fitResults['fitSlice_'+str(ybin)].SetParLimits(p,self.inputConfig['FIT'][str(p)]['MIN'],self.inputConfig['FIT'][str(p)]['MAX'])

                    projX = Rpf.ProjectionX('rebinnedRpf_sliceX_'+str(ybin),ybin,ybin,'e o')
                    # projX.Draw('p e')
                    projX.Write()
                    projX.Fit(fitResults['fitSlice_'+str(ybin)],'EM')
                    
                    # projX.Draw('p e')
                    projX.SetMaximum(1.1)
                    projX.SetMinimum(0.0)
                    projX.SetTitle('fitSlice_'+str(ybin))

                    projXs.append(projX)

                    for ix in range(nxparams):
                        fitResults['xparam_'+str(ix)+'_vs_y'].SetBinContent(ybin,fitResults['fitSlice_'+str(ybin)].GetParameter(ix))
                        fitResults['xparam_'+str(ix)+'_vs_y'].SetBinError(ybin,fitResults['fitSlice_'+str(ybin)].GetParError(ix))
             
                if len(projXs) <= 6:
                    header.makeCan('fitSlices_0-6',self.projPath,projXs,xtitle=self.xVarName)

                else:
                    chunkedProjX = [projXs[i:i + 6] for i in xrange(0, len(projXs), 6)]
                    for i,chunk in enumerate(chunkedProjX):
                        header.makeCan('fitSlices_'+str(i*6)+'-'+str(6*(i+1)),self.projPath,chunk,xtitle=self.xVarName)

                ########################################################
                # And now fit these parameters as a function of y axis #
                ########################################################

                # Build fit for each parameter distribution along y-axis
                drawList = []
                for xparam in range(nxparams):
                    fitResults['fitParam_'+str(xparam)] = TF1('yFunc_'+str(xparam),yFuncString,0,1)
                    # Do the fit
                    fitResults['xparam_'+str(xparam)+'_vs_y'].Fit(fitResults['fitParam_'+str(xparam)],"EM")
                    drawList.append(fitResults['xparam_'+str(xparam)+'_vs_y'])
                    # Get and store parameters found
                    for iy in range(fitResults['fitParam_'+str(xparam)].GetNpar()):
                        fitResults['X'+str(xparam)+'Y'+str(iy)] = fitResults['fitParam_'+str(xparam)].GetParameter(iy)
                        fitResults['X'+str(xparam)+'Y'+str(iy)+'err'] = fitResults['fitParam_'+str(xparam)].GetParError(iy)

                if len(drawList) <= 6:
                    header.makeCan('xparam_v_y',self.projPath,drawList,xtitle=self.yVarName)
                else:
                    chunkedDrawList = [drawList[i:i + 6] for i in xrange(0, len(drawList), 6)]
                    for i,chunk in enumerate(chunkedDrawList):
                        header.makeCan('xparam_v_y_'+str(i),self.projPath,chunk,xtitle=self.yVarName)


                # Remove old fit values and store new ones in inputConfig if PFORM else just make the pseudo2D_Rpf for plotting
                pseudo2D_Rpf = TF2('pseudo2D_Rpf',funcString,0,1,0,1)
                paramIndex = 0

                if 'FORM' in self.inputConfig['FIT'].keys():
                    # self.inputConfig['FIT'] = {'FORM':funcString}
                    for p in range(nxparams):
                        param = 'X'+str(p)+'Y0'
                        pseudo2D_Rpf.SetParameter(paramIndex,fitResults[param])
                        pseudo2D_Rpf.SetParError(paramIndex,fitResults[param+'err'])
                        # self.inputConfig['FIT'][str(p)] = {'NOMINAL':None,'MIN':None,'MAX':None}
                        # self.inputConfig['FIT'][str(p)]['NOMINAL'] = fitResults[param]
                        # self.inputConfig['FIT'][str(p)]['MAX'] = fitResults[param]+sigma*fitResults[param+'err']
                        # self.inputConfig['FIT'][str(p)]['MIN'] = fitResults[param]-sigma*fitResults[param+'err']
                        # self.inputConfig['FIT'][str(p)]['ERROR'] = fitResults[param+'err']

                        paramIndex+=1

                elif 'PFORM' in self.inputConfig['FIT'].keys() or ('XPFORM' in self.inputConfig['FIT'].keys() and 'YPFORM' in self.inputConfig['FIT'].keys()):
                    print 'Resetting fit parameters in input config'
                    self.inputConfig['FIT'] = {'PFORM':funcString}
                    for ix in range(nxparams):
                        for iy in range(nyparams):
                            param = 'X'+str(ix)+'Y'+str(iy)
                            pseudo2D_Rpf.SetParameter(paramIndex,fitResults[param])
                            pseudo2D_Rpf.SetParError(paramIndex,fitResults[param+'err'])
                            if self.fitGuesses != False:
                                self.inputConfig['FIT'][param] = {'NOMINAL':None,'MIN':None,'MAX':None}
                                self.inputConfig['FIT'][param]['NOMINAL'] = fitResults[param]
                                self.inputConfig['FIT'][param]['MAX'] = fitResults[param]+sigma*fitResults[param+'err']
                                self.inputConfig['FIT'][param]['MIN'] = fitResults[param]-sigma*fitResults[param+'err']
                                self.inputConfig['FIT'][param]['ERROR'] = fitResults[param+'err']

                            paramIndex+=1               
                elif 'XFORM' in self.inputConfig['FIT'].keys() and 'YFORM' in self.inputConfig['FIT'].keys():
                    print 'Resetting fit parameters in input config'
                    self.inputConfig['FIT'] = {'FORM':funcString}
                    for ix in range(nxparams):
                        for iy in range(nyparams):
                            param = 'X'+str(ix)+'Y'+str(iy)
                            pseudo2D_Rpf.SetParameter(paramIndex,fitResults[param])
                            pseudo2D_Rpf.SetParError(paramIndex,fitResults[param+'err'])
                            if self.fitGuesses != False:
                                self.inputConfig['FIT'][param] = {'NOMINAL':None,'MIN':None,'MAX':None}
                                self.inputConfig['FIT'][param]['NOMINAL'] = fitResults[param]
                                self.inputConfig['FIT'][param]['MAX'] = fitResults[param]+sigma*fitResults[param+'err']
                                self.inputConfig['FIT'][param]['MIN'] = fitResults[param]-sigma*fitResults[param+'err']
                                self.inputConfig['FIT'][param]['ERROR'] = fitResults[param+'err']

                            paramIndex+=1

                else:
                    print 'Output fit type not supported for fitGuesses. Quittting...'
                    quit()

            # Finally draw the surface
            pseudo2D_results.Close()
            header.makeCan('pseudo2D_Rpf',self.projPath,[pseudo2D_Rpf],xtitle=self.xVarName,ytitle=self.yVarName)

    def _inputOrganizer(self):
        #################################################################################
        # First we need to get histograms from files and store them in a new dictionary #
        #################################################################################
        dict_hists = {}

        # Stores [process,cat] pairs of regions with integral of zero so we can tell the card this
        self.integralZero = []

        # Grab all process names and loop through
        processes = [process for process in self.inputConfig['PROCESS'].keys() if process != "HELP"]
        for process in processes:
            this_process_dict = self.inputConfig['PROCESS'][process]
            
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

            dict_hists[process]['file'] = file_nominal
            dict_hists[process]['pass']['nominal'] = hist_pass
            dict_hists[process]['fail']['nominal'] = hist_fail


            # If there are systematics
            if len(this_process_dict['SYSTEMATICS']) == 0:
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


                    # Only care about syst (right now) if it's a shape (CODE == 2 or 3)
                    if this_syst_dict['CODE'] == 2:   # same file as norm, different hist names
                        dict_hists[process]['pass'][syst+'Up']   = file_nominal.Get(this_syst_dict['HISTPASS_UP'])
                        dict_hists[process]['pass'][syst+'Down'] = file_nominal.Get(this_syst_dict['HISTPASS_DOWN'])
                        dict_hists[process]['fail'][syst+'Up']   = file_nominal.Get(this_syst_dict['HISTFAIL_UP'])
                        dict_hists[process]['fail'][syst+'Down'] = file_nominal.Get(this_syst_dict['HISTFAIL_DOWN'])

                    if this_syst_dict['CODE'] == 3:   # different file as norm and different files for each process if specified, same hist name if not specified in inputConfig
                        # User will most likely have different file for each process but maybe not so check
                        if 'FILE_UP' in this_syst_dict:
                            file_up = TFile.Open(this_syst_dict['FILE_UP'])
                        # DOCUMENT
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
                            dict_hists[process]['pass'][syst+'Up'] = file_up.Get(this_syst_dict['HISTPASS_UP'])            # try to grab hist name from SYSTEMATIC dictionary
                        elif 'HISTPASS' in this_syst_dict:
                            dict_hists[process]['pass'][syst+'Up'] = file_up.Get(this_syst_dict['HISTPASS'])               # else use the same one as nominal distribution
                        elif 'HISTPASS_UP_*' in this_syst_dict:
                            dict_hists[process]['pass'][syst+'Up'] = file_up.Get(this_syst_dict['HISTPASS_UP_*'].replace('*',process))
                        else: 
                            dict_hists[process]['pass'][syst+'Up'] = file_up.Get(this_syst_dict['HISTPASS_UP_'+process])   # or use process specific name

                        if 'HISTPASS_DOWN' in this_syst_dict:
                            dict_hists[process]['pass'][syst+'Down'] = file_down.Get(this_syst_dict['HISTPASS_DOWN'])
                        elif 'HISTPASS' in this_syst_dict:
                            dict_hists[process]['pass'][syst+'Down'] = file_down.Get(this_syst_dict['HISTPASS'])
                        elif 'HISTPASS_DOWN_*' in this_syst_dict:
                            dict_hists[process]['pass'][syst+'Down'] = file_up.Get(this_syst_dict['HISTPASS_DOWN_*'].replace('*',process))
                        else:
                            dict_hists[process]['pass'][syst+'Down'] = file_down.Get(this_syst_dict['HISTPASS_DOWN_' + process])

                        if 'HISTFAIL_UP' in this_syst_dict:
                            dict_hists[process]['fail'][syst+'Up'] = file_up.Get(this_syst_dict['HISTFAIL_UP'])
                        elif 'HISTFAIL' in this_syst_dict:
                            dict_hists[process]['fail'][syst+'Up'] = file_up.Get(this_syst_dict['HISTFAIL'])
                        elif 'HISTFAIL_UP_*' in this_syst_dict:
                            dict_hists[process]['fail'][syst+'Up'] = file_up.Get(this_syst_dict['HISTFAIL_UP_*'].replace('*',process))    
                        else:
                            dict_hists[process]['fail'][syst+'Up'] = file_up.Get(this_syst_dict['HISTFAIL_UP_' + process])

                        if 'HISTFAIL_DOWN' in this_syst_dict:
                            dict_hists[process]['fail'][syst+'Down'] = file_down.Get(this_syst_dict['HISTFAIL_DOWN'])
                        elif 'HISTFAIL' in this_syst_dict:
                            dict_hists[process]['fail'][syst+'Down'] = file_down.Get(this_syst_dict['HISTFAIL'])
                        elif 'HISTFAIL_DOWN_*' in this_syst_dict:
                            dict_hists[process]['fail'][syst+'Down'] = file_up.Get(this_syst_dict['HISTFAIL_DOWN_*'].replace('*',process))
                        else:
                            dict_hists[process]['fail'][syst+'Down'] = file_down.Get(this_syst_dict['HISTFAIL_DOWN_' + process])

        

        #####################################################################
        # With dictionary made, we can split around the signal region and   #
        # start renaming to match the format required by Combine. The       #
        # dictionary key names are conveniently named so we can do this     #
        # with minimumal pain.                                              #
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
            self.organizedDict[process] = {'pass_LOW':{}, 'pass_SIG':{}, 'pass_HIGH':{}, 'fail_LOW':{}, 'fail_SIG':{}, 'fail_HIGH':{}}
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
                    for c in ['LOW','SIG','HIGH']: 
                        # Get new names
                        histname = process + '_' + cat+'_'+c+'_'+self.name
                        if dist != 'nominal':                           # if not nominal dist
                            histname = histname + '_' + dist
                        print 'Making '+histname
                        finalhist = header.copyHistWithNewXbins(temp_hist,self.newXbins[c],histname)#,self.oldYwidth)

                        # Test if histogram is non-zero
                        if finalhist.Integral() <= 0:
                            print 'WARNING: '+process+', '+cat+'_'+c+', '+dist+' has zero or negative events - ' + str(finalhist.Integral())
                            self.integralZero.append([process,cat+'_'+c])
                            # If it is, zero the bins except one to avoid Integral=0 errors in combine
                            for b in range(1,finalhist.GetNbinsX()*finalhist.GetNbinsY()+1):
                                finalhist.SetBinContent(b,1e-10)

                        finalhist.Write()
                        self.organizedDict[process][cat+'_'+c][dist] = finalhist.GetName()#header.copyHistWithNewXbins(temp_hist,self.newXbins[c],histname)


    def _buildFitWorkspace(self):
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
            for cat in ['pass','fail']:
                for c in ['LOW','SIG','HIGH']:
                    for dist in self.organizedDict[process][cat+'_'+c].keys():
                        # For each category
                        Roo_dict[process][cat+'_'+c][dist] = {}
                        var_list = var_lists[c]
                        Roo_dict[process][cat+'_'+c][dist] = {}
                        print 'Making RDH '+self.organizedDict[process][cat+'_'+c][dist]
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
        # Need to build for each category
        for c in ['LOW','SIG','HIGH']:
            bin_list_fail = RooArgList()
            bin_list_pass = RooArgList()

            TH2_data_fail = self.orgFile.Get(self.organizedDict['data_obs']['fail_'+c]['nominal'])
            TH2_data_pass = self.orgFile.Get(self.organizedDict['data_obs']['pass_'+c]['nominal'])

            # Get each bin
            for ybin in range(1,len(self.newYbins)):
                for xbin in range(1,len(self.newXbins[c])):
                    this_full_xbin = self._getFullXbin(xbin,c)
                    # Now that we're in a specific bin, we need to process it

                    # First check if we have an empty pass bin
                    this_pass_bin_zero = True if TH2_data_pass.GetBinContent(xbin,ybin) <= 0 else False 
                    
                    # Now that we know we aren't in the blinded region, make a name for the bin RRV
                    fail_bin_name = 'Fail_'+c+'_bin_'+str(xbin)+'-'+str(ybin)+'_'+self.name
                    pass_bin_name = 'Pass_'+c+'_bin_'+str(xbin)+'-'+str(ybin)+'_'+self.name

                    # Initialize contents by subtracting minor non-QCD components (code 2)
                    bin_content    = TH2_data_fail.GetBinContent(xbin,ybin)
                    bin_range_up   = bin_content + TH2_data_fail.GetBinErrorUp(xbin,ybin)*5
                    bin_range_down = bin_content - TH2_data_fail.GetBinErrorLow(xbin,ybin)*5
                    bin_err_up     = bin_content + TH2_data_fail.GetBinErrorUp(xbin,ybin)
                    bin_err_down   = bin_content - TH2_data_fail.GetBinErrorLow(xbin,ybin)

                    # Now subtract away the parts that we don't want
                    for process in self.organizedDict.keys():
                        this_TH2 = self.orgFile.Get(self.organizedDict[process]['fail_'+c]['nominal'])

                        # Check the code and change bin content and errors accordingly
                        if self.inputConfig['PROCESS'][process]['CODE'] == 0: # signal
                            continue
                        elif self.inputConfig['PROCESS'][process]['CODE'] == 1: # data
                            continue
                        elif self.inputConfig['PROCESS'][process]['CODE'] == 2: # MC
                            bin_content    = bin_content     - this_TH2.GetBinContent(xbin,ybin)
                            bin_range_up   = bin_range_up    - 0.5*this_TH2.GetBinContent(xbin,ybin)
                            bin_range_down = bin_range_down  - 1.5*this_TH2.GetBinContent(xbin,ybin)
                            bin_err_up     = bin_err_up      - this_TH2.GetBinContent(xbin,ybin) + this_TH2.GetBinErrorUp(xbin,ybin)             # Just propagate errors normally
                            bin_err_down   = bin_err_down    - this_TH2.GetBinContent(xbin,ybin) - this_TH2.GetBinErrorLow(xbin,ybin)

                    # If bin content is <= 0, treat this bin as a RooConstVar at value close to 0
                    if (bin_content <= 0) or (this_pass_bin_zero == True):
                        binRRV = RooConstVar(fail_bin_name, fail_bin_name, max(1e-9,bin_content))
                        bin_list_fail.add(binRRV)
                        self.allVars.append(binRRV)

                        # Then get bin center and assign it to a RooConstVar
                        x_center = TH2_data_fail.GetXaxis().GetBinCenter(xbin)
                        y_center = TH2_data_fail.GetYaxis().GetBinCenter(ybin)

                        # Remap to [0,1]
                        x_center_mapped = (x_center - self.newXbins['LOW'][0])/(self.newXbins['HIGH'][-1] - self.newXbins['LOW'][0])
                        y_center_mapped = (y_center - self.newYbins[0])/(self.newYbins[-1] - self.newYbins[0])

                        x_const = RooConstVar("ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,x_center_mapped)
                        y_const = RooConstVar("ConstVar_y_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,y_center_mapped)

                        self.allVars.append(x_const)
                        self.allVars.append(y_const)

                        # # And now get the Rpf function value for this bin 
                        if self.rpf.fitType == 'cheb':
                            this_rpf = self.rpf.evalRpf(x_center_mapped, y_center_mapped,this_full_xbin,ybin) # chebyshev class takes different input
                        else:
                            this_rpf = self.rpf.evalRpf(x_const, y_const,this_full_xbin,ybin)

                        this_bin_pass = RooConstVar(pass_bin_name, pass_bin_name, 1e-9)
                        bin_list_pass.add(this_bin_pass)
                        self.allVars.append(this_bin_pass)

                    else:
                        # Create the fail bin
                        if self.freezeFail:
                            binRRV = RooConstVar(fail_bin_name, fail_bin_name, bin_content)

                        elif bin_content < 1: # Give larger floating to range to bins with fewer events
                            binRRV = RooRealVar(fail_bin_name, fail_bin_name, max(bin_content,0.1), 1e-9, 10)
                            print fail_bin_name + ' < 1'

                        elif bin_content < 10: # Give larger floating to range to bins with fewer events
                            binRRV = RooRealVar(fail_bin_name, fail_bin_name, max(bin_content,1), 1e-9, 50)
                            print fail_bin_name + ' < 10'
                        else:
                            binRRV = RooRealVar(fail_bin_name, fail_bin_name, bin_content, max(1,bin_range_down), bin_range_up)


                        if not self.freezeFail and bin_content >= 10:
                            if bin_content - bin_err_down < 0.001:
                                bin_err_down = bin_content - 0.001#binRRV.getMin()     # For the rare case when bin error is larger than the content
                            
                            binRRV.setAsymError(bin_err_down,bin_err_up)
                            self.floatingBins.append(fail_bin_name)


                        # Store the bin
                        bin_list_fail.add(binRRV)
                        self.allVars.append(binRRV)

                        # Then get bin center and assign it to a RooConstVar
                        x_center = TH2_data_fail.GetXaxis().GetBinCenter(xbin)
                        y_center = TH2_data_fail.GetYaxis().GetBinCenter(ybin)

                        # Remap to [0,1]
                        x_center_mapped = (x_center - self.newXbins['LOW'][0])/(self.newXbins['HIGH'][-1] - self.newXbins['LOW'][0])
                        y_center_mapped = (y_center - self.newYbins[0])/(self.newYbins[-1] - self.newYbins[0])

                        x_const = RooConstVar("ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,x_center_mapped)
                        y_const = RooConstVar("ConstVar_y_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,"ConstVar_x_"+c+'_'+str(xbin)+'-'+str(ybin)+'_'+self.name,y_center_mapped)

                        self.allVars.append(x_const)
                        self.allVars.append(y_const)

                        # And now get the Rpf function value for this bin 
                        if self.rpf.fitType == 'cheb':
                            this_rpf = self.rpf.evalRpf(x_center_mapped, y_center_mapped,this_full_xbin,ybin) # chebyshev class takes different input
                        else:
                            this_rpf = self.rpf.evalRpf(x_const, y_const,this_full_xbin,ybin)

                        formula_arg_list = RooArgList(binRRV,this_rpf)
                        this_bin_pass = RooFormulaVar(pass_bin_name, pass_bin_name, "@0*@1",formula_arg_list)
                        bin_list_pass.add(this_bin_pass)
                        self.allVars.append(this_bin_pass)


            print "Making RPH2Ds"
            Roo_dict['qcd']['fail_'+c] = {}
            Roo_dict['qcd']['pass_'+c] = {}

            Roo_dict['qcd']['fail_'+c]['RPH2D'] = RooParametricHist2D('qcd_fail_'+c+'_'+self.name,'qcd_fail_'+c+'_'+self.name,x_vars[c], y_var, bin_list_fail, TH2_data_fail)
            Roo_dict['qcd']['fail_'+c]['norm']  = RooAddition('qcd_fail_'+c+'_'+self.name+'_norm','qcd_fail_'+c+'_'+self.name+'_norm',bin_list_fail)
            Roo_dict['qcd']['pass_'+c]['RPH2D'] = RooParametricHist2D('qcd_pass_'+c+'_'+self.name,'qcd_pass_'+c+'_'+self.name,x_vars[c], y_var, bin_list_pass, TH2_data_fail)
            Roo_dict['qcd']['pass_'+c]['norm']  = RooAddition('qcd_pass_'+c+'_'+self.name+'_norm','qcd_pass_'+c+'_'+self.name+'_norm',bin_list_pass)

        print "Making workspace..."
        # Make workspace to save in
        self.workspace = RooWorkspace("w_"+self.name)
        for process in Roo_dict.keys():
            for cat in [k for k in Roo_dict[process].keys() if 'file' not in k]:
                if process == 'qcd':
                    rooObj = Roo_dict[process][cat]
                    # if type(rooObj) != dict:
                    #     print "Importing " + rooObj.GetName() + ' from ' + process + ', ' + cat + ', ' +c
                    #     getattr(self.workspace,'import')(rooObj,RooFit.RecycleConflictNodes(),RooFit.Silence())
                    # else:
                    for itemkey in rooObj.keys():
                        print "Importing " + rooObj[itemkey].GetName() + ' from ' + process + ', ' + cat + ', ' + itemkey
                        getattr(self.workspace,'import')(rooObj[itemkey],RooFit.RecycleConflictNodes(),RooFit.Silence())
                
                else:
                    for dist in Roo_dict[process][cat].keys():
                        rooObj = Roo_dict[process][cat][dist]
                        # if type(rooObj) != dict: 
                        #     print "Importing " + rooObj.GetName() + ' from ' + process + ', ' + cat + ', ' +dist+', '+ c
                        #     getattr(self.workspace,'import')(rooObj,RooFit.RecycleConflictNodes(),RooFit.Silence())
                        # else:
                        for itemkey in rooObj.keys():
                            print "Importing " + rooObj[itemkey].GetName() + ' from ' + process + ', ' + cat  +', ' +dist+ ', ' + itemkey
                            getattr(self.workspace,'import')(rooObj[itemkey],RooFit.RecycleConflictNodes(),RooFit.Silence())

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
        kmax = str(len([syst for syst in self.inputConfig['SYSTEMATIC'].keys() if syst != 'HELP']))

        card_new.write('imax '+imax+'\n')      
        card_new.write('jmax '+jmax+'\n')
        card_new.write('kmax '+kmax+'\n')
        card_new.write('-'*120+'\n')

        ##########
        # Shapes #
        ##########
        procs_with_systs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and len(self.inputConfig['PROCESS'][proc]['SYSTEMATICS']) != 0]
        procs_without_systs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and len(self.inputConfig['PROCESS'][proc]['SYSTEMATICS']) == 0]
        procs_without_systs.append('qcd')   # Again, qcd not in the input JSON but needs to be in the combine card!

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
            if self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 0 or self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 1:        # lnN
                syst_type = 'lnN'
            elif self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 2 or self.inputConfig['SYSTEMATIC'][syst]['CODE'] == 3:      # shape
                syst_type = 'shape'
            else:
                print 'Systematic ' + syst + ' does not have one of the four allowed codes (0,1,2,3). Quitting.'
                quit()
            
            syst_lines[syst] = syst + ' ' + syst_type + ' '

        signal_procs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs' and self.inputConfig['PROCESS'][proc]['CODE'] == 0]
        MC_bkg_procs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs' and (self.inputConfig['PROCESS'][proc]['CODE'] == 2 or self.inputConfig['PROCESS'][proc]['CODE'] == 3)]

        all_procs = [proc for proc in self.inputConfig['PROCESS'].keys() if proc != 'HELP' and proc != 'data_obs']
        all_procs.append('qcd')

        for chan in channels:
            for proc in all_procs:
                # Start lines
                bin_line += (chan+' ')
                processName_line += (proc+' ')

                # If we have processes with integral of 0 in certain regions, we need to tell combine this here
                # if proc in [p[0] for p in self.integralZero]:


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
                    # If the systematic is applicable to the process
                    if proc != 'qcd':
                        if syst_line_key in self.inputConfig['PROCESS'][proc]['SYSTEMATICS']:
                            # If symmetric lnN...
                            if self.inputConfig['SYSTEMATIC'][syst_line_key]['CODE'] == 0:
                                thisVal = str(self.inputConfig['SYSTEMATIC'][syst_line_key]['VAL'])
                            # If asymmetric lnN...
                            elif self.inputConfig['SYSTEMATIC'][syst_line_key]['CODE'] == 1:
                                thisVal = str(self.inputConfig['SYSTEMATIC'][syst_line_key]['VALDOWN']) + '/' + str(self.inputConfig['SYSTEMATIC'][syst_line_key]['VALUP'])
                            # If shape...
                            else:
                                thisVal = str(self.inputConfig['SYSTEMATIC'][syst_line_key]['SCALE'])
                        # Otherwise place a '-'
                        else:
                            thisVal = '-'  
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
        for p in self.rpf.rpfVars.keys():
            card_new.write(header.colliMate(self.rpf.rpfVars[p].GetName()+' flatParam\n',22))

        for b in self.floatingBins:
            card_new.write(header.colliMate(b+' flatParam\n',22))
           
        card_new.close() # need to add masking if blinded

    def plotFitResults(self,fittag,simfit=False): # fittag means 'b' or 's'
        allVars = []

        #####################
        #   Get everything  #
        #####################

        # File with histograms and RooFitResult parameters
        if simfit == False:
            post_file = TFile.Open(self.projPath+'/postfitshapes_'+fittag+'.root')
            fd_file = TFile.Open(self.projPath+'/fitDiagnostics.root')
        else:
            post_file = TFile.Open(self.tag+'/postfitshapes_'+fittag+'.root')
            fd_file = TFile.Open(self.tag+'/fitDiagnostics.root')

        fit_result = fd_file.Get('fit_'+fittag)

        x_low = self.newXbins['LOW'][0]
        x_high = self.newXbins['HIGH'][-1]
        y_low = self.newYbins[0]
        y_high = self.newYbins[-1]
        y_nbins = len(self.newYbins)-1

        # Define low, middle, high projection regions for y (x regions defined already via signal region bounds)
        y_turnon_endBin = post_file.Get('pass_LOW_'+self.name+'_prefit/data_obs').ProjectionY().GetMaximumBin()
        y_tail_beginningBin = int((y_nbins - y_turnon_endBin)/2.0 + y_turnon_endBin)
        print 'Finding start and end bin indexes of signal range. Looking for '+str(self.sigStart)+', '+str(self.sigEnd)
        for ix,xwall in enumerate(self.fullXbins):
            if xwall == self.sigStart:
                print 'Assigning start bin as '+str(ix+1)
                x_sigstart_bin = ix+1
            if xwall == self.sigEnd:
                print 'Assigning end bin as '+str(ix)
                x_sigend_bin = ix

        if y_turnon_endBin > y_nbins/2.0:  # in case this isn't a distribution with a turn-on
            y_turnon_endBin = int(round(y_nbins/3.0))
            y_tail_beginningBin = 2*y_turnon_endBin
        y_turnon_endVal = str(int(post_file.Get('pass_LOW_'+self.name+'_prefit/data_obs').GetYaxis().GetBinUpEdge(y_turnon_endBin)))
        y_tail_beginningVal = str(int(post_file.Get('pass_LOW_'+self.name+'_prefit/data_obs').GetYaxis().GetBinLowEdge(y_tail_beginningBin)))
     

        # Final fit signal strength
        if fittag == 's':
            tree_fit_sb = fd_file.Get('tree_fit_sb')
            tree_fit_sb.GetEntry(0)
            signal_strength = tree_fit_sb.r
        else:
            signal_strength = 0

        #####################
        #    Data vs Bkg    #
        #####################

        hist_dict = {}

        # Organize and make any projections or 2D distributions
        for process in [process for process in self.inputConfig['PROCESS'] if process != 'HELP']+['qcd']:
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
                hist_dict[process][cat]['prefit_projx2'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+y_turnon_endVal+'-'+y_tail_beginningVal,     y_turnon_endBin+1,      y_tail_beginningBin,'e')
                hist_dict[process][cat]['prefit_projx3'] = hist_dict[process][cat]['prefit_2D'].ProjectionX(base_proj_name_pre+'projx_'+y_tail_beginningVal+'-'+str(y_high),         y_tail_beginningBin+1,  y_nbins,'e')

                hist_dict[process][cat]['prefit_projy1'] = hist_dict[process][cat]['prefit_LOW'].ProjectionY(base_proj_name_pre+'projy_'+str(x_low)+'-'+str(self.sigStart),          1,                      hist_dict[process][cat]['prefit_LOW'].GetNbinsX(),'e')
                if self.blindedPlots:
                    hist_dict[process][cat]['prefit_projy2'] = hist_dict[process][cat]['prefit_2D'].ProjectionY(base_proj_name_pre+'projy_'+str(self.sigStart)+'-'+str(self.sigEnd), x_sigstart_bin,           x_sigend_bin,'e')
                else:
                    hist_dict[process][cat]['prefit_projy2'] = hist_dict[process][cat]['prefit_SIG'].ProjectionY(base_proj_name_pre+'projy_'+str(self.sigStart)+'-'+str(self.sigEnd),    1,                      hist_dict[process][cat]['prefit_SIG'].GetNbinsX(),'e')
                hist_dict[process][cat]['prefit_projy3'] = hist_dict[process][cat]['prefit_HIGH'].ProjectionY(base_proj_name_pre+'projy_'+str(self.sigEnd)+'-'+str(x_high),          1,                      hist_dict[process][cat]['prefit_HIGH'].GetNbinsX(),'e')

                hist_dict[process][cat]['postfit_projx1'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+str(y_low)+'-'+y_turnon_endVal,           1,                      y_turnon_endBin, 'e')
                hist_dict[process][cat]['postfit_projx2'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+y_turnon_endVal+'-'+y_tail_beginningVal,  y_turnon_endBin+1,      y_tail_beginningBin,'e')
                hist_dict[process][cat]['postfit_projx3'] = hist_dict[process][cat]['postfit_2D'].ProjectionX(base_proj_name_post+'projx_'+y_tail_beginningVal+'-'+str(y_high),      y_tail_beginningBin+1,  y_nbins,'e')

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

        # Add together processes that we want to see as one
        if self.plotTogether != False:
            hist_dict = self.plotProcessesTogether(hist_dict)
            
        process_list = hist_dict.keys()

        # Create lists for the 2D projections (ones you want to see together)
        for process in hist_dict.keys():    # Canvas
            isSignal = (process != 'qcd' and self.inputConfig['PROCESS'][process]['CODE'] == 0)
            twoDList = []         
            for cat in ['fail','pass']:
                for fit in ['prefit', 'postfit']:
                    if isSignal and fittag == 's' and fit == 'postfit':
                        hist_dict[process][cat][fit+'_2D'].Scale(signal_strength)
                        twoDList.append(hist_dict[process][cat][fit+'_2D'])
                    else:
                        twoDList.append(hist_dict[process][cat][fit+'_2D'])

            if isSignal and fittag != 's':
                continue
            else:
                header.makeCan('fit_'+fittag+'/'+process+'_fit'+fittag+'_2D',self.projPath,twoDList,xtitle=self.xVarTitle,ytitle=self.yVarTitle)

        # Invert the last two items (unique to b*) - customize as needed
        process_list[-1],process_list[-2] = process_list[-2],process_list[-1]

        # Get the colors
        colors = []
        for process in process_list:
            if process != 'data_obs':
                if process not in self.inputConfig['PROCESS'].keys():
                    continue
                if (process != 'qcd' and self.inputConfig['PROCESS'][process]['CODE'] != 0):
                    if (process in self.inputConfig['PROCESS'].keys()) and ('COLOR' in self.inputConfig['PROCESS'][process].keys()):
                        colors.append(self.inputConfig['PROCESS'][process]["COLOR"])
                    else:
                        colors.append(None)

        # Put QCD on bottom of stack since it's smooth
        colors = [kYellow]+colors

        # Create lists for makeCan of the projections
        for plotType in ['postfit_projx','postfit_projy']:   # Canvases
            bkgList = []
            dataList = []
            signal_list = []
            
            for cat in ['fail','pass']: # Row 
                for regionNum in range(1,4):    # Column
                    bkg_process_list = []
                    for process in process_list:
                        if process != 'data_obs':
                            if (process != 'qcd' and self.inputConfig['PROCESS'][process]['CODE'] != 0):
                                bkg_process_list.append(hist_dict[process][cat][plotType+str(regionNum)])
                            elif (process != 'qcd' and self.inputConfig['PROCESS'][process]['CODE'] == 0):
                                # if fittag == 'fit_s':
                                hist_dict[process][cat][plotType+str(regionNum)].Scale(signal_strength)
                                signal_list.append(hist_dict[process][cat][plotType+str(regionNum)])
                                # else:
                                    # signal_list.append(hist_dict[process][cat][plotType+str(regionNum)]) #.replace('post','pre')

                        else:
                            dataList.append(hist_dict[process][cat][plotType+str(regionNum)])

                    
                    ## Need to do this ##
                    # if 'x' in plotType:
                    #     hist_dict['qcd'][cat][plotType+str(regionNum)].GetXaxis().Set(len(self.fullXbins)-1,array.array('d',self.fullXbins))
                    # elif 'y' in plotType:
                    #     hist_dict['qcd'][cat][plotType+str(regionNum)].GetXaxis().Set(len(self.newYbins)-1,array.array('d',self.newYbins))
                    ## otherwise a bunch of dumb warnings are thrown about bin limits ##

                    # Put QCD on bottom of stack since it's smooth
                    bkg_process_list = [hist_dict['qcd'][cat][plotType+str(regionNum)]]+bkg_process_list


                    bkgList.append(bkg_process_list)

            # An attempt to debug the warning about adding histograms with different bin limits - seems there's no reason the warning should be coming up
            # for i,data in enumerate(dataList):
            #     print data.GetName()
            #     for b in bkgList[i]:
            #         for i in range(data.GetXaxis().GetXbins().fN):
            #             if not TMath.AreEqualRel(data.GetXaxis().GetXbins().GetAt(i),b.GetXaxis().GetXbins().GetAt(i),0.00000000001):
            #                 print b.GetName() + ' does not agree. Data: '+str(data.GetXaxis().GetXbins().GetAt(i)) +' vs Bkg: '+str(b.GetXaxis().GetXbins().GetAt(i))


            if 'x' in plotType:
                header.makeCan('fit_'+fittag+'/'+plotType+'_fit'+fittag,self.projPath,dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=self.xVarTitle)
                header.makeCan('fit_'+fittag+'/'+plotType+'_fit'+fittag+'_log',self.projPath,dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=self.xVarTitle,logy=True)
            elif 'y' in plotType:
                header.makeCan('fit_'+fittag+'/'+plotType+'_fit'+fittag,self.projPath,dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=self.yVarTitle)
                header.makeCan('fit_'+fittag+'/'+plotType+'_fit'+fittag+'_log',self.projPath,dataList,bkglist=bkgList,signals=signal_list,colors=colors,xtitle=self.yVarTitle,logy=True)

        # Make comparisons for each background process of pre and post fit projections
        for plotType in ['projx','projy']:
            for process in process_list:
                if process != 'data_obs':
                    pre_list = []
                    post_list = []
                    for cat in ['fail','pass']: # Row 
                        for regionNum in range(1,4):    # Column
                            pre_list.append([hist_dict[process][cat]['prefit_'+plotType+str(regionNum)]])  # in terms of makeCan these are "bkg hists"
                            post_list.append(hist_dict[process][cat]['postfit_'+plotType+str(regionNum)])   # and these are "data hists"
                            if process != 'qcd':
                                if 'COLOR' in self.inputConfig['PROCESS'][process].keys():
                                    prepostcolors = [self.inputConfig['PROCESS'][process]['COLOR']]
                                else:
                                    prepostcolors = [0]
                            else:
                                prepostcolors = [kYellow]

                    if 'x' in plotType: header.makeCan('fit_'+fittag+'/'+process+'_'+plotType+'_fit'+fittag,self.projPath,post_list,bkglist=pre_list,colors=prepostcolors,xtitle=self.xVarTitle,datastyle='histe')
                    if 'y' in plotType: header.makeCan('fit_'+fittag+'/'+process+'_'+plotType+'_fit'+fittag,self.projPath,post_list,bkglist=pre_list,colors=prepostcolors,xtitle=self.yVarTitle,datastyle='histe')


        ##############
        #    Rp/f    #
        ##############
        # Need to sample the space to get the Rp/f with proper errors (1000 samples)
        rpf_xnbins = len(self.fullXbins)-1
        rpf_ynbins = len(self.newYbins)-1
        rpf_zbins = [i/1000. for i in range(0,1001)]
        rpf_samples = TH3F('rpf_samples','rpf_samples',rpf_xnbins, array.array('d',self.fullXbins), rpf_ynbins, array.array('d',self.newYbins), len(rpf_zbins)-1, array.array('d',rpf_zbins))# TH3 to store samples
        sample_size = 500

        # Collect all final parameter values
        param_final = fit_result.floatParsFinal()
        coeffs_final = RooArgSet()
        for v in self.rpf.rpfVars.keys():
            coeffs_final.add(param_final.find(v))

        if self.rpf.fitType != 'cheb':
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
                    if coeff_sample.GetName() in self.rpf.rpfVars.keys():
                        self.rpf.setRpfParam(coeff_sample.GetName(), coeff_sample.getValV())
                    coeff_sample = coeffIter_sample.Next()

                # Loop over bins and fill
                for xbin in range(1,rpf_xnbins+1):
                    for ybin in range(1,rpf_ynbins+1):
                        bin_val = 0

                        thisXCenter = rpf_samples.GetXaxis().GetBinCenter(xbin)
                        thisYCenter = rpf_samples.GetYaxis().GetBinCenter(ybin)

                        thisXMapped = (thisXCenter - self.newXbins['LOW'][0])/(self.newXbins['HIGH'][-1] - self.newXbins['LOW'][0])
                        thisYMapped = (thisYCenter - self.newYbins[0])/(self.newYbins[-1] - self.newYbins[0])

                        # Determine the category
                        if thisXCenter > self.newXbins['LOW'][0] and thisXCenter < self.newXbins['LOW'][-1]: # in the LOW category
                            thisxcat = 'LOW'
                        elif thisXCenter > self.newXbins['SIG'][0] and thisXCenter < self.newXbins['SIG'][-1]: # in the SIG category
                            thisxcat = 'SIG'
                        elif thisXCenter > self.newXbins['HIGH'][0] and thisXCenter < self.newXbins['HIGH'][-1]: # in the HIGH category
                            thisxcat = 'HIGH'

                        bin_val = self.rpf.getRpfBinVal(thisxcat,xbin,ybin)

                        rpf_samples.Fill(thisXCenter,thisYCenter,bin_val)

        elif self.rpf.fitType == 'cheb':
            # Import the basis shapes
            cheb_shapes = TFile.Open(globalDir+'basis_plots/basis_shapes.root')
            first_shape_name = cheb_shapes.GetListOfKeys().First().GetName()
            first_shape = cheb_shapes.Get(first_shape_name) # just used to grab binning and such
            cheb_xnbins = first_shape.GetNbinsX()
            cheb_xmin = first_shape.GetXaxis().GetXmin()
            cheb_xmax = first_shape.GetXaxis().GetXmax()
            cheb_ynbins = first_shape.GetNbinsY()
            cheb_ymin = first_shape.GetYaxis().GetXmin()
            cheb_ymax = first_shape.GetYaxis().GetXmax()

            # Loop over samples
            for i in range(sample_size):
                sys.stdout.write('\rSampling '+str(100*float(i)/float(sample_size)) + '%')
                sys.stdout.flush()

                # Randomize the parameters
                param_sample = fit_result.randomizePars()

                # Make TH2 for this sample
                chebSum = TH2F('chebSum','chebSum',cheb_xnbins,cheb_xmin,cheb_xmax,cheb_ynbins,cheb_ymin,cheb_ymax)

                # Grab relevant coefficients and loop over them to sum over the shapes
                chebCoeffs = param_sample.selectByName('ChebCoeff_*x*y*'+suffix)        # Another trick here - if suffix='', this will grab everything including those
                if suffix == '':                                                        # polyCoeffs with the suffix. Need to remove those by explicitely grabbing them
                    chebCoeffsToRemove = param_sample.selectByName('ChebCoeff_*x*y*_*') # and using .remove(collection)
                    chebCoeffs.remove(chebCoeffsToRemove)

                # Looping...
                chebIter = chebCoeffs.createIterator()
                chebCoeff = chebIter.Next()
                while chebCoeff:
                    chebName = chebCoeff.GetName()
                    xLabel = chebName[len('chebCoeff_'):len('chebCoeff_')+2] 
                    yLabel = chebName[len('chebCoeff_'+xLabel):len('chebCoeff_'+xLabel)+2]

                    # Grab and scale the basis shape
                    tempScaled = cheb_shapes.Get('cheb_Tx'+xLabel[-1]+'_Ty'+yLabel[-1]).Clone()
                    tempScaled.Scale(chebCoeffs.find(chebName).getValV())

                    # Add to the sum
                    chebSum.Add(tempScaled)
                    chebCoeff = chebIter.Next()

                for xbin in range(1,chebSum.GetNbinsX()+1):
                    for ybin in range(1,chebSum.GetNbinsY()+1):
                        thisXCenter = rpf_samples.GetXaxis().GetBinCenter(xbin)
                        thisYCenter = rpf_samples.GetYaxis().GetBinCenter(ybin)
                        rpf_samples.Fill(thisXCenter,thisYCenter,chebSum.GetBinContent(xbin,ybin))

                del chebSum

        print '\n'
        rpf_final = TH2F('rpf_final','rpf_final',rpf_xnbins, array.array('d',self.fullXbins), rpf_ynbins, array.array('d',self.newYbins))
        # Now loop over all x,y bin in rpf_samples, project onto Z axis, 
        # get the mean and RMS and set as the bin content and error in rpf_final
        for xbin in range(1,rpf_final.GetNbinsX()+1):
            for ybin in range(1,rpf_final.GetNbinsY()+1):
                temp_projz = rpf_samples.ProjectionZ('temp_projz',xbin,xbin,ybin,ybin)
                rpf_final.SetBinContent(xbin,ybin,temp_projz.GetMean())
                rpf_final.SetBinError(xbin,ybin,temp_projz.GetRMS())

        rpf_c = TCanvas('rpf_c','Post-fit R_{P/F}',800,700)
        rpf_final.Draw('lego')
        rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_lego.pdf','pdf')
        rpf_final.Draw('pe')
        rpf_c.Print(self.projPath+'plots/fit_'+fittag+'/postfit_rpf_errs.pdf','pdf')

        if os.path.isfile(self.projPath+'plots/rpf_comparisons.root'):

            rpf_file = TFile.Open(globalDir+'/plots/rpf_comparisons.root','UPDATE')

            # Do a ratio and diff with pre-fit
            prefit_rpf = rpf_file.Get('RpfToRemap_unit')

            # Ratio
            rpf_ratio = rpf_final.Clone('rpf_ratio_fit'+fittag)
            rpf_ratio.Divide(prefit_rpf)

            # Difference
            rpf_diff = rpf_final.Clone('rpf_diff_fit'+fittag)
            rpf_diff.Add(prefit_rpf,-1)


            # rpf_ratio_c = TCanvas('rpf_ratio_c','Ratio of post-fit to pre-fit R_{P/F}',800,700)
            # rpf_ratio.Draw('surf')
            # rpf_ratio_c.Print(globalDir+'/plots/rpf_post-to-pre_ratio.pdf','pdf')

            rpf_file.cd()
            rpf_final.Write()
            rpf_ratio.Write()
            rpf_diff.Write()

            rpf_file.Close()

    def plotProcessesTogether(self,histDict):
        process_list = histDict.keys()

        for summation in self.plotTogether.keys():  # For each set we're add together
            process_list.append(summation)   # add the name to a list so we can keep track
            hist_dict[summation] = {}
            first_process = self.plotTogether[summation][0]
            self.inputConfig["PROCESS"][summation] = {"COLOR":inputConfig["PROCESS"][first_process]["COLOR"],
                                                     "CODE":inputConfig["PROCESS"][first_process]["CODE"] }

            for cat in hist_dict[first_process].keys():   # for each pass/fail
                hist_dict[summation][cat] = {}
                for reg in hist_dict[first_process][cat].keys():  # and each region
                    first_hist = hist_dict[first_process][cat][reg].Clone(summation+'_'+cat+'_'+reg)    # Clone the "first" histogram and give it the totalProcName
                    for proc in self.plotTogether[summation]:   # for each process in the list of ones we're adding together
                        if proc != first_process: # not the first one since we've cloned that
                            first_hist.Add(hist_dict[proc][cat][reg])    # add it
                    hist_dict[summation][cat][reg] = first_hist    # Put it in the hist_dict


            for proc in self.plotTogether[summation]:
                process_list.remove(proc)

        return histDict

# WRAPPER FUNCTIONS
def runMLFit(twoDs,rMin,rMax,skipPlots=False,plotOn=''):
    # Set verbosity - chosen from first of configs
    verbose = ''
    if twoDs[0].verbosity != False:
        verbose = ' -v '+twoDs[0].verbosity
    
    # Set systematics
    #syst_option = ' -S 0' # Default 0 systematics
    #for twoD in twoDs:
    #    for proc in twoD.inputConfig['PROCESS'].keys():
    #        if type(twoD.inputConfig['PROCESS'][proc]) == dict:
    #            if twoD.inputConfig['PROCESS'][proc]['CODE'] == 0:
    #                if len(twoD.inputConfig['PROCESS'][proc]['SYSTEMATICS']) != 0: # If at any point there's a process
    #                    syst_option = ''                                           # with non-zero systematics, turn them on

    # Set signal strength range - chosen from first of configs
    sig_option = ' --rMin '+rMin+' --rMax '+rMax
    # if twoDs[0].signalOff:
    #     sig_option = ' --rMin 0 --rMax 0'

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
    if plotOn == '':
        if len(twoDs) > 1:
            card_name = 'card_'+twoDs[0].tag+'.txt'
            projDir = twoDs[0].tag
        else:
            card_name = 'card_'+twoDs[0].name+'.txt'
            projDir = twoDs[0].projPath
    else:
        if len(twoDs) > 1:
            card_name = 'card_'+plotOn.split('/')[1]+'.txt'
            projDir = twoDs[0].tag

        else:
            card_name = 'card_'+plotOn.split('/')[1]+'.txt'
            projDir = twoDs[0].projPath

        if plotOn[-1] != '/': plotOn_depth = '../'*(len(plotOn.count('/'))+1)
        else: plotOn_depth = '../'*(len(plotOn.count('/')))

    # Run Combine
    print 'cd '+projDir
    FitDiagnostics_command = 'combine -M FitDiagnostics -d '+card_name+' '+blind_option+' --saveWorkspace --cminDefaultMinimizerStrategy 0' + sig_option +verbose #+syst_option -> now out of date
    print 'Executing '+FitDiagnostics_command

    with header.cd(projDir):
        command_saveout = open('FitDiagnostics_command.txt','w')
        command_saveout.write(FitDiagnostics_command)
        command_saveout.close()

        if os.path.isfile('fitDiagnostics.root'):
            header.executeCmd('rm fitDiagnostics.root')

        header.executeCmd(FitDiagnostics_command)

        if not os.path.isfile('fitDiagnostics.root'):
            print "Combine failed and never made fitDiagnostics.root. Quitting..."
            for i in twoDs:
                del i
            quit()

        diffnuis_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnostics.root --abs -g nuisance_pulls.root'
        header.executeCmd(diffnuis_cmd)

        systematic_analyzer_cmd = 'python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/systematicsAnalyzer.py '+card_name+' --all -f html > systematics_table.html'
        header.executeCmd(systematic_analyzer_cmd)

    # Save out Rp/f to a text file and make a re-run config
    for twoD in twoDs:
        for fittag in ['b','s']:
            param_out = open(twoD.projPath+'rpf_params_'+twoD.name+'fit'+fittag+'.txt','w')
            rerun_config = header.dictCopy(twoD.inputConfig)

            try:
                coeffs_final = TFile.Open(projDir+'/fitDiagnostics.root').Get('fit_'+fittag).floatParsFinal()
                coeffIter_final = coeffs_final.createIterator()
                coeff_final = coeffIter_final.Next()
                while coeff_final:
                    if coeff_final.GetName() in twoD.rpf.rpfVars.keys():
                        # Text file
                        param_out.write(coeff_final.GetName()+': ' + str(coeff_final.getValV()) + ' +/- ' + str(coeff_final.getError())+'\n')
                        # Re run config
                        for k in rerun_config['FIT'].keys():
                            if k in coeff_final.GetName():
                                rerun_config['FIT'][k]['ERROR'] = coeff_final.getError()
                                rerun_config['FIT'][k]['NOMINAL'] = coeff_final.getValV()
                                # if coeff_final.getValV() >= 0:
                                rerun_config['FIT'][k]['MIN'] = coeff_final.getValV()-3*coeff_final.getError()
                                rerun_config['FIT'][k]['MAX'] = coeff_final.getValV()+3*coeff_final.getError()
                                # else:
                                #     rerun_config['FIT'][k]['MAX'] = coeff_final.getValV()*0.5
                                #     rerun_config['FIT'][k]['MIN'] = coeff_final.getValV()*1.5

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
        if plotOn == '':
            with header.cd(projDir):
                bshapes_cmd = 'PostFit2DShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -o postfitshapes_b.root -f fitDiagnostics.root:fit_b --postfit --sampling --print 2> PostFitShapes2D_stderr_b.txt'
                header.executeCmd(bshapes_cmd)
                sshapes_cmd = 'PostFit2DShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -o postfitshapes_s.root -f fitDiagnostics.root:fit_s --postfit --sampling --print 2> PostFitShapes2D_stderr_s.txt'
                header.executeCmd(sshapes_cmd)

        else:
            with header.cd(plotOn):
                print 'Plotting fit result from '+projDir+' onto workspace from '+card_name
                bshapes_cmd = 'PostFit2DShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -o postfitshapes_b.root -f '+plotOn_depth+'fitDiagnostics.root:fit_b --postfit --sampling --print 2> PostFitShapes2D_stderr_b.txt'
                header.executeCmd(bshapes_cmd)
                sshapes_cmd = 'PostFit2DShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -o postfitshapes_s.root -f '+plotOn_depth+'fitDiagnostics.root:fit_s --postfit --sampling --print 2> PostFitShapes2D_stderr_s.txt'
                header.executeCmd(sshapes_cmd)

def runLimit(twoDs,postfitWorkspaceDir,blindData=True,location=''):
    # Set verbosity - chosen from first of configs
    verbose = ''
    if twoDs[0].verbosity != False:
        verbose = ' -v '+twoDs[0].verbosity

    # Set systematics
    syst_option = ''
    for twoD in twoDs:
        for proc in twoD.inputConfig['PROCESS'].keys():
            if type(twoD.inputConfig['PROCESS'][proc]) == dict:
                if twoD.inputConfig['PROCESS'][proc]['CODE'] == 0:
                    if len(twoD.inputConfig['PROCESS'][proc]['SYSTEMATICS']) != 0: # If at any point there's a process
                        syst_option = ''    

    # Set signal strength range
    sig_option = ' --rMin 0 --rMax 5'

    # Run blind (turns off data everywhere) but don't mask (examines signal region)
    if blindData:
        blind_option = ' --run blind'

    # Set the project directory
    if len(twoDs) > 1:
        identifier = twoDs[0].tag
        projDir = twoDs[0].tag
    else:
        identifier = twoDs[0].name
        projDir = twoDs[0].projPath

    card_name = 'card_'+identifier+'.txt'
    # Check if we can import post-fit result made during MLfit step
    if not os.path.isfile(postfitWorkspaceDir+'/fitDiagnostics.root'):
        print 'ERROR: '+postfitWorkspaceDir+'/fitDiagnostics.root does not exist. Please check that run_MLfit.py finished correctly. Quitting...'
        quit()

    # Make a prefit workspace from the data card
    print 'cd '+projDir
    with header.cd(projDir):
        t2w_cmd = 'text2workspace.py -b '+card_name+' -o limitworkspace.root' 
        header.executeCmd(t2w_cmd)

    # Morph workspace according to imported fit result
    prefit_file = TFile(projDir+'/limitworkspace.root','update')
    postfit_w = prefit_file.Get('w')
    fit_result_file = TFile.Open(postfitWorkspaceDir+'/fitDiagnostics.root')
    fit_result = fit_result_file.Get("fit_b")
    postfit_vars = fit_result.floatParsFinal()

    for idx in range(postfit_vars.getSize()):
        par_name = postfit_vars[idx].GetName()
        if postfit_w.var(par_name):
            print 'Setting '+par_name+' to '+str(postfit_vars[idx].getValV())+' +/- '+str(postfit_vars[idx].getError())
            var = postfit_w.var(par_name)
            var.setVal(postfit_vars[idx].getValV())
            var.setError(postfit_vars[idx].getError())

    prefit_file.Close()

    current_dir = os.getcwd()

    aL_cmd = 'combine -M AsymptoticLimits limitworkspace.root -w w --saveWorkspace' +blind_option + syst_option# + sig_option 

    # Run combine if not on condor
    if location == 'local':    
        print 'cd '+projDir
        with header.cd(projDir):
            header.executeCmd(aL_cmd)
    # If on condor, write a script to run (python will then release its memory usage)
    elif location == 'condor':
        # Do all of the project specific shell scripting here
        shell_finisher = open('shell_finisher.sh','w')
        shell_finisher.write('#!/bin/sh\n')
        shell_finisher.write('cd '+projDir+'\n')
        shell_finisher.write(aL_cmd+'\n')
        shell_finisher.write('cd '+current_dir+'\n')
        shell_finisher.write('tar -czvf '+identifier+'.tgz '+projDir+'/\n')
        shell_finisher.write('cp '+identifier+'.tgz ../../../')
        shell_finisher.close()
