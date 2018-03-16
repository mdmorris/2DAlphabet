#####################################################################################################
# input_organizer.py - written by Lucas Corcodilos, 3/7/18                                          #
# --------------------------------------------------------                                          #
# This script uses the input JSON as a map to all of the relevant histograms to the analysis.       #
# It grabs each of the histograms, stores them in nested dictionaries, and then renames them and    #
# writes them out to a single rootfile and sends dict to build_workspace.py.                        #
#####################################################################################################

import ROOT
from ROOT import *
import header
from header import copyHistWithNewXbounds


# inputConfig is a dictionary from the input JSON
def main(inputConfig, blinded):

    # Initialize an output root file
    outfile = TFile('organized_TH2s.root','recreate')

    #################################################################################
    # First we need to get histograms from files and store them in a new dictionary #
    #################################################################################
    dictHists = {}

    # Grab all process names and loop through
    processes = [process for process in inputConfig['PROCESS'].keys() if process != "HELP"]
    for process in processes:
        thisProcessDict = inputConfig['PROCESS'][process]
        
        
        dictHists[process] = {  
            'file': 0,
            'pass': {},
            'fail': {}
        }
       
        if blinded:
            dictHists[process]['passLow'] = {}
            dictHists[process]['passHigh'] = {}
            dictHists[process]['failLow'] = {}
            dictHists[process]['failHigh'] = {}
            

        # Grab nominal pass and fail distributions
        file_nominal = TFile.Open(thisProcessDict['FILE'])
        hist_pass = file_nominal.Get(thisProcessDict['HISTPASS'])
        hist_fail = file_nominal.Get(thisProcessDict['HISTFAIL'])

        dictHists[process]['file'] = file_nominal
        dictHists[process]['pass']['nominal'] = hist_pass
        dictHists[process]['fail']['nominal'] = hist_fail

        if blinded:
            dictHists[process]['passLow']['nominal'] = hist_pass
            dictHists[process]['failLow']['nominal'] = hist_fail
            dictHists[process]['passHigh']['nominal'] = hist_pass
            dictHists[process]['failHigh']['nominal'] = hist_fail

        # If there are systematics
        if len(thisProcessDict['SYSTEMATICS']) == 0:
            print 'No systematics for process ' + process
        else:
            # Loop through them and grab info from inputConfig['SYSTEMATIC']
            for syst in thisProcessDict['SYSTEMATICS']:
                try:
                    thisSystDict = inputConfig['SYSTEMATIC'][syst]

                # Quit if syst does not exist and user does not want to skip
                except:
                    skip = raw_input('No entry named "' + syst + '" exists in the SYSTEMATIC section of the input JSON. Skip it? (y/n)')
                    if skip == 'y' or skip == 'Y':
                        print 'Skipping ' + syst
                    else: 
                        print 'Quiting'
                        quit()


                # Only care about syst (right now) if it's a shape (CODE == 2 or 3)
                if thisSystDict['CODE'] == 2:   # same file as norm, different hist names
                    dictHists[process]['pass'][syst+'Up']   = file_nominal.Get(thisSystDict['HISTPASS_UP'])
                    dictHists[process]['pass'][syst+'Down'] = file_nominal.Get(thisSystDict['HISTPASS_DOWN'])
                    dictHists[process]['fail'][syst+'Up']   = file_nominal.Get(thisSystDict['HISTFAIL_UP'])
                    dictHists[process]['fail'][syst+'Down'] = file_nominal.Get(thisSystDict['HISTFAIL_DOWN'])

                if thisSystDict['CODE'] == 3:   # different file as norm, same hist name if not specified in inputConfig
                    file_up = TFile.Open(thisSystDict['FILE_UP'])
                    file_down = TFile.Open(thisSystDict['FILE_DOWN'])

                    if 'HISTPASS_UP' in thisSystDict:
                        dictHists[process]['pass'][syst+'Up'] = file_up.Get(thisSystDict['HISTPASS_UP'])        # try to grab hist name from SYSTEMATIC dictionary
                    else:
                        dictHists[process]['pass'][syst+'Up'] = file_up.Get(thisProcessDict['HISTPASS'])        # else use the same one as nominal distribution

                    if 'HISTPASS_DOWN' in thisSystDict:
                        dictHists[process]['pass'][syst+'Down'] = file_down.Get(thisSystDict['HISTPASS_DOWN'])
                    else:
                        dictHists[process]['pass'][syst+'Down'] = file_down.Get(thisProcessDict['HISTPASS'])

                    if 'HISTFAIL_UP' in thisSystDict:
                        dictHists[process]['fail'][syst+'Up'] = file_up.Get(thisSystDict['HISTFAIL_UP'])
                    else:
                        dictHists[process]['fail'][syst+'Up'] = file_up.Get(thisProcessDict['HISTFAIL'])

                    if 'HISTFAIL_DOWN' in thisSystDict:
                        dictHists[process]['fail'][syst+'Down'] = file_down.Get(thisSystDict['HISTFAIL_DOWN'])
                    else:
                        dictHists[process]['fail'][syst+'Down'] = file_down.Get(thisProcessDict['HISTFAIL'])


    #####################################################################
    # With dictionary made, we can split around the signal region and   #
    # start renaming to match the format required by Combine. The       #
    # dictionary key names are  conveniently named so we can do this    #
    # with minimumal pain.                                              #
    #####################################################################
    
    # Quickly grab out axis names and binning
    xVarName = inputConfig['BINNING']['X']['NAME']
    yVarName = inputConfig['BINNING']['Y']['NAME']

    newXmin = inputConfig['BINNING']['X']['LOW']
    newXmax = inputConfig['BINNING']['X']['HIGH']
    sigStart = inputConfig['BINNING']['X']['SIGSTART']
    sigEnd = inputConfig['BINNING']['X']['SIGEND']
    newXnbins = inputConfig['BINNING']['X']['NBINS']
    newXwidth = float(newXmax-newXmin)/float(newXnbins)

    newYmin = inputConfig['BINNING']['Y']['LOW']
    newYmax = inputConfig['BINNING']['Y']['HIGH']
    newYnbins = inputConfig['BINNING']['Y']['NBINS']



    # For each process, category, and dist (nominal, systUp, etc)
    for process in dictHists.keys():
        for cat in ['pass','fail']:
            thisProcessCatDict = dictHists[process][cat]
            for dist in thisProcessCatDict.keys():

                # Check if the user has changed the binning in the config file
                oldXmin = thisProcessCatDict[dist].GetXaxis().GetXmin()
                oldXmax = thisProcessCatDict[dist].GetXaxis().GetXmax()
                oldYmin = thisProcessCatDict[dist].GetYaxis().GetXmin()
                oldYmax = thisProcessCatDict[dist].GetYaxis().GetXmax()
                
                # If the edges have changed, print an error and quit
                if (oldXmin != newXmin) or (oldXmax != newXmax) or (oldYmin != newYmin) or (oldYmax != newYmax):
                    print "Error! Bin edges in input histogram " + histname+" are different than those in the input JSON. This is not currently supported. Exiting."
                    quit()

                # If the number of bins have changed, rebin
                if thisProcessCatDict[dist].GetNbinsX() != newXnbins or thisProcessCatDict[dist].GetNbinsY() != newYnbins:
                    print "Applying new bins: ["+str(thisProcessCatDict[dist].GetNbinsX())+","+str(thisProcessCatDict[dist].GetNbinsY())+"] -> ["+str(newXnbins)+","+str(newYnbins)+"]"
                    ratio_X = int(float(thisProcessCatDict[dist].GetNbinsX())/float(newXnbins))
                    ratio_Y = int(float(thisProcessCatDict[dist].GetNbinsY())/float(newYnbins))
                    thisProcessCatDict[dist].Rebin2D(ratio_X,ratio_Y)

                # Get new names
                histname = process + '_' + cat
                if dist != 'nominal':                           # if not nominal dist
                    histname = histname + '_' + dist
                

                # Name the nominal histograms
                thisProcessCatDict[dist].SetName(histname)
                thisProcessCatDict[dist].SetTitle(histname)

                outfile.cd()
                thisProcessCatDict[dist].Write()

                if blinded:
                    low_histname = process + '_' + cat + 'Low'
                    high_histname = process + '_' + cat + 'High'
                    if dist != 'nominal':                           # if not nominal dist
                        low_histname = low_histname + '_' + dist
                        high_histname = high_histname + '_' + dist

                    # Create the split histograms (do the naming for them in this step)
                    hist_to_split = dictHists[process][cat][dist]
                    low_hist = copyHistWithNewXbounds(hist_to_split,low_histname,newXwidth,newXmin,sigStart)
                    high_hist = copyHistWithNewXbounds(hist_to_split,high_histname,newXwidth,sigEnd,newXmax)
                    # Save them
                    dictHists[process][cat+'Low'][dist] = low_hist
                    dictHists[process][cat+'High'][dist] = high_hist

                    outfile.cd()
                    dictHists[process][cat+'Low'][dist].Write()
                    dictHists[process][cat+'High'][dist].Write()

    outfile.Close()
    return dictHists