import ROOT
from ROOT import *

gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)

def main(inputConfig, organizedDict, blinded):
    allVars = []

    #####################
    #   Get everything  #
    #####################

    # Open up our files and workspaces
    new_file = TFile.Open('MaxLikelihoodFitResult.root')
    new_w = new_file.Get('MaxLikelihoodFitResult')

    old_file = TFile.Open('base_QCDMC.root')    # Need to get the data_obs pass and fail
    old_w = old_file.Get('w_2D')                # which is not stored in the output of Combine
                                                # (the sum of the two is stored there)

    # Build another dictionary that can sort and save our pdfs and RDHs
    # based upon info from the config
    new_roo_dict = dictStructureCopy(organizedDict)

    if blinded:
        catagories = ['passLow','passHigh','failLow', 'failHigh']
    else:
        catagories = ['pass','fail']

    # For each process, catagory, and distribution...
    for proc in new_roo_dict.keys():
        for cat in new_roo_dict[proc].keys():
            if cat not in catagories:                   # Remove any keys that aren't part of catagories
                del new_roo_dict[proc][cat]             # so we don't accidentally access them later
            for dist in new_roo_dict[proc][cat]:
                # Grab the correct Combine output

                if inputConfig['PROCESS'][proc]['CODE'] == 0:       # If signal
                    new_roo_dict[proc][cat][dist]['PDF'] = new_w.pdf('shapeSig_'+cat+'_'+proc+'_morph')
                    new_roo_dict[proc][cat][dist]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_'+proc)     # normalization

                elif inputConfig['PROCESS'][proc]['CODE'] == 1:     # If data
                    new_roo_dict[proc][cat][dist]['RDH'] = new_w.data('data_obs')

                elif inputConfig['PROCESS'][proc]['CODE'] == 2:     # If unchanged MC bkg
                    new_roo_dict[proc][cat][dist]['PDF'] = new_w.pdf('shapeBkg_'+cat+'_'+proc+'_morph')
                    new_roo_dict[proc][cat][dist]['NORM'] = new_w.function()

                elif inputConfig['PROCESS'][proc]['CODE'] == 3:     # If renormalized MC bkg
                    new_roo_dict[proc][cat][dist]['PDF'] = new_w.pdf('shapeSignal_'+cat+'_'+proc+'_morph')
                    new_roo_dict[proc][cat][dist]['NORM'] = new_w.function('n_exp_final_bin'+cat+'_proc_'+proc)     # normalization

                else: 
                    print 'Process ' + proc + ' has code ' + str(inputConfig['PROCESS'][proc]['CODE']) + ' in the input configuration which is not valid. Quitting...'
                    quit()

