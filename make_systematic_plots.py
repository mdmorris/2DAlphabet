#################################################################################
# A fairly simple script that just grabs the nominal, up, and down shapes for   #
# each systematic, plots them together and saves out the plots                  #
#################################################################################

import ROOT
from ROOT import *

def main(inputConfig,tag):
    for proc in inputConfig['PROCESS'].keys():
        if proc != 'HELP':
            for syst in inputConfig['PROCESS'][proc]['SYSTEMATICS']:
                # For each systematic in each process

                print proc + '_'+syst

                # If code 2 or 3 (shape based), grab the up and down shapes for passing and project onto the Y axis
                if inputConfig['SYSTEMATIC'][syst]['CODE'] == 2:
                    thisFile = TFile.Open(inputConfig['PROCESS'][proc]['FILE'])
                    passUp2D = thisFile.Get(inputConfig['SYSTEMATIC'][syst]['HISTPASS_UP'])
                    passDown2D = thisFile.Get(inputConfig['SYSTEMATIC'][syst]['HISTPASS_DOWN'])
                    passUp = passUp2D.ProjectionY(proc + '_'+syst,21,28)
                    passDown = passDown2D.ProjectionY(proc + '_'+syst+'_down',21,28)

                elif inputConfig['SYSTEMATIC'][syst]['CODE'] == 3:
                    if 'FILE_UP_'+proc in inputConfig['SYSTEMATIC'][syst]:
                        thisFileUp = TFile.Open(inputConfig['SYSTEMATIC'][syst]['FILE_UP_'+proc])
                        thisFileDown = TFile.Open(inputConfig['SYSTEMATIC'][syst]['FILE_DOWN_'+proc])

                    elif 'FILE_UP_*' in inputConfig['SYSTEMATIC'][syst]:
                        thisFileUp = TFile.Open(inputConfig['SYSTEMATIC'][syst]['FILE_UP_*'].replace('*',proc))
                        thisFileDown = TFile.Open(inputConfig['SYSTEMATIC'][syst]['FILE_DOWN_*'].replace('*',proc))

                    else:
                        print 'Could not identify file for ' + proc +', '+syst

                    passUp = thisFileUp.Get(inputConfig['SYSTEMATIC'][syst]['HISTPASS']).ProjectionY(proc + '_'+syst,21,28)
                    passDown = thisFileDown.Get(inputConfig['SYSTEMATIC'][syst]['HISTPASS']).ProjectionY(proc + '_'+syst+'_down',21,28)

                else:
                    continue

                # Setup the nominal shape (consistent across all of the pltos)
                fileNom = TFile.Open(inputConfig['PROCESS'][proc]['FILE'])
                passNom = fileNom.Get(inputConfig['PROCESS'][proc]['HISTPASS']).ProjectionY(proc + '_'+syst+'_nom',21,28)

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

                passNom.SetXTitle(inputConfig['BINNING']['Y']['TITLE'])
                passUp.SetXTitle(inputConfig['BINNING']['Y']['TITLE'])
                passDown.SetXTitle(inputConfig['BINNING']['Y']['TITLE'])
                passNom.SetTitle(proc + ' - ' + syst + ' uncertainty')


                passNom.Draw('hist')
                passUp.Draw('same hist')
                passDown.Draw('same hist')
                

                thisCan.Print(tag+'/UncertPlots/Uncertainty_'+proc+'_'+syst+'.pdf','pdf')

                


