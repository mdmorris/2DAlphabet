import sys, os, time, decimal, pickle
import header
import random
import ROOT
from ROOT import *

gStyle.SetOptStat(0)
gStyle.SetOptFit(1)
gStyle.SetOptTitle(0)
gROOT.SetBatch()

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--projDir", dest="projDir",
                help="Home of the project - has the cards, fit results, etc")
parser.add_option("-a", "--altDir", dest="altDir",
                help="Home of the alternative model that you'd like to compare against the one in projDir")
parser.add_option("-t", "--toys", dest="toys",default='100',
                help="Number of toys to generate - fed to combine")
parser.add_option('--rMin', metavar='F', type='string', action='store',
                default =   '0',
                dest    =   'rMin',
                help    =   'Minimum bound on r (signal strength)')
parser.add_option('--rMax', metavar='F', type='string', action='store',
                default =   '5',
                dest    =   'rMax',
                help    =   'Minimum bound on r (signal strength)')
# parser.add_option("-b", "--blind", 
#                 action="store_false", dest="blind", default=True,
#                 help="Blinds the signal region in the fit")
parser.add_option("--plotOnly", 
                action="store_true", dest="plotOnly", default=False,
                help="Only plot")
parser.add_option("--dryrun", 
                action="store_true", dest="dryrun", default=False,
                help="Dry run the combine commands to console")
parser.add_option("--gof", 
                action="store_true", dest="gof", default=False,
                help="Perform goodness of fit test")
parser.add_option("--signalInjection", 
                action="store", dest="signalInjection", default='',
                help="Perform signal injection test")
parser.add_option("--biasStudy", 
                action="store_true", dest="biasStudy", default=False,
                help="Perform bias study")
parser.add_option("--ftest", 
                action="store_true", dest="ftest", default=False,
                help="Perform F-test")
parser.add_option("--diagnosticsWithToys", 
                action="store_true", dest="diagnosticsWithToys", default=False,
                help="Perform diagnostics with toys")
parser.add_option("--post", 
                action="store_true", dest="post", default=False,
                help="Run in conjunction with diagnosticsWithToys or signalInjection to create plots")
(options, args) = parser.parse_args()

projDir = options.projDir # home of the workspace - has the cards, fit results, etc
tag = projDir.split('/')[0]
if projDir.split('/')[-1] != '': card_tag = projDir.split('/')[-1]
else: card_tag = projDir.split('/')[-2]

if tag == '':
    print 'ERROR in project directory name (where your workspace and data card lives). Did you accidentally provide a leading slash? (ie /projDir/) Quitting...'
    quit()

if not os.path.isdir(projDir): 
    print projDir +' is not a directory. Quitting...'
    quit()

# Start with tests that only require the projDir and don't do a comparison
with header.cd(projDir):
    #######
    # GOF #
    #######
    if options.gof:
        if not options.plotOnly:
            commands = []

            # Determine if there are channel masks from the fit in the projDir
            masked_regions = []
            f = TFile.Open('higgsCombineTest.FitDiagnostics.mH120.root')
            w = f.Get('w')
            allVars = RooArgList(w.allVars())

            # Loop over all vars in original workspace and add any masked region names to list
            for i in range(allVars.getSize()):
                if 'mask_pass_SIG_' in allVars[i].GetName():
                    if allVars[i].getValV() == 1:
                        masked_regions.append(allVars[i].GetName())

            f.Close()

            # Set the string to specify the blinding
            blind_string = "--setParametersForFit "
            for m in masked_regions:
                blind_string+=m+'=1,'
            blind_string = blind_string[:-1]+' '
            blind_string = blind_string+blind_string.replace('setParametersForFit','setParametersForEval')

            seed = random.randint(100000,999999)        
            commands.append('combine -M GoodnessOfFit card_'+card_tag+'.txt --text2workspace "--channel-masks" -m 120 --algo=saturated '+blind_string)
            commands.append('combine -M GoodnessOfFit card_'+card_tag+'.txt --text2workspace "--channel-masks" -m 120 --algo=saturated '+blind_string+' -t '+options.toys+' -s '+str(seed))
            
            # Run commands
            for c in commands:
                header.executeCmd(c,options.dryrun)

        # Now to analyze the output

        # Get observation
        gofOutput = TFile.Open('higgsCombineTest.GoodnessOfFit.mH120.root')
        gofLimitTree = gofOutput.Get('limit')
        gofLimitTree.GetEntry(0)
        gofLimit = gofLimitTree.limit

        # Get toys
        toyOutput = TFile.Open('higgsCombineTest.GoodnessOfFit.mH120.'+str(seed)+'.root')
        toyLimitTree = toyOutput.Get('limit')
        toyLimitTree.Draw('limit>>hlimit') 
        toyLimits = gDirectory.Get('hlimit')
        time.sleep(1) # if you don't sleep the code moves too fast and won't perform the fit
        toyLimits.Fit("gaus")

        # Fit toys and derive p-value
        gaus = toyLimits.GetFunction("gaus")
        pvalue = 1-(1/gaus.Integral(-float("inf"),float("inf")))*gaus.Integral(-float("inf"),gofLimit)

        # Write out for reference
        out = open('gof_results.txt','w')
        out.write('Test statistic in data = '+str(gofLimit))
        out.write('Mean from toys = '+str(gaus.GetParameter(1)))
        out.write('Width from toys = '+str(gaus.GetParameter(2)))
        out.write('p-value = '+str(pvalue))

        # Extend the axis if needed
        if toyLimits.GetXaxis().GetXmax() < gofLimit:
            print 'Axis limit greater than GOF t value'
            binwidth = toyLimits.GetXaxis().GetBinWidth(1)
            xmin = toyLimits.GetXaxis().GetXmin()
            new_xmax = int(gofLimit*1.1)
            new_nbins = int((new_xmax-xmin)/binwidth)
            toyLimitTree.Draw('limit>>hlimitrebin('+str(new_nbins)+', '+str(xmin)+', '+str(new_xmax)+')') 
            toyLimits = gDirectory.Get('hlimitrebin')
            toyLimits.Fit("gaus")
            gaus = toyLimits.GetFunction("gaus")

        # Arrow for observed
        arrow = TArrow(gofLimit,0.25*toyLimits.GetMaximum(),gofLimit,0)
        arrow.SetLineWidth(2)

        # Legend
        leg = TLegend(0.6,0.6,0.89,0.89)
        leg.SetLineColor(kWhite)
        leg.SetLineWidth(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.AddEntry(toyLimits,"toy data","lep")
        leg.AddEntry(arrow,"observed = %.1f"%gofLimit,"l")
        leg.AddEntry(0,"p-value = %.2E"%decimal.Decimal(pvalue),"")

        # Draw
        cout = TCanvas('cout','cout',800,700)
        toyLimits.Draw('pez')
        arrow.Draw()
        leg.Draw()

        cout.Print('gof_plot.pdf','pdf')

    #########################################
    # Toy diagnositics and signal injection #
    #########################################
    if options.diagnosticsWithToys or options.signalInjection != '':
        if options.diagnosticsWithToys:
            expectSignal = '0'
            run_name = 'diagnosticsWithToys'
        elif options.signalInjection != '':
            expectSignal = options.signalInjection
            run_name = 'signalInjection'+expectSignal
        if not options.plotOnly:
            projtag = projDir.split('/')[0]

            ########################################################################
            # First morph the base workspace to post-fit according to MLfit result #
            ########################################################################
            # Check if we can import post-fit result made during MLfit step
            if not os.path.isfile('fitDiagnostics.root'):
                print 'ERROR: '+projDir+'/fitDiagnostics.root does not exist. Please check that run_MLfit.py finished correctly. Quitting...'
                quit()

            # Make a prefit workspace from the data card
            t2w_cmd = 'text2workspace.py -b card_'+card_tag+'.txt -o '+run_name+'workspace.root' 
            header.executeCmd(t2w_cmd,options.dryrun)

            # Morph workspace according to imported fit result
            prefit_file = TFile(run_name+'workspace.root','update')
            postfit_w = prefit_file.Get('w')
            fit_result_file = TFile.Open('fitDiagnostics.root')
            fit_result = fit_result_file.Get("fit_b")
            postfit_vars = fit_result.floatParsFinal()

            for idx in range(postfit_vars.getSize()):
                par_name = postfit_vars[idx].GetName()
                if postfit_w.var(par_name):
                    
                    var = postfit_w.var(par_name)
                    var.setVal(postfit_vars[idx].getValV())
                    # for pn in ['fullPoly','splitPoly','cheb','generic']:
                    #     if pn in par_name:
                    #         var.setError(abs(min(postfit_vars[idx].getError(),postfit_vars[idx].getValV())))
                    #         print 'Setting '+par_name+' to '+str(postfit_vars[idx].getValV())+' +/- '+str(abs(min(postfit_vars[idx].getError(),postfit_vars[idx].getValV())))
                    #     else:
                    if 'Fail_' in par_name and (postfit_vars[idx].getValV() - postfit_vars[idx].getError()) < 0:
                        var.setError(postfit_vars[idx].getValV()-0.001)
                    else:
                        var.setError(postfit_vars[idx].getError())
                    print 'Setting '+par_name+' to '+str(var.getValV())+' +/- '+str(var.getError())


            prefit_file.Close()

            #########################################
            # Now generate toys from this workspace #
            #########################################
            seed = random.randint(100000,999999)
            # gen_command = 'combine card_toygen.txt -M GenerateOnly -t '+options.toys+' --expectSignal '+expectSignal+'  --saveToys -m 120 -s '+str(seed)+' -n '+run_name
            # header.executeCmd(gen_command,options.dryrun)
            fit_command = 'combine '+run_name+'workspace.root -M FitDiagnostics --cminDefaultMinimizerStrategy 0 --expectSignal '+expectSignal+' -t '+options.toys+' --rMin '+options.rMin+' --rMax '+options.rMax+' --noErrors --minos none --skipBOnlyFit -n '+run_name
            header.executeCmd(fit_command,options.dryrun)

        ################################
        # Plot and save out the result #
        ################################
        result_can = TCanvas('sigpull_can','sigpull_can',800,700)
        fitdaig_out = TFile.Open('fitDiagnostics'+run_name+'.root')
        tree_fit_sb = fitdaig_out.Get('tree_fit_sb')
        tree_fit_sb.Draw("(r-"+expectSignal+")/rErr>>sigpull(20,-5,5)")
        tree_fit_sb.Draw("(r-"+expectSignal+")>>sigstrength(20,-5,5)")
        hsigpull = gDirectory.Get('sigpull')
        hsignstrength = gDirectory.Get('sigstrength')

        hsigpull.Fit("gaus")
        hsigpull.SetTitle(run_name)
        hsigpull.GetXaxis().SetTitle('(r-'+expectSignal+')/rErr')
        result_can.cd()
        hsigpull.Draw('pe')
        result_can.Print(run_name+'_sigpull.pdf','pdf')

        hsignstrength.Fit("gaus")
        hsignstrength.SetTitle(run_name)
        hsignstrength.GetXaxis().SetTitle('r-'+expectSignal)
        result_can.cd()
        hsignstrength.Draw('pe')
        result_can.Print(run_name+'_sigstrength.pdf','pdf')


            # command_to_diagnose = open('FitDiagnostics_command.txt','r').readline()
            # command_to_diagnose += ' --toysFile -n '+run_name+' -t '+options.toys+' --toysFrequentist --expectSignal 0 --noErrors --minos none'
        # Need to write bit to plot
        # else:

 

# Can run against other models or one model against itself to do a b-only vs s+b comparison
if options.biasStudy or options.ftest:
    altDir = options.altDir # home of the alternate workspace - has the cards, fit results, etc
    alttag = altDir.split('/')[0]
    if altDir.split('/')[-1] != '': altcard_tag = altDir.split('/')[-1]
    else: altcard_tag = altDir.split('/')[-2]

    if not os.path.isdir(altDir): 
        print altDir +' is not a directory. Quitting...'
        quit()

    ##############
    # Bias study #
    ##############
    if options.biasStudy:
        ########################################################################
        # First morph the base workspace to post-fit according to MLfit result #
        ########################################################################
        # Check if we can import post-fit result made during MLfit step
        if not os.path.isfile(projDir+'/fitDiagnostics.root'):
            print 'ERROR: '+projDir+'/fitDiagnostics.root does not exist. Please check that run_MLfit.py finished correctly. Quitting...'
            quit()

        # Make a prefit workspace from the data card
        print 'cd '+projDir
        with header.cd(projDir):
            t2w_cmd = 'text2workspace.py -b card_'+card_tag+'.txt -o biasworkspace.root' 
            header.executeCmd(t2w_cmd,options.dryrun)

            # Morph workspace according to imported fit result
            prefit_file = TFile('biasworkspace.root','update')
            postfit_w = prefit_file.Get('w')
            fit_result_file = TFile.Open('fitDiagnostics.root')
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

            #########################################
            # Now generate toys from this workspace #
            #########################################
            seed = random.randint(100000,999999)
            gen_command = 'combine -M biasworkspace.root GenerateOnly -t '+options.toys+' --expectSignal 0 --bypassFrequentistFit --saveToys -m 120 -s '+str(seed)
            header.executeCmd(gen_command)

        ######################################################################
        # Go into the alternative model directory and run the toy fits there #
        ######################################################################
        if altDir[-1] != '/': altDir_depth = '../'*(len(altDir.count('/'))+1)
        else: altDir_depth = '../'*(len(altDir.count('/')))

        print 'cd '+altDir
        with header.cd(altDir):
            fit_command = 'combine card_'+altcard_tag+'.txt -M FitDiagnostics --toysFile '+altDir_depth+projDir+'/higgsCombineTest.GenerateOnly.mH120.'+str(seed)+'.root  -t '+options.toys+' --rMin '+options.rMin+' --rMax '+options.rMax
            header.executeCmd(fit_command,options.dryrun)

            ################################
            # Plot and save out the result #
            ################################
            fitdaig_out = TFile.Open('fitDiagnostics.root')
            tree_fit_sb = fitdaig_out.Get('fitdaig_out')
            tree_fit_sb.Draw("(r-1)/rErr>>sigpull(20,-4,4)")
            hout = gDirectory.Get('sigpull')
            hout.Fit("gaus")

            hout.SaveAs('biasAgainst_'+tag+'.pdf','pdf')


    ##########
    # F test #
    ##########
    elif options.ftest:
        seed = random.randint(100000,999999)  
        if altDir[-1] != '/': altDir_depth = '../'*(len(altDir.count('/'))+1)
        else: altDir_depth = '../'*(len(altDir.count('/')))

        # Do some basic checks before wasting compute time
        base_nrpf_params, base_nbins = header.ftestInfoLookup(header.projInfoLookup(projDir,card_tag))
        alt_nrpf_params, alt_nbins = header.ftestInfoLookup(header.projInfoLookup(altDir,altcard_tag))
        toy_fit_filename = 'higgsCombineTest.GoodnessOfFit.mH120.'+seed+'.root'
        base_fit_filename = 'higgsCombineTest.GoodnessOfFit.mH120.root'

        # If the number of bins in the two models doesn't match, specify which to use or quit
        if base_nbins != alt_nbins:
            error_input = raw_input('ERROR: number of bins in the two models does not match (%i vs %i). Please repeat back which number you would like to calculate for or enter any other string to abort.')
            if int(error_input) == base_nbins:
                ftest_nbins = base_nbins
            elif int(error_input) == alt_nbins:
                ftest_nbins = alt_nbins
            else:
                print 'Quitting...'
                quit()
        else:
            ftest_nbins = base_nbins

        # Steps:
        # 1) Run GOF in base and alt models
        # 2) Run frequentist toy generation (from data) from base model
        # 3) Run GOF for both models over toys
        # 4) Compare results from steps 1 and 3
        # Base
        with header.cd(projDir):
            base_fit_cmd = 'combine card_'+card_tag+'.txt -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 --rMax 10.0 --rMin -10.0'
            header.executeCmd(base_fit_cmd)
        # Alt
        with header.cd(altDir):
            alt_fit_cmd = 'combine card_'+altcard_tag+'.txt -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 --rMax 20.0 --rMin -20.0'
            header.executeCmd(alt_fit_cmd)
        # Generation
        with header.cd(projDir):
            gen_cmd = 'combine card_'+card_tag+'.txt -M GenerateOnly --rMax 10.0 --rMin -10.0 --toysFrequentist -t '+options.toys+' --expectSignal 0 --bypassFrequentistFit --saveToys --seed '+seed
            header.executeCmd(gen_cmd)
        # Base fit to toys
            base_toy_fit_cmd = 'combine card_'+card_tag+'.txt -M GoodnessOfFit --rMax 10.0 --rMin -10.0 --algo saturated -t '+options.toys+' --toysFile higgsCombineTest.GenerateOnly.mH120.'+seed+'.root'
            header.executeCmd(toy_fit_cmd)
        # Alt fit to toys
        with header.cd(altDir):
            alt_toy_fit_cmd = 'combine card_'+altcard_tag+'.txt -M GoodnessOfFit --rMax 20.0 --rMin -20.0 --algo saturated -t '+options.toys+' --toysFile '+altDir_depth+projDir+'/higgsCombineTest.GenerateOnly.mH120.'+seed+'.root'
            header.executeCmd(alt_toy_fit_cmd)


        print 'Analyzing F-test results...'
        base_fstat = header.FStatCalc(projDir+"/"+base_fit_filename, altDir+"/"+base_fit_filename, base_nrpf_params, alt_nrpf_params, ftest_nbins)
        toys_fstat = header.FStatCalc(projDir+"/"+toy_fit_filename, altDir+"/"+toy_fit_filename, base_nrpf_params, alt_nrpf_params, ftest_nbins)

        if len(base_fstat) == 0: base_fstat = [0.0]
        print "base_fstat: ",base_fstat

        toy_pass = 0
        for toyval in toys_fstat:
            print 'toys_fstat vs base_fstat:',toyval,base_fstat[0]
            if base_fstat[0] > toyval:
                toy_pass+=1

        pval = 1
        if len(toys_fstat) > 0:
            pval = float(toy_pass)/float(len(toys_fstat))
            print "F Test p-value",pval

        print "passing toys/number of toys = " + str(float(toy_pass)/float(len(toys_fstat)))

        # Now we plot
        c = TCanvas('c','c',800,600)    
        c.SetLeftMargin(0.12) 
        c.SetBottomMargin(0.12)
        c.SetRightMargin(0.1)
        c.SetTopMargin(0.1)
        ftestHist_nbins = 70.
        ftestHist = TH1F(tag+"Fhist",tag+"Fhist",ftestHist_nbins,0,max(max(toys_fstat),iCentral)+1)
        ftestHist_cut = TH1F(tag+"Fhist_cut",tag+"Fhist cut",ftestHist_nbins,0,max(max(toys_fstat),iCentral)+1)
        ftestHist.GetXaxis().SetTitle("F = #frac{-2log(#lambda_{1}/#lambda_{2})/(p_{2}-p_{1})}{-2log#lambda_{2}/(n-p_{2})}")
        ftestHist.GetXaxis().SetTitleSize(0.025)
        ftestHist.GetXaxis().SetTitleOffset(2)
        ftestHist.GetYaxis().SetTitle("Pseudodatasets")
        ftestHist.GetYaxis().SetTitleOffset(0.85)
        
        for toyval in toys_fstat:
            ftestHist.Fill(toyval)
            if toyval > base_fstat[0]:
                ftestHist_cut.Fill(toyval)

        ftestHist.SetMarkerStyle(20)
        ftestHist.Draw("pez")
        ftestobs  = TArrow(base_fstat[0],0.25*ftestHist.GetMaximum(),base_fstat[0],0)
        ftestobs.SetLineColor(kBlue+1)
        ftestobs.SetLineWidth(2)
        ftestHist_cut.SetLineColor(kViolet-10)
        ftestHist_cut.SetFillColor(kViolet-10)
        ftestHist_cut.Draw("histsame")

        ftest_p1 = min(base_nrpf_params,alt_nrpf_params)
        ftest_p2 = max(base_nrpf_params,alt_nrpf_params)
        ftest_nbins = base_nbins
        fdist = TF1("fDist", "[0]*TMath::FDist(x, [1], [2])", 0,max(max(toys_fstat),base_fstat[0])+1)
        fdist.SetParameter(0,ftestHist.Integral()*((max(max(toys_fstat),base_fstat[0])+1)/ftestHist_nbins))
        fdist.SetParameter(1,ftest_p2-ftest_p1)
        fdist.SetParameter(2,ftest_nbins-ftest_p2)
        fdist.Draw('same')
      
        ftestHist.Draw("pezsame")
        ftestobs.Draw()
        tLeg = TLegend(0.6,0.6,0.89,0.89)
        tLeg.SetLineColor(kWhite)
        tLeg.SetLineWidth(0)
        tLeg.SetFillStyle(0)
        tLeg.SetTextFont(42)
        tLeg.AddEntry(ftestHist,"toy data","lep")
        tLeg.AddEntry(ftestobs,"observed = %.1f"%iCentral,"l")
        tLeg.AddEntry(ftestHist_cut,"p-value = %.2f"%(1-prob),"f")
        tLeg.AddEntry(fdist,"F-dist, ndf = (%.0f, %.0f) "%(fdist.GetParameter(1),fdist.GetParameter(2)),"l")
        tLeg.Draw("same")
        latex = TLatex()
        latex.SetTextAlign(11)
        latex.SetTextSize(0.06)
        latex.SetTextFont(62)
        latex.SetNDC()
        latex.DrawLatex(0.12,0.91,"CMS")
        latex.SetTextSize(0.05)
        latex.SetTextFont(52)
        # if options.isData:
        l.DrawLatex(0.23,0.91,"Preliminary")
        # else:
        #     l.DrawLatex(0.23,0.91,"Simulation")
        latex.SetTextFont(42)
        # latex.DrawLatex(0.76,0.91,"%.1f fb^{-1}"%options.lumi)
        latex.SetTextFont(52)
        latex.SetTextSize(0.045)
        c.SaveAs(tag+"_ftest.pdf")
        c.SaveAs(tag+"_ftest.png")


