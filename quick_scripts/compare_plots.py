import ROOT
from ROOT import *

gen = 'flat100kv1kBlind'

c = TFile.Open('../'+gen+'/postfitshapes.root')
m = TFile.Open('../unit_tests/distributions/flat100kv1k.root')
p = TFile.Open('../../2DAlphabet/'+gen+'/scaled_hists.root')
fd = TFile.Open('../'+gen+'/fitDiagnostics.root')
ub = TFile.Open('../'+gen+'/postfitshapes_full.root')


# Pre-fit comparisons
for thisSet in ['data_obs','signal']:
    for thisCat in ['pass','fail']:
        can = TCanvas('can','can',1600,700)
        can.Divide(3,1)

        c_prefit = c.Get(thisCat + '_prefit/'+thisSet)
        m_prefit = m.Get(thisSet+'_'+thisCat).Rebin2D(2,2)

        c_prefit.Sumw2()
        m_prefit.Sumw2()

        diff_prefit = c_prefit.Clone(thisSet+'_'+thisCat+'_prefit_diff')
        diff_prefit.Add(m_prefit,-1)


        can.cd(1)
        c_prefit.SetTitle('Combine Pre-fit Plotting - ' + thisSet + ', ' + thisCat)
        c_prefit.SetMinimum(0)
        c_prefit.Draw('lego')
        can.cd(2)
        m_prefit.SetTitle('Custom Pre-fit Plotting - ' + thisSet + ', ' + thisCat)
        m_prefit.SetMinimum(0)
        m_prefit.Draw('lego')
        can.cd(3)
        diff_prefit.SetTitle('Difference')
        diff_prefit.Draw('lego')

        can.Print('2DValidation_'+gen+'_prefit_'+thisSet+'_'+thisCat+'.pdf(','pdf')

        can2 = TCanvas('can2','can2',800,700)
        c_prefit_projy = c_prefit.ProjectionY(c_prefit.GetName()+'_py',0,-1,'e')
        c_prefit_projy.SetTitle('Combine Pre-fit Plotting Projection onto Y axis - ' + thisSet + ', ' + thisCat)
        c_prefit_projy.Draw('pe')
        can2.Print('2DValidation_'+gen+'_prefit_'+thisSet+'_'+thisCat+'.pdf','pdf')

        can3 = TCanvas('can3','can3',800,700)
        m_prefit_projy = m_prefit.ProjectionY(m_prefit.GetName()+'_py',0,-1,'e') 
        m_prefit_projy.SetTitle('Custom Pre-fit Plotting Projection onto Y axis - ' + thisSet + ', ' + thisCat)
        m_prefit_projy.Draw('pe')
        can3.Print('2DValidation_'+gen+'_prefit_'+thisSet+'_'+thisCat+'.pdf)','pdf')

        # blindCan = TCanvas('blindCan','blindCan',1600,700)
        # blindCan.Divide(3,1)
        # blindCan.cd(1)

        # unblinded_prefit = ub.Get(thisCat + '_prefit/'+thisSet)
        # unblinded_prefit.Draw('lego')
        # blindCan.cd(2)
        # c_prefit.Draw('lego')
        # blindCan.cd(3)
        # diff_blind = unblinded_prefit.Clone()
        # diff_blind.Add(c_prefit,-1)
        # diff_blind.Draw('lego')


        # blindCan.Print('2DValidation_'+gen+'_postfit_'+thisSet+'_'+thisCat+'.pdf)','pdf')


# Post-fit comparisons
for thisSet in ['qcd','signal']:
    for thisCat in ['pass','fail']:
        can = TCanvas('can','can',1600,700)
        can.Divide(3,1)

        c_postfit = c.Get(thisCat + '_postfit/'+thisSet)
        p_postfit = p.Get(thisSet+'_'+thisCat+'__xaxis_yaxis')

        c_postfit.Sumw2()
        p_postfit.Sumw2()

        diff_postfit = c_postfit.Clone(thisSet+'_'+thisCat+'_postfit_diff')
        diff_postfit.Divide(p_postfit)

        can.cd(1)
        c_postfit.SetTitle('Combine Post-fit Plotting - ' + thisSet + ', ' + thisCat)
        c_postfit.SetMinimum(0)
        c_postfit.Draw('lego')
        can.cd(2)
        p_postfit.SetTitle('Custom Post-fit Plotting - ' + thisSet + ', ' + thisCat)
        p_postfit.SetMinimum(0)
        p_postfit.Draw('lego')
        can.cd(3)
        diff_postfit.SetTitle('Difference')
        diff_postfit.Draw('lego')

        # if thisSet == 'signal':
        #     for i in range(1,13):
        #         for j in range(1,11):
        #             print diff_postfit.GetBinContent(i,j)

        can.Print('2DValidation_'+gen+'_postfit_'+thisSet+'_'+thisCat+'.pdf(','pdf')

        can2 = TCanvas('can2','can2',800,700)
        c_postfit_projy = c_postfit.ProjectionY(c_postfit.GetName()+'_py',0,-1,'e')
        c_postfit_projy.SetTitle('Combine Post-fit Plotting Projection onto Y axis - ' + thisSet + ', ' + thisCat)
        c_postfit_projy.Draw('pe')
        can2.Print('2DValidation_'+gen+'_postfit_'+thisSet+'_'+thisCat+'.pdf','pdf')

        can3 = TCanvas('can3','can3',800,700)
        p_postfit_projy = p_postfit.ProjectionY(p_postfit.GetName()+'_py',0,-1,'e') 
        p_postfit_projy.SetTitle('Custom Post-fit Plotting Projection onto Y axis - ' + thisSet + ', ' + thisCat)
        p_postfit_projy.Draw('pe')
        can3.Print('2DValidation_'+gen+'_postfit_'+thisSet+'_'+thisCat+'.pdf','pdf')


        blindCan = TCanvas('blindCan','blindCan',1600,1500)
        blindCan.Divide(2,2)
        blindCan.cd(1)

        unblinded_postfit = ub.Get(thisCat + '_postfit/'+thisSet)
        unblinded_postfit.SetTitle('Unblinded morph using blinded fit results')
        unblinded_postfit.SetMinimum(0)
        unblinded_postfit.Draw('lego')
        blindCan.cd(2)
        c_postfit.Draw('lego')
        blindCan.cd(3)
        diff_blind = unblinded_postfit.Clone()
        diff_blind.Add(c_postfit,-1)
        diff_blind.SetTitle('Unblinded minus Blinded Plots (sidebands should be 0)')
        diff_blind.Draw('lego')
        blindCan.cd(4)
        diff_blind_zeroed = diff_blind.Clone()
        # for i in range(1,13):
        for j in range(1,11):
            diff_blind_zeroed.SetBinContent(8,j,0)
        diff_blind_zeroed.SetTitle('Difference again but signal region bins set to 0')
        diff_blind_zeroed.Draw('lego')
        blindCan.Print('2DValidation_'+gen+'_postfit_'+thisSet+'_'+thisCat+'.pdf)','pdf')



        if thisSet == 'signal':
            can4 = TCanvas('can4','can4',1600,1700)
            can4.Divide(2,2)

            can4.cd(1)
            c_postfit.Draw('lego')
            can4.cd(2)
            c_scaled = c_postfit.Clone(c_postfit.GetName()+'_scaled')
            tree = fd.Get('tree_fit_sb')
            tree.GetEntry(0)
            c_scaled.Scale(abs(tree.r))
            c_scaled.Draw('lego')
            can4.cd(3)
            p_postfit.Draw('lego')
            can4.cd(4)
            diff_scale = c_scaled.Clone(p_postfit.GetName()+'_scale_diff')
            diff_scale.Divide(p_postfit)
            

            for i in range(1,13):
                for j in range(1,11):
                    if diff_scale.GetBinContent(i,j) == 0:
                        diff_scale.SetBinContent(i,j,1)
            diff_scale.Draw('surf')

            can4.Print('2DValidation_'+gen+'_scale_'+thisSet+'_'+thisCat+'.pdf','pdf')
