import ROOT
from ROOT import *

c = TFile.Open('../flat1M/postfitplots.root')
m = TFile.Open('../unit_tests/distributions/flat1M.root')
p = TFile.Open('../flat1M/scaled_hists.root')


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
        c_prefit.Draw('lego')
        can.cd(2)
        m_prefit.SetTitle('Custom Pre-fit Plotting - ' + thisSet + ', ' + thisCat)
        m_prefit.Draw('lego')
        can.cd(3)
        diff_prefit.SetTitle('Difference')
        diff_prefit.Draw('surf')

        can.Print('../flat1M/2Dvalidation/prefit_'+thisSet+'_'+thisCat+'.pdf(','pdf')

        can2 = TCanvas('can2','can2',800,700)
        c_prefit_projy = c_prefit.ProjectionY('e')
        c_prefit_projy.SetTitle('Combine Pre-fit Plotting Projection onto Y axis - ' + thisSet + ', ' + thisCat)
        c_prefit_projy.Draw('pe')
        can2.Print('../flat1M/2Dvalidation/prefit_'+thisSet+'_'+thisCat+'.pdf','pdf')

        can3 = TCanvas('can3','can3',800,700)
        m_prefit_projy = m_prefit.ProjectionY('e') 
        m_prefit_projy.SetTitle('Custom Pre-fit Plotting Projection onto Y axis - ' + thisSet + ', ' + thisCat)
        m_prefit_projy.Draw('pe')
        can3.Print('../flat1M/2Dvalidation/prefit_'+thisSet+'_'+thisCat+'.pdf)','pdf')


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
        c_postfit.Draw('lego')
        can.cd(2)
        p_postfit.SetTitle('Custom Post-fit Plotting - ' + thisSet + ', ' + thisCat)
        p_postfit.Draw('lego')
        can.cd(3)
        diff_postfit.SetTitle('Difference')
        diff_postfit.Draw('lego')

        if thisSet == 'signal':
            for i in range(1,13):
                for j in range(1,11):
                    print diff_postfit.GetBinContent(i,j)

        can.Print('../flat1M/2Dvalidation/postfit_'+thisSet+'_'+thisCat+'.pdf(','pdf')

        can2 = TCanvas('can2','can2',800,700)
        c_postfit_projy = c_postfit.ProjectionY('e')
        c_postfit_projy.SetTitle('Combine Post-fit Plotting Projection onto Y axis - ' + thisSet + ', ' + thisCat)
        c_postfit_projy.Draw('pe')
        can2.Print('../flat1M/2Dvalidation/postfit_'+thisSet+'_'+thisCat+'.pdf','pdf')

        can3 = TCanvas('can3','can3',800,700)
        p_postfit_projy = p_postfit.ProjectionY('e') 
        p_postfit_projy.SetTitle('Custom Post-fit Plotting Projection onto Y axis - ' + thisSet + ', ' + thisCat)
        p_postfit_projy.Draw('pe')
        can3.Print('../flat1M/2Dvalidation/postfit_'+thisSet+'_'+thisCat+'.pdf)','pdf')


