from TwoDAlphabet.helpers import roofit_form_to_TF1
from numpy.lib.function_base import piecewise
from ROOT import RooRealVar, RooFormulaVar, RooArgList, RooParametricHist2D, RooConstVar, TFormula
from binning import copy_hist_with_new_bins
import numpy as np

class Generic2D:
    '''Wraps various input distributions in a common type so that
    distributions can easily be manipulated and compared.
    '''
    def __init__(self,name,binning):
        self.name = name
        self.binning = binning
        self.nuisances = []
        self.binVars = {}
        self.obj = {}
        self.allVars = []

    def __add__(self,other):
        pass
    def __sub__(self,other):
        pass
    def __mul__(self,other):
        pass
    def __div__(self,other):
        pass

    def __radd__(self,other):
        return self+other
    def __rsub__(self,other):
        raise NotImplementedError
    def __rmul__(self,other):
        return self*other
    def __rdiv__(self,other):
        raise NotImplementedError

class ParametricFunction(Generic2D):
    def __init__(self,name,formula,constriants={}):
        '''[summary]

        Args:
            formula ([type]): Must reference by ordinal with @. Use "x" and "y" to represent
                the "x" and "y" axes of the space. All other terms are indexed starting at 0. Ex. "@0 + x*@1 +y*@2".
            constriants (dict, optional): Map of formula parameters to constraint information. Defaults to {} in which
                case the constraint will be flat, the starting value of the parameter will be 0 with a step size of 0.1,
                and the range of the parameter will be [-1000,1000].
        '''
        super().__init__(name,binning)
        self.funcVars = {}
        self.npars = self.getNparams()
        self.formula = self.replace_xy(formula)
        self.func = RooFormulaVar()

    def replace_xy(self,f):
        nCoeffs = self.npars
        f = f.replace('+x','+@'+str(nCoeffs)).replace('+y','+@'+str(nCoeffs+1))
        f = f.replace('*x','*@'+str(nCoeffs)).replace('*y','*@'+str(nCoeffs+1))
        f = f.replace('-x','-@'+str(nCoeffs)).replace('-y','-@'+str(nCoeffs+1))
        f = f.replace('/x','/@'+str(nCoeffs)).replace('/y','/@'+str(nCoeffs+1))
        f = f.replace('(x','(@'+str(nCoeffs)).replace('(y','(@'+str(nCoeffs+1))
        f = f.replace(' x',' @'+str(nCoeffs)).replace(' y',' @'+str(nCoeffs+1))
        return f

    def getNparams(self,f):
        return TFormula('t',roofit_form_to_TF1(f)).GetNpars()

class BinnedDistribution(Generic2D):
    def __init__(self,name,inhist,binning,constant=False):
        '''[summary]

        Args:
            inhist (TH2): Input 2D histogram to build RooParametricHist2D.
            binning (Binning): Binning object used to create LOW, SIG, HIGH regions along X axis.
            constant (bool, optional): If true, use RooConstVars for bins. Defaults to False and RooRealVars are used.
        '''
        super().__init__(name,binning)
        for cat in ['LOW','SIG','HIGH']:
            cat_name = name+'_'+cat
            cat_hist = copy_hist_with_new_bins(cat_name,'X',inhist,self.binning.xbinByCat[cat])
            bin_pars = RooArgList()
            for ybin in range(1,inhist.GetNbinsY()+1):
                for xbin in range(1,inhist.GetNbinsX()+1):
                    bin_name = '%s_%s_x%sy%s'%(name,cat,xbin,ybin)
                    if constant:
                        self.binVars[bin_name] = RooConstVar(bin_name, bin_name, inhist.GetBinContent(xbin,ybin))
                    else:
                        self.binVars[bin_name] = RooRealVar(bin_name, bin_name, inhist.GetBinContent(xbin,ybin))
                    self.allVars.append(self.binVars[bin_name])
                    bin_pars.add(self.binVars[bin_name])

            self.obj[cat] = RooParametricHist2D(cat_name, cat_name, self.binning.xVars[cat], self.binning.yVar, bin_pars, cat_hist)
                     
    def AddShapeTemplates(self,nuis_name,up_shape,down_shape,constraint="param"):
        nuisance_par = RooRealVar(nuis_name,nuis_name,0,-5,5)
        self.nuisances.append({'name':nuis_name, 'constraint':constraint, 'obj': nuisance_par})

        for cat in ['LOW','SIG','HIGH']:
            rph_name = self.obj[cat].GetName()
            cat_hist_up =   copy_hist_with_new_bins(up_shape.GetName()+'_'+cat,  'X', up_shape,   self.binning.xbinByCat[cat])
            cat_hist_down = copy_hist_with_new_bins(down_shape.GetName()+'_'+cat,'X', down_shape, self.binning.xbinByCat[cat])
            bin_pars = RooArgList()
            for ybin in range(1,cat_hist_up.GetNbinsY()+1):
                for xbin in range(1,cat_hist_up.GetNbinsX()+1):
                    bin_name = '%s_%s_x%sy%s'%(rph_name,nuis_name,xbin,ybin)
                    self.binVar[bin_name] = singleBinInterpCombine( # change to singleBinInterpQuad to change interpolation method
                                                bin_name, self.binVars[bin_name], nuisance_par,
                                                cat_hist_up.GetBinContent(xbin,ybin),
                                                cat_hist_down.GetBinContent(xbin,ybin)
                    )
                    bin_pars.add(self.binVar[bin_name])
            # replace existing RPH
            self.obj[cat] = RooParametricHist2D(rph_name, rph_name, self.binning.xVars[cat], self.binning.yVar, bin_pars, cat_hist_up)

class SmoothedDistribtuion(BinnedDistribution):
    def __init__(self,name,binning,inhist):
        super().__init__(name,binning)

# ------------------
def singleBinInterpQuad(name, nuis, binVar, upVal, downVal):
    nomVal = binVar.getValV()
    a,b,c = solve_quad([(-1,downVal),(0,nomVal),(1,upVal)])
    m_0, b_0 = linear_extrap((-1,downVal),[a,b,c])
    m_1, b_1 = linear_extrap((1,upVal),[a,b,c])
    pieces = [
        '(@0 < -1)*(%s*@0+%s)'%(m_0,b_0),
        '(@0 > -1 && @0 < 1)*(%s*@0**2+%s*@0+%s)'%(a,b,c),
        '(@0 > 1)*(%s*@0+%s)'%(m_1,b_1)]
    
    return RooFormulaVar(name, name, '@1*(({0})-{1})/{1}'.format('+'.join(pieces),nomVal), RooArgList(nuis,binVar))

def singleBinInterpCombine(name, nuis, binVar, upVal, downVal):
    nomVal = binVar.getValV()
    pieces = '(abs(@0) < 1)*(0.125*@0*(pow(@0,2)*(3*pow(@0,2) - 10) + 15))-(@0 < -1)+(@0 > 1)'
    alpha = '@1*( (@0/2) * (({up}-{down})+({up}+{down}-2*{nom})*({piecewise})) )/{nom}'.format(up=upVal,down=downVal,nom=nomVal,piecewise=pieces)
    return RooFormulaVar(name,name,alpha,RooArgList(nuis,binVar))

class quadratic:
    def __init__(self, params):
        self.a = params[0]
        self.b = params[1]
        self.c = params[2]
    def y(self,x):
        return self.a * x**2 + self.b * x + self.c
    def slope(self,x):
        return 2*self.a*x + self.b

class line:
    def __init__(self, slope, intercept):
        self.m = slope
        self.b = intercept
    def y(self,x):
        return self.m*x + self.b

def solve_quad(threepts):
    # A = [[x_1^2, x_1, 1],
    #      [x_2^2, x_2, 1],
    #      [x_3^2, x_3, 1]]
    # X = [[a],[b],[c]]
    # B = [[y_1], [y_2], [y_3]]
    A, B = [], []
    for pt in threepts:
        x, y = pt[0], pt[1]
        A.append([x**2, x, 1])
        B.append([y])
    A = np.array(A)
    B = np.array(B)
    out = [x[0] for x in np.linalg.solve(A,B)]
    return out[0], out[1], out[2]

def linear_extrap(point,quad_params):
    q = quadratic(quad_params)
    slope = q.slope(point[0])
    const = -1*slope*point[0]+point[1]
    return slope,const

print (solve_quad([(0,0),(-1,1),(1,1)]))

def InitQCDHist(data,bkgList):
    qcd = data.Clone(data.GetName().replace('data_obs','qcd'))
    for bkg in bkgList:
        qcd.Add(bkg,-1)
    return qcd
