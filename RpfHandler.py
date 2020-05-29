import ROOT
from ROOT import *

import header, os, subprocess, copy

class RpfHandler():
    def __init__ (self,fitDict,name,dummy_TH2=None,tag=''):
        self.fitDict = fitDict
        self.projDir = name+'/'+tag
        self.fitType = self.fitType()
        self.name = name
        self.tag = tag
        self.dummy_TH2 = dummy_TH2
        self.rpfVars = {}
        self.binVars = {} # RooFormulaVars for each bin for later evaluation
        self.allVars = []

        # Initialize rpfVars according to form of function
        if self.fitType == 'splitPoly':
            # Do some quick checks to make sure these are formatted correctly
            header.checkFitForm(self.fitDict['XPFORM'],self.fitDict['YPFORM'])

            nxparams = max([int(param[1:]) for param in self.fitDict.keys() if param.find('X') != -1 and param != 'XPFORM'])
            nyparams = max([int(param[1:]) for param in self.fitDict.keys() if param.find('Y') != -1 and param != 'YPFORM'])
            print 'Total number of x fit parameters is ' + str(nxparams)
            print 'Total number of y fit parameters is ' + str(nyparams)

            # Make and store RRVs for each param
            for this_var in {'X':nxparams, 'Y':nyparams}.keys():
                nparams = {'X':nxparams, 'Y':nyparams}[this_var]
                for ip in range(1,nparams+1):
                    input_param_vals = self.fitDict[this_var+str(ip)]
                    this_nom = input_param_vals['NOMINAL']

                    varname = self.fitType+this_var+'_'+str(ip)+'_'+self.name

                    # Set ranges if they're defined
                    if 'MIN' in input_param_vals.keys() and 'MAX' in input_param_vals.keys():
                        this_min = input_param_vals['MIN']
                        this_max = input_param_vals['MAX']
                        self.rpfVars[varname] = RooRealVar(varname,varname,this_nom,this_min,this_max)

                        if 'ERROR_UP' in input_param_vals.keys() and 'ERROR_DOWN' in input_param_vals.keys():
                            self.rpfVars[varname].setAsymError(input_param_vals['ERROR_DOWN'],input_param_vals['ERROR_UP'])
                        elif 'ERROR' in input_param_vals.keys():
                            self.rpfVars[varname].setError(input_param_vals['ERROR'])

                    else:
                        input_break = raw_input('WARNING: Upper and lower bounds on global fit parameter ' +this_var+str(ip) + ' not defined. Please type "N" if you would like to quit or enter if you would like to treat the parameter as constant')
                        if input_break == 'N':
                            quit()
                        else:
                            self.rpfVars[varname] = RooConstVar(varname,varname,this_nom)

                    # self.allVars.append(self.rpfVars[varname])

        elif self.fitType == 'fullPoly':
            # Polynomial Order
            polXO = 0
            polYO = 0
            for param_name in [key for key in self.fitDict.keys() if key != 'HELP' and key != 'PFORM']:
                # Assuming poly order is a single digit (pretty reasonable I think...)
                tempXorder = int(param_name[param_name.find('X')+1])
                tempYorder = int(param_name[param_name.find('Y')+1])
                if tempXorder > polXO:
                    polXO = tempXorder
                if tempYorder > polYO:
                    polYO = tempYorder

            self.polXO = polXO
            self.polYO = polYO

            for yi in range(polYO+1):
                for xi in range(polXO+1):

                    input_param_vals = self.fitDict['X'+str(xi)+'Y'+str(yi)]
                    this_nom = input_param_vals['NOMINAL']
                    
                    varname = self.fitType+'X'+str(xi)+'Y'+str(yi)+'_'+self.name

                    if 'MIN' in input_param_vals.keys() and 'MAX' in input_param_vals.keys():
                        this_min = input_param_vals['MIN']
                        this_max = input_param_vals['MAX']
                        self.rpfVars[varname] = RooRealVar(varname,varname,this_nom,this_min,this_max)

                        if 'ERROR_UP' in input_param_vals.keys() and 'ERROR_DOWN' in input_param_vals.keys():
                            self.rpfVars[varname].setAsymError(input_param_vals['ERROR_DOWN'],input_param_vals['ERROR_UP'])
                        elif 'ERROR' in input_param_vals.keys():
                            self.rpfVars[varname].setError(input_param_vals['ERROR'])

                    else:
                        input_break = raw_input('WARNING: Upper and lower bounds on fit parameter ' + 'x'+str(xi)+'y'+str(yi) + ' not defined. Please type "N" if you would like to quit or enter if you would like to treat the parameter as constant')
                        if input_break == 'N':
                            quit()
                        else:
                            self.rpfVars[varname] = RooConstVar(varname,varname,this_nom)
                    # self.allVars.append(self.rpfVars[varname])

        elif self.fitType == 'generic':
            for c in sorted(self.fitDict.keys()):
                input_param_vals = self.fitDict[c]
                if c.isdigit():
                    this_nom = input_param_vals['NOMINAL']
                    if 'MIN' in input_param_vals.keys() and 'MAX' in input_param_vals.keys():
                        this_min = input_param_vals['MIN']
                        this_max = input_param_vals['MAX']
                        varname = self.fitType+c+'_'+self.name
                        self.rpfVars[varname] = RooRealVar(varname, varname, this_nom, this_min, this_max)

                        if 'ERROR_UP' in input_param_vals.keys() and 'ERROR_DOWN' in input_param_vals.keys():
                            self.rpfVars[varname].setAsymError(input_param_vals['ERROR_DOWN'],input_param_vals['ERROR_UP'])
                        elif 'ERROR' in self.fitDict[c].keys():
                            self.rpfVars[varname].setError(input_param_vals['ERROR'])

                    else:
                        input_break = raw_input('WARNING: Upper and lower bounds on fit parameter ' + c + ' not defined. Please type "N" if you would like to quit or enter if you would like to treat the parameter as constant')
                        if input_break == 'N':
                            quit()
                        else:
                            self.rpfVars[varname] = RooConstVar(varname,varname,this_nom)
                    # self.allVars.append(self.rpfVars[varname])

            
        elif self.fitType == 'cheb':
            # Make a dummy TH2 for binning to be used by the chebyshevBasis class
            self.chebBasis = chebyshevBasis(name,'chebBasis',self.fitDict['CHEBYSHEV']['XORDER'],self.fitDict['CHEBYSHEV']['YORDER'],dummy_TH2)
            if os.path.isdir(self.projDir+'/basis_plots'): header.executeCmd('rm -rf %s/basis_plots'%self.projDir)
            header.executeCmd('mkdir %s/basis_plots'%self.projDir)

            self.chebBasis.drawBasis(self.projDir+'/basis_plots')

            # This class handles the RRV management itself so we just keep track of the class instance
            cheb_coeffs = self.chebBasis.getCoeffs()
            i = 0
            while cheb_coeffs.at(i):
                c = cheb_coeffs.at(i)
                self.rpfVars[c.GetName()] = c
                i += 1

    def fitType(self):
        if 'XPFORM' in self.fitDict.keys() and 'YPFORM' in self.fitDict.keys():
            return 'splitPoly'
        elif 'PFORM' in self.fitDict.keys():
            return 'fullPoly'
        elif 'CHEBYSHEV' in self.fitDict.keys():
            return 'cheb'
        elif 'FORM' in self.fitDict.keys():
            return 'generic'
        else:
            print 'ERROR: Fit type failed. Check that your FIT section is defined correctly. Quitting...'
            quit()

    def getRooFunctionForm(self):
        if self.fitType == 'splitPoly':
            return '('+self.fitDict['XPFORM'] +')*('+self.fitDict['YPFORM']+')'
        elif self.fitType == 'fullPoly':
            return self.fitDict['PFORM']
        elif self.fitType == 'generic':
            params = [int(param) for param in self.fitDict.keys() if param != 'FORM' and param != 'HELP']
            if len(params) > 0:
                nCoeffs = max([int(param) for param in self.fitDict.keys() if param != 'FORM' and param != 'HELP'])+1
            else: nCoeffs = 0
            
            gFormula = self.fitDict['FORM'].replace('+x','+@'+str(nCoeffs)).replace('+y','+@'+str(nCoeffs+1))
            gFormula = gFormula.replace('*x','*@'+str(nCoeffs)).replace('*y','*@'+str(nCoeffs+1))
            gFormula = gFormula.replace('-x','-@'+str(nCoeffs)).replace('-y','-@'+str(nCoeffs+1))
            gFormula = gFormula.replace('/x','/@'+str(nCoeffs)).replace('/y','/@'+str(nCoeffs+1))
            gFormula = gFormula.replace('(x','(@'+str(nCoeffs)).replace('(y','(@'+str(nCoeffs+1))
            gFormula = gFormula.replace(' x',' @'+str(nCoeffs)).replace(' y',' @'+str(nCoeffs+1))
            return gFormula

        elif self.fitType == 'cheb':
            return False    # No string form for this. Dependencies already built into chebyshevBasis class

    def evalRpf(self,xConst,yConst,xbin,ybin):
        # Get category name for naming
        c = xConst.GetName().split('_')[2]

        full_formula_name = 'formula_'+c+'_bin_'+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name

        if self.fitType == 'splitPoly':
            # Make var lists
            x_param_list = RooArgList(xConst)
            y_param_list = RooArgList(yConst)
            for v in self.getRpfVarNames():
                if 'splitPolyX' in v:
                    x_param_list.add(self.rpfVars[v])
                elif 'splitPolyY' in v:
                    y_param_list.add(self.rpfVars[v])

            # Make formulas
            x_formula_name = 'xPol_'+c+'_bin_'+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name
            x_formula_var = RooFormulaVar(x_formula_name, x_formula_name, self.fitDict['XPFORM'].replace('x','@0'), x_param_list)
            self.allVars.extend([x_formula_var, x_param_list])

            y_formula_name = 'yPol_'+c+'_bin_'+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name
            y_formula_var = RooFormulaVar(y_formula_name, y_formula_name, self.fitDict['YPFORM'].replace('y','@0'), y_param_list)
            self.allVars.extend([y_formula_var, y_param_list])

            rpf_val = RooFormulaVar(full_formula_name, full_formula_name, "max(0.000001,@0*@1)", RooArgList(x_formula_var, y_formula_var))


        elif self.fitType == 'fullPoly':
            x_poly_list = RooArgList()
            for y_coeff in range(self.polYO+1):
                x_coeff_list = RooArgList()

                # Get each x coefficient for this y
                for x_coeff in range(self.polXO+1):
                    v = self.fitType+'X'+str(x_coeff)+'Y'+str(y_coeff)+'_'+self.name                  
                    x_coeff_list.add(self.rpfVars[v])

                # Make the polynomial and save it to the list of x polynomials
                this_x_polyvar_label = "xPol_y"+str(y_coeff)+'_'+c+"_bin_"+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name
                x_poly_var = RooPolyVar(this_x_polyvar_label,this_x_polyvar_label,xConst,x_coeff_list)
                x_poly_list.add(x_poly_var)
                self.allVars.append(x_poly_var)

            # Now make a polynomial out of the x polynomials
            this_x_polyvar_label = "yPol_y"+str(y_coeff)+'_'+c+"_bin_"+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name
            rpf_val_poly = RooPolyVar(this_x_polyvar_label,this_x_polyvar_label,yConst,x_poly_list)
            self.allVars.append(rpf_val_poly)
            rpf_val = RooFormulaVar(full_formula_name,full_formula_name,"max(0.000001,@0)",RooArgList(rpf_val_poly))

        elif self.fitType == 'generic':
            formula_list = RooArgList()
            for c in self.getRpfVarNames():
                formula_list.add(self.rpfVars[c])

            generic_formula = self.getRooFunctionForm()
            formula_list.add(xConst)
            formula_list.add(yConst)
            rpf_val = RooFormulaVar(full_formula_name,full_formula_name,"max(0.000001,"+generic_formula+")",formula_list)

        elif self.fitType == 'cheb':
            # chebBasis.getBinVal() returns a RooAddition with the proper construction to be the sum of the shapes with floating coefficients
            cheb_val = self.chebBasis.getBinVal(xConst.getValV(), yConst.getValV())
            rpf_val = RooFormulaVar(full_formula_name,full_formula_name,'(1+@0)/2',RooArgList(cheb_val)) # translate the [-1,1] range to [0,1] to guarantee positivity
            self.allVars.append(cheb_val)

        self.storeRpfBin(rpf_val,full_formula_name)         

        return rpf_val

    def storeRpfBin(self,rpf_val,name):
        self.binVars[name] = rpf_val

    def setRpfParam(self,varname,value):
        self.rpfVars[varname].setVal(value)

    def getRpfBinVal(self,c,xbin,ybin):
        formula_name = 'formula_'+c+'_bin_'+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name
        return self.binVars[formula_name].getValV()

    def getRpfBinRRV(self,c,xbin,ybin):
        formula_name = 'formula_'+c+'_bin_'+str(int(xbin))+"-"+str(int(ybin))+'_'+self.name
        return self.binVars[formula_name]

    def getRpfVarNames(self):
        return sorted(self.rpfVars.keys())

    # Removes everything but basic values for storage
    def getReducedCopy(self):
        newcopy = RpfHandler(self.fitDict,self.name,self.dummy_TH2,self.tag)
        newcopy.rpfVars = copy.deepcopy(self.rpfVars)
        return newcopy
