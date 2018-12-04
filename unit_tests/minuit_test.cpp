#include <stdio.h>
#include <iostream>

#include "TFitter.h"
#include "TMath.h"


// fcn passes back f = - 2*ln(L), (where L is the likelihood of the event)
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    double lnL = 0.0;
    double failbin = 10000;
    double passbin = 1000;
    double rpf = par[0];

    double passest = failbin*rpf;

    lnL = passest - passbin*TMath::Log(passest) + passbin*TMath::Log(passbin)-passbin; //TMath::Log(TMath::Factorial(passbin));//
    // std::cout << "TMath::Log(passest) = " << TMath::Log(passest) << std::endl;
    // std::cout << "TMath::Factorial(passbin) = " << TMath::Factorial(passbin) << std::endl;
    std::cout << lnL << std::endl;
    std::cout << "\t 10000 * " << rpf << " = " << passest << std::endl;

    f = 2.0 * lnL;

}


int minuit_test(){
    
    TVirtualFitter * minuit = TVirtualFitter::Fitter(0,1);
    minuit->SetFCN(fcn);
    minuit->SetParameter(0,"rpf", 0.11, 0.01, 0, 0);
    Double_t arglist[10];
    arglist[0] = 10000.;
    minuit->ExecuteCommand("MIGRAD", arglist,0);
    Double_t up = 1.0;
    minuit->SetErrorDef(up);
    arglist[0] = 0.;
    minuit->ExecuteCommand("MINOS", arglist, 0);

    return 1;
}


