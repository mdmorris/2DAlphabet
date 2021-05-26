#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

#include "TObject.h"
#include "TSystem.h"
#include "TFitter.h"
#include "TMath.h"
#include "TFile.h"
#include "TH2F.h"
#include "TObject.h"

TFile* infile = TFile::Open("distributions/lin100kv10k.root");
TH2F* data_failHist = (TH2F*)infile->Get("data_obs_fail");
TH2F* data_passHist = (TH2F*)infile->Get("data_obs_pass");

// TH2F *rpf = data_passHist->Clone('rpf')
// rpf->Divide(data_failHist)

int orderX = 1;
int orderY = 1;

double passbin;
double failbin;
double passest;
double lnL;
double rpf;

int nxbins;
int nybins;
float xval;
float yval;
float xunit;
float yunit;

// fcn passes back f = - 2*ln(L), (where L is the likelihood of the event)
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
    lnL = 0.0;
    // double failbin = 10000;
    // double passbin = 1000;
    // double rpf = par[0];
    nxbins = data_failHist->GetNbinsX();
    nybins = data_failHist->GetNbinsY();

    for (int xbinf = 1; xbinf <= nxbins; xbinf++) {
        for (int ybinf = 1; ybinf <= nybins; ybinf++) {
            // std::cout << xbinf << ", " << ybinf << std::endl;
            xval = data_failHist->GetXaxis()->GetBinCenter(xbinf);
            yval = data_failHist->GetYaxis()->GetBinCenter(ybinf);

            // xunit = (xval-data_failHist->GetXaxis()->GetXmin())/(data_failHist->GetXaxis()->GetXmax()-data_failHist->GetXaxis()->GetXmin());
            // yunit = (yval-data_failHist->GetYaxis()->GetXmin())/(data_failHist->GetYaxis()->GetXmax()-data_failHist->GetYaxis()->GetXmin());

            rpf = 0.0;
            for (int oy = 0; oy <= orderY; oy++) {
                for (int ox = 0; ox <= orderX; ox++) {
                    rpf += par[(orderX+1)*oy+ox]*pow(xval,ox)*pow(yval,oy);
                }
            }

            failbin = data_failHist->GetBinContent(xbinf,ybinf);
            passbin = data_passHist->GetBinContent(xbinf,ybinf);
            passest = failbin*rpf;

            if (rpf <= 0) passest = 0.000000001;

            if (passbin <= 0) passbin = 0.00000001;

            // Ran an actual test of this. Get infs at passbin > 170. But % difference between Sterling's approximation and true value is <1% for passbin > 90
            if (passbin < 90.0) {
                lnL += passest - passbin*TMath::Log(passest) + TMath::Log(TMath::Factorial(passbin));
            } else {
                lnL += passest - passbin*TMath::Log(passest) + passbin*TMath::Log(passbin)-passbin;
            }
            
        }
    }
    
    // std::cout << "TMath::Log(passest) = " << TMath::Log(passest) << std::endl;
    // std::cout << "TMath::Factorial(passbin) = " << TMath::Factorial(passbin) << std::endl;
    // std::cout << lnL << std::endl;
    // std::cout << "\t 10000 * " << rpf << " = " << passest << std::endl;

    f = 2.0 * lnL;

}


int minuit_test(){

    TVirtualFitter * minuit = TVirtualFitter::Fitter(0,1);
    minuit->SetFCN(fcn);

    string sx;
    string sy;
    string name;
    for (int oy = 0; oy <= orderY; oy++) {
        for (int ox = 0; ox <= orderX; ox++) {
            sx = std::to_string(ox);
            sy = std::to_string(oy);
            name = "polX"+sx+"Y"+sy;
            std::cout<<name<<" "<<(orderX+1)*oy+ox<<std::endl;
            minuit->SetParameter((orderX+1)*oy+ox, name.c_str(), 0.01, 0.005, 0, 1.0);
        }
    }

    
    Double_t arglist[10];
    arglist[0] = 10000.;
    minuit->ExecuteCommand("MIGRAD", arglist,0);
    Double_t up = 1.0;
    minuit->SetErrorDef(up);
    // arglist[0] = 0.;
    // minuit->ExecuteCommand("MINOS", arglist, 0);

    
    //"([0]+[1]*x+[2]*x**2)+([3]+[4]*x+[5]*x**2)*y+([6]+[7]*x+[8]*x**2)*y**2"
    //"[0]+[1]*x+[2]*y+[3]*x*y"
    //"[0]"


    TF2 *newrpf = new TF2("newrpf", "[0]+[1]*x+[2]*y+[3]*x*y", data_failHist->GetXaxis()->GetXmin(), data_failHist->GetXaxis()->GetXmax(), data_failHist->GetYaxis()->GetXmin(), data_failHist->GetYaxis()->GetXmax()); 
    for (int i = 0; i < 4; i++){
        newrpf->SetParameter(i,minuit->GetParameter(i));
        newrpf->SetParError(i,minuit->GetParError(i));
    }

    TH2F *final_est = (TH2F*)data_failHist->Clone("final_est");

    nxbins = data_failHist->GetNbinsX();
    nybins = data_failHist->GetNbinsY();
    float rpf_val;
    float xval;
    float yval;

    for (int xbinf = 1; xbinf <= nxbins; xbinf++) {
        for (int ybinf = 1; ybinf <= nybins; ybinf++) {
            xval = data_failHist->GetXaxis()->GetBinCenter(xbinf);
            yval = data_failHist->GetYaxis()->GetBinCenter(ybinf);
            rpf_val = newrpf->Eval(xval,yval);
            final_est->SetBinContent(xbinf,ybinf, final_est->GetBinContent(xbinf,ybinf)*rpf_val);
        }
    }

    // final_est->Multiply(newrpf);
    
    TFile* out = TFile::Open("result_minuit_lin.root","RECREATE");
    out->cd();
    data_failHist->Write();
    data_passHist->Write();
    newrpf->Write();
    final_est->Write();


    TH2F *diff = (TH2F*)final_est->Clone("diff");
    diff->Add(data_passHist,-1);
    diff->Write();

    return 1;

}


