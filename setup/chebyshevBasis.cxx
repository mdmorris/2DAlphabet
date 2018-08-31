#include "HiggsAnalysis/CombinedLimit/interface/chebyshevBasis.h"
#include <map>

ClassImp(chebyshevBasis)

chebyshevBasis::chebyshevBasis( const char *name,
                                const char *title,
                                const Int_t &_orderX,
                                const Int_t &_orderY,
                                const TH2 &_binning ) :
    orderX(_orderX), orderY(_orderY)
{
    // Make TH2Fs with basis Chebyshevs and RooRealVars
    std::string polyKey;
    for (Int_t oX=0; oX<=orderX; ++oX){
        for (Int_t oY=0; oY<=orderY; ++oY) {
            if (oX == 0 and oY == 0){
                polyKey = "x"+std::to_string(oX)+"y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,true,true);
            } else if (oX == 0) {
                polyKey = "x"+std::to_string(oX)+"+y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,true,true);
                polyKey = "x"+std::to_string(oX)+"-y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,true,false);
            } else if (oY == 0) {
                polyKey = "+x"+std::to_string(oX)+"y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,true,true);
                polyKey = "-x"+std::to_string(oX)+"y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,false,true);
            } else {
                // +x, +y
                polyKey = "+x"+std::to_string(oX)+"+y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,true,true);
                // +x, -y
                polyKey = "+x"+std::to_string(oX)+"-y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,true,false);
                // -x, +y
                polyKey = "-x"+std::to_string(oX)+"+y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,false,true);
                // -x, -y
                polyKey = "-x"+std::to_string(oX)+"-y"+std::to_string(oY);
                polyHists[polyKey] = Make2DChebyshev(oX,oY,false,false);  
            }
        }
    }


    // Make the var map variables
    xmin = _binning.GetXaxis()->GetXmin();
    xmax = _binning.GetXaxis()->GetXmax();
    ymin = _binning.GetYaxis()->GetXmin();
    ymax = _binning.GetYaxis()->GetXmax();
    slope_x = 2/(xmax-xmin);
    slope_y = 2/(ymax-ymin);

}

float chebyshevBasis::Eval2DChebyshev(float x, float y, int thisOrderX, int thisOrderY,bool positive_x, bool positive_y) const{
    // Initialize polynomials
    float Tx_n2 = 0; // T_n-2
    float Tx_n1 = 0; // T_n-1
    float Tx_n = 0;  // T_n
    float Ty_n2 = 0; // T_n-2
    float Ty_n1 = 0; // T_n-1
    float Ty_n = 0;  // T_n
    float Tx = 0;
    float Ty = 0;

    // Construct polynomial in x
    for(Int_t i=0; i<=thisOrderX; ++i){
        // T_0 == 1
        if (i == 0) {
            Tx_n = 1;
        }
        // T_1 = x
        else if (i == 1) {
            Tx_n1 = 1;
            Tx_n = x;
        }
        // T_n = 2*x*T_n-1 - T_n-2
        else {
            Tx_n2 = Tx_n1;
            Tx_n1 = Tx_n;
            Tx_n = 2*x*Tx_n1-Tx_n2;
        }

    }
    // Construct polynomial in y
    for (Int_t j=0; j<=thisOrderY; ++j){
        // T_0 == 1
        if (j == 0) {
            Ty_n = 1;
        }
        // T_1 = x
        else if (j == 1) {
            Ty_n1 = 1;
            Ty_n = y;
        }
        // T_n = 2*x*T_n-1 - T_n-2
        else {
            Ty_n2 = Ty_n1;
            Ty_n1 = Ty_n;
            Ty_n = 2*y*Ty_n1-Ty_n2;
        }
    }

    // std::cout << "X: " + std::to_string(Tx_n) + ", " + std::to_string(Tx_n1) + ", " + std::to_string(Tx_n2) << std::endl;
    // std::cout << "Y: " + std::to_string(Ty_n) + ", " + std::to_string(Ty_n1) + ", " + std::to_string(Ty_n2) << std::endl;

    // Only want to shift if order is > 0
    if (thisOrderX == 0) {
        Tx = Tx_n;
    } else if (positive_x) {
        Tx = (Tx_n+1.0)/2.0;
    } else if (!positive_x) {
        Tx = (Tx_n-1.0)/-2.0;
    }
    
    if (thisOrderY == 0) {
        Ty = Ty_n;
    } else if (positive_y) {
        Ty = (Ty_n+1.0)/2.0;
    } else if (!positive_y) {
        Ty = (Ty_n-1.0)/-2.0;
    }
    // std::cout << "oX = " << thisOrderX << ", oY = " << thisOrderY << std::endl;
    // std::cout << "Tx = " << Tx << ", Ty = " << Ty << std::endl;

    return Tx*Ty;

}

TH2F chebyshevBasis::Make2DChebyshev(const Int_t& thisOrderX, const Int_t& thisOrderY, const bool positive_x, const bool positive_y) const {
    // What gets made is a little tricky since we don't want duplicates or 'zero' distributions
    // Since T0 = 1, doing (T0+/-1) makes no sense so we'll keep it at T0
    // That means though that if we have Tx0 and Ty0, we only need one polynomial
    // and for one of either Tx0 or Ty0 (say Tx0), we only need two polynomials (plus and minus for TyN)

    std::string cheb_name;
    if (thisOrderX == 0 and thisOrderY == 0){
        cheb_name = "cheb_Tx"+std::to_string(thisOrderX)+"_Ty"+std::to_string(thisOrderY);
    } else if (thisOrderX == 0) {
        if (positive_y) {
                cheb_name = "cheb_Tx"+std::to_string(thisOrderX)+"_pTy"+std::to_string(thisOrderY);
            } else {
                cheb_name = "cheb_Tx"+std::to_string(thisOrderX)+"_nTy"+std::to_string(thisOrderY);
            }
    } else if (thisOrderY == 0) {
        if (positive_x) {
                cheb_name = "cheb_pTx"+std::to_string(thisOrderX)+"_Ty"+std::to_string(thisOrderY);
            } else {
                cheb_name = "cheb_nTx"+std::to_string(thisOrderX)+"_Ty"+std::to_string(thisOrderY);
            }
    }
    else{
        if (positive_x){
            if (positive_y) {
                cheb_name = "cheb_pTx"+std::to_string(thisOrderX)+"_pTy"+std::to_string(thisOrderY);
            } else {
                cheb_name = "cheb_pTx"+std::to_string(thisOrderX)+"_nTy"+std::to_string(thisOrderY);
            }
        } else {
            if (positive_y) {
                cheb_name = "cheb_nTx"+std::to_string(thisOrderX)+"_pTy"+std::to_string(thisOrderY);
            } else {
                cheb_name = "cheb_nTx"+std::to_string(thisOrderX)+"_nTy"+std::to_string(thisOrderY);
            }
        }
    }
    const char *cheb_name_cstr = cheb_name.c_str();

    // std::cout << std::to_string(orderX) + ", "+std::to_string(orderY) << std::endl;
    // Initialize the TH2F
    TH2F* cheb2D = new TH2F(cheb_name_cstr, cheb_name_cstr, 100, -1.0, 1.0, 100, -1.0, 1.0);

    // Loop over TH2F bins and evaluate
    for (Int_t xbin=1; xbin<=cheb2D->GetNbinsX(); xbin++){
        for (Int_t ybin=1; ybin<=cheb2D->GetNbinsY(); ybin++){
            float xcenter = cheb2D->GetXaxis()->GetBinCenter(xbin);
            float ycenter = cheb2D->GetYaxis()->GetBinCenter(ybin);
            float val = Eval2DChebyshev(xcenter,ycenter,thisOrderX,thisOrderY,positive_x,positive_y);
           
            // std::cout << "x "+std::to_string(xcenter) + ",y "+std::to_string(ycenter)+": "+std::to_string(val) << std::endl;

            cheb2D->SetBinContent(xbin,ybin,val);
        }
    }

    return *cheb2D;
}

std::pair<float, float> chebyshevBasis::mapToChebyshev(float ix,float iy) const {
    
    float xp = slope_x*(ix-xmax)+1;  // Map x
    float yp = slope_y*(iy-ymax)+1;  // Map y
    
    return std::make_pair(xp,yp);
}

float chebyshevBasis::getChebyshevVal(float xval, float yval, float thisOrderX, float thisOrderY, bool positive_x, bool positive_y) const {
    // Get Histogram of interest
    std::string polyToGet;
    if (thisOrderX == 0 && thisOrderY == 0) {
        polyToGet = "x"+std::to_string(static_cast<int>(thisOrderX))+"y"+std::to_string(static_cast<int>(thisOrderY));
    } else if (thisOrderX == 0) {
        if (positive_y){
            polyToGet = "x"+std::to_string(static_cast<int>(thisOrderX))+"+y"+std::to_string(static_cast<int>(thisOrderY));
        } else {
            polyToGet = "x"+std::to_string(static_cast<int>(thisOrderX))+"-y"+std::to_string(static_cast<int>(thisOrderY));
        }
    } else if (thisOrderY == 0) {
        if (positive_x){
            polyToGet = "+x"+std::to_string(static_cast<int>(thisOrderX))+"y"+std::to_string(static_cast<int>(thisOrderY));
        } else {
            polyToGet = "-x"+std::to_string(static_cast<int>(thisOrderX))+"y"+std::to_string(static_cast<int>(thisOrderY));
        }
    } else {
        if (positive_x){
            if (positive_y){
                polyToGet = "+x"+std::to_string(static_cast<int>(thisOrderX))+"+y"+std::to_string(static_cast<int>(thisOrderY));
            } else {
                polyToGet = "+x"+std::to_string(static_cast<int>(thisOrderX))+"-y"+std::to_string(static_cast<int>(thisOrderY));
            }
        } else {
            if (positive_y){
                polyToGet = "-x"+std::to_string(static_cast<int>(thisOrderX))+"+y"+std::to_string(static_cast<int>(thisOrderY));
            } else {
                polyToGet = "-x"+std::to_string(static_cast<int>(thisOrderX))+"-y"+std::to_string(static_cast<int>(thisOrderY));
            }
        }
    }
    TH2F thisHist = polyHists.at(polyToGet);

    // Map "real" axis values to [-1,1] range
    std::pair<float, float> mappedvals = mapToChebyshev(xval,yval);

    // Find the bins corrsponding to those values
    Int_t thisXbin = thisHist.GetXaxis()->FindBin(mappedvals.first);
    Int_t thisYbin = thisHist.GetYaxis()->FindBin(mappedvals.second);

    // Return the bin content
    return thisHist.GetBinContent(thisXbin,thisYbin);
}

void chebyshevBasis::drawBasis(TFile* outfile) {

    TCanvas* basisCan;
    std::string canName;

    outfile->cd();
    
    for (Int_t oX=0; oX<=orderX; ++oX){
        for (Int_t oY=0; oY<=orderY; ++oY) {
            if (oX == 0 and oY == 0){
                basisCan = new TCanvas("basisCan","basisCan",800,700);
                polyHists.at("x"+std::to_string(oX)+"y"+std::to_string(oY)).Draw("surf");
                polyHists.at("x"+std::to_string(oX)+"y"+std::to_string(oY)).Write();

            } else if (oX == 0) {
                basisCan = new TCanvas("basisCan","basisCan",1600,700);
                basisCan->Divide(2,1);
                basisCan->cd(1);
                polyHists.at("x"+std::to_string(oX)+"+y"+std::to_string(oY)).Draw("surf");
                polyHists.at("x"+std::to_string(oX)+"+y"+std::to_string(oY)).Write();
                basisCan->cd(2);
                polyHists.at("x"+std::to_string(oX)+"-y"+std::to_string(oY)).Draw("surf");
                polyHists.at("x"+std::to_string(oX)+"-y"+std::to_string(oY)).Write();


            } else if (oY == 0) {
                basisCan = new TCanvas("basisCan","basisCan",1600,700);
                basisCan->Divide(2,1);
                basisCan->cd(1);
                polyHists.at("+x"+std::to_string(oX)+"y"+std::to_string(oY)).Draw("surf");
                polyHists.at("+x"+std::to_string(oX)+"y"+std::to_string(oY)).Write();
                basisCan->cd(2);
                polyHists.at("-x"+std::to_string(oX)+"y"+std::to_string(oY)).Draw("surf");
                polyHists.at("-x"+std::to_string(oX)+"y"+std::to_string(oY)).Write();
            
            } else {
                basisCan = new TCanvas("basisCan","basisCan",800,700);
                basisCan->Divide(2,2);
                basisCan->cd(1);
                polyHists.at("+x"+std::to_string(oX)+"+y"+std::to_string(oY)).Draw("surf");
                polyHists.at("+x"+std::to_string(oX)+"+y"+std::to_string(oY)).Write();
                basisCan->cd(2);
                polyHists.at("-x"+std::to_string(oX)+"+y"+std::to_string(oY)).Draw("surf");
                polyHists.at("-x"+std::to_string(oX)+"+y"+std::to_string(oY)).Write();
                basisCan->cd(3);
                polyHists.at("+x"+std::to_string(oX)+"-y"+std::to_string(oY)).Draw("surf");
                polyHists.at("+x"+std::to_string(oX)+"-y"+std::to_string(oY)).Write();
                basisCan->cd(4);
                polyHists.at("-x"+std::to_string(oX)+"-y"+std::to_string(oY)).Draw("surf");
                polyHists.at("-x"+std::to_string(oX)+"-y"+std::to_string(oY)).Write();
            }
            canName = "basis_plots/x"+std::to_string(oX)+"y"+std::to_string(oY)+".png";
            basisCan->Print(canName.c_str(),"png");
        }
    }
}