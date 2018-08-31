#ifndef CHEBYSHEVBASIS
#define CHEBYSHEVBASIS

#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "RooFit.h"
#include "Riostream.h" 
#include "RooRealVar.h"
#include <map>

class chebyshevBasis : public TObject {
public:
    chebyshevBasis() {};
    chebyshevBasis( const char *name,
                    const char *title,
                    const int &_orderX,
                    const int &_orderY,
                    const TH2 &_binning );
    inline virtual ~chebyshevBasis (){};
    float getChebyshevVal(float xval, float yval, float thisOrderX, float thisOrderY, bool positive_x, bool positive_y) const; 
    void drawBasis(TFile* outfile);

protected:
    // variables
    Int_t orderX;
    Int_t orderY;
    Double_t xmin;
    Double_t xmax;
    Double_t ymin;
    Double_t ymax;
    float slope_x;
    float slope_y;
    std::map<std::string, TH2F> polyHists;  

    // functions
    std::pair<float, float> mapToChebyshev(float ix,float iy) const;
    float Eval2DChebyshev(float x, float y, int thisOrderX, int thisOrderY, bool positive_x, bool positive_y) const;
    TH2F Make2DChebyshev(const int& thisOrderX, const int& thisOrderY, const bool positive_x, const bool positive_y) const;
    

private:
    ClassDef(chebyshevBasis, 1)
};

#endif