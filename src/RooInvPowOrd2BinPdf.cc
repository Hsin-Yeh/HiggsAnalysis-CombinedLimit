//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <math.h>

#include "../interface/RooInvPowOrd2BinPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"

using namespace std;
using namespace RooFit;

ClassImp(RooInvPowOrd2BinPdf)
//---------------------------------------------------------------------------
RooInvPowOrd2BinPdf::RooInvPowOrd2BinPdf(const char *name, const char *title,
				   RooAbsReal& _th1x,  
				   RooAbsReal& _p1, RooAbsReal& _p2, RooAbsReal& _p3) : RooAbsPdf(name, title), 
//TH3* _Hnominal) : RooAbsPdf(name, title), 
  th1x("th1x", "th1x Observable", this, _th1x),
  p1("p1", "p1", this, _p1),
  p2("p2", "p2", this, _p2),
  p3("p3", "p3", this, _p3),
  xBins(0),
  xMax(0),
  xMin(0),
  relTol(1E-12),
  absTol(1E-12)
{
  memset(&xArray, 0, sizeof(xArray));
}
//---------------------------------------------------------------------------
RooInvPowOrd2BinPdf::RooInvPowOrd2BinPdf(const RooInvPowOrd2BinPdf& other, const char* name) :
   RooAbsPdf(other, name), 
   th1x("th1x", this, other.th1x),  
   p1("p1", this, other.p1),
   p2("p2", this, other.p2),
   p3("p3", this, other.p3),
   xBins(other.xBins),
   xMax(other.xMax),
   xMin(other.xMin),
   relTol(other.relTol),
   absTol(other.absTol)
{
  //memset(&xArray, 0, sizeof(xArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] = other.xArray[i];
  }
}
//---------------------------------------------------------------------------
void RooInvPowOrd2BinPdf::setTH1Binning(TH1* _Hnominal){
  xBins = _Hnominal->GetXaxis()->GetNbins();
  xMin = _Hnominal->GetXaxis()->GetBinLowEdge(1);
  xMax = _Hnominal->GetXaxis()->GetBinUpEdge(xBins);
  memset(&xArray, 0, sizeof(xArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] =  _Hnominal->GetXaxis()->GetBinLowEdge(i+1);
  }
}
//---------------------------------------------------------------------------
void RooInvPowOrd2BinPdf::setRelTol(double _relTol){
  relTol = _relTol;
}
//---------------------------------------------------------------------------
void RooInvPowOrd2BinPdf::setAbsTol(double _absTol){
  absTol = _absTol;
}
//---------------------------------------------------------------------------
Double_t RooInvPowOrd2BinPdf::evaluate() const
{
  Double_t integral = 0.0;
  

  Int_t iBin = (Int_t) th1x;
  if(iBin < 0 || iBin >= xBins) {
    //cout << "in bin " << iBin << " which is outside of range" << endl;
    return 0.0;
  }

  
  Double_t xLow = xArray[iBin];
  Double_t xHigh = xArray[iBin+1];
    
  // define the function to be integrated numerically
  InvPowOrd2Function func;
  double params[3];
  params[0] = p1;
  params[1] = p2;
  params[2] = p3;
  func.SetParameters(params);

  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
  ig.SetFunction(func,false);

  
  integral = ig.Integral(xLow,xHigh);
  //Double_t total_integral = ig.Integral(xMin,xMax);

  if (integral>0.0) {
    return integral;
  } else return 0;

}

// //---------------------------------------------------------------------------
Int_t RooInvPowOrd2BinPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if (matchArgs(allVars, analVars, th1x)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooInvPowOrd2BinPdf::analyticalIntegral(Int_t code, const char* rangeName) const{

   Double_t th1xMin = th1x.min(rangeName); Double_t th1xMax = th1x.max(rangeName);
   Int_t iBinMin = (Int_t) th1xMin; Int_t iBinMax = (Int_t) th1xMax;


   Double_t integral = 0.0;
      
   //cout <<  "iBinMin = " << iBinMin << ",iBinMax = " << iBinMax << endl;

   
   // define the function to be integrated numerically  
   InvPowOrd2Function func;
   double params[3];
   params[0] = p1;
   params[1] = p2;
   params[2] = p3;
   func.SetParameters(params);

   
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,absTol,relTol);
  ig.SetFunction(func,false);
    

   if (code==1 && iBinMin<=0 && iBinMax>=xBins){
     integral = ig.Integral(xMin,xMax);
     
   }
   else if(code==1) { 
     for (Int_t iBin=iBinMin; iBin<iBinMax; iBin++){
       
       if(iBin < 0 || iBin >= xBins) {
	 integral += 0.0;
       }
       else{	 
	 Double_t xLow = xArray[iBin];
	 Double_t xHigh = xArray[iBin+1];    
	 integral += ig.Integral(xLow,xHigh);
       }
     }
   } else {
     cout << "WARNING IN RooInvPowOrd2BinPdf: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.0;
   }

   if (integral>0.0) {
     
     return integral;
   } else return 1.0;
}
// //---------------------------------------------------------------------------

