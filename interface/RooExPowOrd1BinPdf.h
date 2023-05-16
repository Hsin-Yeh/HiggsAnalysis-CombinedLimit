//---------------------------------------------------------------------------
#ifndef HiggsAnalysis_CombinedLimit_RooExPowOrd1BinPdf
#define HiggsAnalysis_CombinedLimit_RooExPowOrd1BinPdf
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH1.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

//---------------------------------------------------------------------------
class RooExPowOrd1BinPdf : public RooAbsPdf
{
public:
   RooExPowOrd1BinPdf() {} ;
   RooExPowOrd1BinPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,	RooAbsReal& _p2);
   RooExPowOrd1BinPdf(const RooExPowOrd1BinPdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooExPowOrd1BinPdf(*this,newname); }
   inline virtual ~RooExPowOrd1BinPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooExPowOrd1BinPdf,1) // RazorModExpBinPdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class ExPowOrd1Function: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = exp(p[0]*x)*pow(x, -1.*p[1]*p[1]);
     return pdf;
   }
   double DoEval(double x) const
   {
     double pdf = exp(pars[0]*x)*pow(x, -1.*pars[1]*pars[1]);
     return pdf;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new ExPowOrd1Function();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 2;
   }
};
