//---------------------------------------------------------------------------
#ifndef HiggsAnalysis_CombinedLimit_RooModDijetOrd1ForDiPhotonsBinPdf_h
#define HiggsAnalysis_CombinedLimit_RooModDijetOrd1ForDiPhotonsBinPdf_h
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
class RooModDijetOrd1ForDiPhotonsBinPdf : public RooAbsPdf
{
public:
   RooModDijetOrd1ForDiPhotonsBinPdf() {} ;
   RooModDijetOrd1ForDiPhotonsBinPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,RooAbsReal& _p2, RooAbsReal& _p3, RooAbsReal& _p4);
   RooModDijetOrd1ForDiPhotonsBinPdf(const RooModDijetOrd1ForDiPhotonsBinPdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooModDijetOrd1ForDiPhotonsBinPdf(*this,newname); }
   inline virtual ~RooModDijetOrd1ForDiPhotonsBinPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy p3;        // p2
   RooRealProxy p4;        // p2
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooModDijetOrd1ForDiPhotonsBinPdf,1) // RazorDijetBinPdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class ModDijetOrd1ForDiPhotonsFunction: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = pow(x,p[0]+p[1]*log(x))*pow(1.-x*p[2],-1.*p[3]*p[3]);
     return pdf;
   }
   
   double DoEval(double x) const
   {
     double pdf = pow(x,pars[0]+pars[1]*log(x))*pow(1.-x*pars[2],-1.*pars[3]*pars[3]);
     return pdf;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new ModDijetOrd1ForDiPhotonsFunction();
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
      return 4;
   }
};
