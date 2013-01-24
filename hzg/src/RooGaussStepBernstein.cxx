/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooBernstein.cxx 45779 2012-08-31 15:44:51Z moneta $
 * Authors:                                                                  *
 *   Kyle Cranmer
 *                                                                           *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Bernstein basis polynomials are positive-definite in the range [0,1].
// In this implementation, we extend [0,1] to be the range of the parameter.
// There are n+1 Bernstein basis polynomials of degree n.
// Thus, by providing N coefficients that are positive-definite, there 
// is a natural way to have well bahaved polynomail PDFs.
// For any n, the n+1 basis polynomials 'form a partition of unity', eg.
//  they sum to one for all values of x. See
// http://www.idav.ucdavis.edu/education/CAGDNotes/Bernstein-Polynomials.pdf
// Step function piece means that for x < value f(x) = 0 and for
// x >= value: f(x) = Gauss(x)RooBernstein
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include "TMath.h"
#include "Math/SpecFunc.h"
#include <UWHiggs/hzg/interface/RooGaussStepBernstein.h>
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

using namespace std;

ClassImp(RooGaussStepBernstein)


//_____________________________________________________________________________
RooGaussStepBernstein::RooGaussStepBernstein()
{
}


//_____________________________________________________________________________
RooGaussStepBernstein::RooGaussStepBernstein(const char* name, 
					     const char* title, 
					     RooAbsReal& x, 
					     RooAbsReal& mean,
					     RooAbsReal& sigma,
					     const RooArgList& coefList): 
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _mean("mean","mean of gaussian convolved with polyn",this,mean),
  _sigma("sigma","sigma of gaussian convolved with polyn",this,sigma),
  _coefList("coefficients","List of coefficients",this)
{
  // Constructor
  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RooGaussStepBernstein::ctor(" << GetName() 
	   << ") ERROR: coefficient " << coef->GetName() 
	   << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _coefList.add(*coef) ;
  }
  delete coefIter ;
}



//_____________________________________________________________________________
RooGaussStepBernstein::RooGaussStepBernstein(const RooGaussStepBernstein& other, const char* name) :
  RooAbsPdf(other, name), 
  _x("x", this, other._x), 
  _mean("mean",this,other._mean),
  _sigma("sigma",this,other._sigma),
  _coefList("coefList",this,other._coefList)
{
}


//_____________________________________________________________________________
Double_t RooGaussStepBernstein::evaluate() const 
{
  // something that's positive definite!!!!!!
  Double_t xmin = _x.min(); // old  

  //Double_t xmin = _stepThresh;
  Double_t x = (_x - xmin) / (_x.max() - xmin);   

  Int_t degree = _coefList.getSize() - 1; // n+1 polys of degree n  
  RooFIter iter;
  
  Double_t mean  = ( _mean - xmin )/( _x.max() - xmin ); // scale to [0,1]
  Double_t sigma = _sigma ;// scale to [0,1]
  Double_t gaus_param = 0.5*(mean - x)*(mean - x)/(sigma*sigma);
  /*
  std::cout << mean << ' ' << sigma << ' ' 
	    << x << ' ' << gaus_param << std::endl;
  */
    
  Double_t prefactor,gamma,hyperg_one,hyperg_two,multifact_two,beta,coef;
  
  // iterate through each 
  Double_t result = 0.0;
  for(Int_t i = 0; i <= degree; ++i) {    
    // beta is the sum of bernstein coefficients for a particular polynomial
    // order
    beta = 0.0;
    iter = _coefList.fwdIterator();
    for(Int_t k=i; k <= degree; ++k) {
      // coef is the bernstein polynomial coefficient
      coef = ((RooAbsReal *)iter.next())->getVal();      

      beta += (coef*TMath::Binomial(degree,k)*
	       TMath::Binomial(k,i)*std::pow(-1,k-i));
    }

    prefactor = std::pow(2.0,0.5*(i-1));
    prefactor *= std::exp(-gaus_param);
    prefactor *= std::pow(1.0/(sigma*sigma),-0.5*(i+1));
    
    // gamma function multiplicative piece  
    gamma = ROOT::Math::tgamma(0.5*(1+i));
    
    // the hypergeometric function terms
    hyperg_one = ROOT::Math::conf_hyperg(0.5*(1+i),
					 0.5,
					 gaus_param);
    hyperg_two = ROOT::Math::conf_hyperg(0.5*(2+i),
					 1.5,
					 gaus_param);
    
    multifact_two = std::sqrt(2.0)*std::sqrt(1.0/(sigma*sigma))*(-mean + x);

    result += beta*prefactor*gamma*(hyperg_one + multifact_two*hyperg_two);
    
    
  }

  //std::cout << result << std::endl;

  return result;
}


//_____________________________________________________________________________
Int_t RooGaussStepBernstein::getAnalyticalIntegral(RooArgSet& allVars, 
						   RooArgSet& analVars, 
						   const char* rangeName) const
{
  // No analytical calculation available (yet) of integrals over subranges
  if (rangeName && strlen(rangeName)) {
    return 0 ;
  }
  
  // no analytical calculation of integrals yet, this is evil.
  if (matchArgs(allVars, analVars, _x)) return 0;
  return 0;
}


//_____________________________________________________________________________
Double_t RooGaussStepBernstein::analyticalIntegral(Int_t code, 
						   const char* rangeName) const
{
  assert(code==1) ;
  return 1.0;
}
