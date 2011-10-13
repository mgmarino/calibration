#include "TApplication.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <vector>
#include "fit2.h"

namespace FIT
{
  std::vector<double> data;
  double Emin,Emax;
  double PARA[4] = {0.0};
  double ERRS[4] = {0.0};
}

void FIT::SetEmin(double min)
{
  FIT::Emin = min;
}

void FIT::SetEmax(double max)
{
  FIT::Emax = max;
}

void FIT::AddDataPoint(double point)
{ 
  if(point < FIT::Emin || point > FIT::Emax){
    std::cout << "Warning: trying to add point which is not within energy limits" << std::endl;
    std::cout << "(Emin|Emax|Epoint)" << std::endl;
    std::cout << FIT::Emin << "|" << FIT::Emax << "|" << point << std::endl << std::endl;
  }
  FIT::data.push_back(point);
}

void FIT::ClearData()
{
  FIT::data.clear();
}

void FIT::SetParams(double *params)
{
  for(int i=0; i<3; i++){
    FIT::PARA[i] = params[i];
  }
}

double* FIT::GetParams()
{
  return FIT::PARA;
}

double* FIT::GetParamErrors()
{
  return FIT::ERRS;
}
  
TH1D* FIT::GetHist(std::string name)
{
  std::cout << "size of data: " << FIT::data.size() << std::endl;
  TH1D* hist = new TH1D(("h"+name).c_str(),("h"+name).c_str(),50,FIT::Emin,FIT::Emax);
  int size = FIT::data.size();
  for(int i=0; i<size; i++){
    hist->Fill(FIT::data[i]);
  }
  int lower = hist->FindBin(FIT::Emin);
  int upper = hist->FindBin(FIT::Emax);
  FIT::PARA[3] = hist->Integral(lower,upper) * (FIT::Emax-FIT::Emin)/50.;
  TF1 *f = GetFunction(FIT::PARA,name);
  hist->GetListOfFunctions()->Add(f);
  return hist;
}

void FIT::fit()
{
  std::cout << "size of data: " << FIT::data.size() << std::endl;
  double arglist[2] = {0.0};
  TMinuit minuit(3);
  minuit.SetFCN(likelihood);
  arglist[0] = 0.5;
  int ierflg = 1;
  minuit.mnexcm("SET ERR",arglist,1,ierflg);
  std::cout << "SET ERR errflag: " << ierflg << std::endl;

  minuit.mnparm(0,"A1",FIT::PARA[0],0.1,0,0,ierflg);
  std::cout << "A1 errflag: " << ierflg << std::endl;
  minuit.mnparm(1,"E",FIT::PARA[1],0.5,0,0,ierflg);
  std::cout << "E errflag: " << ierflg << std::endl;
  minuit.mnparm(2,"s",FIT::PARA[2],0.5,0,0,ierflg);
  std::cout << "s errflag: " << ierflg << std::endl;

  arglist[0] = 500.;
  arglist[1] = 1.0;
  double amin = 0.0,edm = 0.0,errdef = 0.0;
  int nvpar = 0,nparx = 0,icstat = 0;
  minuit.mnexcm("SIMPLEX",arglist,2,ierflg);
  std::cout << "SIMPLEX errflag: " << ierflg << std::endl;
  minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit.mnprin(3,amin);
  minuit.mnexcm("MIGRAD",arglist,2,ierflg);
  std::cout << "MIGRAD errflag: " << ierflg << std::endl;
  minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit.mnprin(3,amin);
  for(int i=0; i<3; i++){
    minuit.GetParameter(i,FIT::PARA[i],FIT::ERRS[i]);
  }
}

void FIT::likelihood(int &npar, double *gin, double &f, double *parameters, int iflag)
{
  std::cout << "trying energy = " << parameters[1] << ", sigma = " << parameters[2] << ", A1 = " << parameters[0] << std::endl;
  double logl = 0.0;
  double integral = FuncIntegral(FIT::Emax,parameters) - FuncIntegral(FIT::Emin,parameters);
  unsigned int ndata = FIT::data.size();
  for(int i=0; i<ndata; i++){
    logl += TMath::Log(function(FIT::data[i],parameters));
  }
	logl -= ndata * TMath::Log(integral);
  f = -logl;
}

double FIT::function(double x, double *parameters)
{
  double value = parameters[0]*TMath::Exp(-(x-parameters[1])*(x-parameters[1])/(2*parameters[2]*parameters[2]));
  value += TMath::Erfc((x-parameters[1])/(TMath::Sqrt2()*parameters[2]));
  //std::cout << "function called, return value = " << value << std::endl;
  return value;
}

double FIT::drawfunction(double *x, double *parameters)
{
  double integral = FuncIntegral(FIT::Emax,parameters) - FuncIntegral(FIT::Emin,parameters);
  return parameters[3]*function(*x,parameters)/integral;
}

double FIT::ErfcIntegral(double x, double a, double b)  //gives the integral over erf(b+a*x)
{
  double value = -b/a*TMath::Erf(b+a*x);
  value += x*TMath::Erfc(b+a*x);
  value += -TMath::Exp(-a*a*x*x-2*a*b*x-b*b) / (a*TMath::Sqrt(TMath::Pi()));
  return value;
}

double FIT::FuncIntegral(double x, double *parameters)
{
  double integral = parameters[0]*parameters[2]*TMath::Sqrt(TMath::Pi()/2.)*TMath::Erf((x-parameters[1])/(TMath::Sqrt2()*parameters[2]));
  integral += ErfcIntegral(x,1/(TMath::Sqrt2()*parameters[2]),-parameters[1]/(TMath::Sqrt2()*parameters[2]));
  return integral;
}

TF1* FIT::GetFunction(double *params, std::string name)
{
  TF1 *myfunc = new TF1(("f"+name).c_str(),drawfunction,FIT::Emin,FIT::Emax,4);
  myfunc->SetParameters(params);
  myfunc->SetParNames("A1","E","s","N");
  myfunc->Draw();
  std::cout << "integral = " << myfunc->Integral(FIT::Emin,FIT::Emax) << std::endl;
  return myfunc;
}
