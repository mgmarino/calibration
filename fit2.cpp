#include "TApplication.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "fit2.h"

std::vector<double> data;
double Emin,Emax;
double para[4] = {0.0};
double err[4] = {0.0};

void SetEmin(double min)
{
  Emin = min;
}

void SetEmax(double max)
{
  Emax = max;
}

void AddDataPoint(double point)
{ 
  if(point >= Emin && point <= Emax){
    data.push_back(point);
  }
}

void SetParams(double *params)
{
  for(int i=0; i<3; i++){
    para[i] = params[i];
  }
}

double* GetParams()
{
  return para;
}

double* GetParamErrors()
{
  return err;
}
  
TH1D* GetHist()
{
  TH1D* hist = new TH1D("hist","hist",50,Emin,Emax);
  int size = data.size();
  for(int i=0; i<size; i++){
    hist->Fill(data[i]);
  }
  int lower = hist->FindBin(Emin);
  int upper = hist->FindBin(Emax);
  para[3] = hist->Integral(lower,upper) * (Emax-Emin)/50.;
  TF1 *f = GetFunction(para);
  hist->GetListOfFunctions()->Add(f);
  return hist;
}

void fit()
{
  std::cout << "size of data: " << data.size() << std::endl;
  double arglist[2] = {0.0};
  TMinuit minuit(3);
  minuit.SetFCN(likelihood);
  arglist[0] = 0.5;
  int ierflg = 1;
  minuit.mnexcm("SET ERR",arglist,1,ierflg);
  std::cout << "SET ERR errflag: " << ierflg << std::endl;

  minuit.mnparm(0,"A1",para[0],0.1,0,0,ierflg);
  std::cout << "A1 errflag: " << ierflg << std::endl;
  minuit.mnparm(1,"E",para[1],0.5,0,0,ierflg);
  std::cout << "E errflag: " << ierflg << std::endl;
  minuit.mnparm(2,"s",para[2],0.5,0,0,ierflg);
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
  /*
  for(int i=0; i<3; i++){
    minuit.GetParameter(i,para[i],err[i]);
  }
  */
}

void likelihood(int &npar, double *gin, double &f, double *parameters, int iflag)
{
  std::cout << "trying energy = " << parameters[1] << ", sigma = " << parameters[2] << ", A1 = " << parameters[0] << std::endl;
  double logl = 0.0;
  double integral = FuncIntegral(Emax,parameters) - FuncIntegral(Emin,parameters);
  unsigned int ndata = data.size();
  for(int i=0; i<ndata; i++){
    logl += TMath::Log(function(data[i],parameters));
  }
	logl -= ndata * TMath::Log(integral);
  f = -logl;
}

double function(double x, double *parameters)
{
  double value = parameters[0]*TMath::Exp(-(x-parameters[1])*(x-parameters[1])/(2*parameters[2]*parameters[2]));
  value += TMath::Erfc((x-parameters[1])/(TMath::Sqrt2()*parameters[2]));
  //std::cout << "function called, return value = " << value << std::endl;
  return value;
}

double drawfunction(double *x, double *parameters)
{
  double integral = FuncIntegral(Emax,parameters) - FuncIntegral(Emin,parameters);
  return parameters[3]*function(*x,parameters)/integral;
}

double ErfcIntegral(double x, double a, double b)  //gives the integral over erf(b+a*x)
{
  double value = -b/a*TMath::Erf(b+a*x);
  value += x*TMath::Erfc(b+a*x);
  value += -TMath::Exp(-a*a*x*x-2*a*b*x-b*b) / (a*TMath::Sqrt(TMath::Pi()));
  return value;
}

double FuncIntegral(double x, double *parameters)
{
  double integral = parameters[0]*parameters[2]*TMath::Sqrt(TMath::Pi()/2.)*TMath::Erf((x-parameters[1])/(TMath::Sqrt2()*parameters[2]));
  integral += ErfcIntegral(x,1/(TMath::Sqrt2()*parameters[2]),-parameters[1]/(TMath::Sqrt2()*parameters[2]));
  return integral;
}

TF1* GetFunction(double *params)
{
  TF1 *myfunc = new TF1("drawfunction",drawfunction,Emin,Emax,4);
  myfunc->SetParameters(params);
  myfunc->SetParNames("A1","E","s","N");
  myfunc->Draw();
  std::cout << "integral = " << myfunc->Integral(Emin,Emax) << std::endl;
  return myfunc;
}
