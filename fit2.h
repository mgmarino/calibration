#ifndef FIT2_H
#define FIT2_H

#include "TF1.h"
#include "TH1.h"

void SetEmin(double min);
void SetEmax(double max);
void SetParams(double *params);
double* GetParams();
double* GetParamErrors();
void AddDataPoint(double point);
void fit();
void likelihood(int &npar, double *gin, double &f, double *parameters, int iflag);
double function(double x, double *parameters);
double ErfcIntegral(double x, double a, double b);
double FuncIntegral(double x, double *parameters);
TF1* GetFunction(double *params);
TH1D* GetHist();
double drawfunction(double *x, double *parameters);

#endif
