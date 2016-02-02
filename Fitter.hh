#ifndef Fitter_h
#define Fitter_h 1

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "Foil.hh"

class Fitter
{
public:
  Fitter();
  virtual ~Fitter();

public:
  static const int npar = 4;//number of fitting parameters
  ROOT::Math::Minimizer* min;
  std::vector<Foil> foils;
  void AddFoil(Foil foil);
  void Fit();
  void ShowSummary();
  void DrawSpectrum(bool err);
  void DrawCovariance();
  void DrawCorrelation();

private:
  static const double NA = 6.02214e23;//Avogadro constant
  static const double kb = 8.61733e-5;//Boltzmann constant
  double chi_sq(const double *p);
  double spectrum(double E, const double *p);
  double spect(double_t *E, double *p);
  double cover(double E, Foil foil);
  void cholcov_conv(double covmat[npar][npar], double cholcovmat[npar][npar]);
};

#endif
