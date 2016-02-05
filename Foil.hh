#ifndef Foil_h
#define Foil_h 1

#include "Material.hh"
#include "Cover.hh"

class Foil : public Material
{
public:
  Foil(std::string NAME, std::string XSECFILE, double DENSITY, double a, double ABUNDANCE, double THICKNESS, double AREA, double ri, double ri_ERR);
  virtual ~Foil();

public:
  double RI;//Activation rate [/sec.]
  double RI_err;//Uncertainty of activation rate
  double area;//cm^2
  std::vector<Cover> covers;
  void AddCover(Cover cover);
};

#endif
