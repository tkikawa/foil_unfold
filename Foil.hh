#ifndef Foil_h
#define Foil_h 1

#include "Material.hh"
#include "Cover.hh"

class Foil : public Material
{
public:
  Foil(std::string NAME, std::string XSECFILE, double DENSITY, double a, double ABUNDANCE, double THICKNESS, double ri, double ri_ERR);
  virtual ~Foil();

public:
  double RI;//number of produced RI
  double RI_err;//uncertainty of number of produced RI
  std::vector<Cover> covers;
  void AddCover(Cover cover);
};

#endif
