#include <iostream>
#include "Foil.hh"

Foil::Foil(std::string NAME, std::string XSECFILE, double DENSITY, double a, double ABUNDANCE, double THICKNESS, double ri, double ri_ERR)
{
  name = NAME;
  xsecfile = XSECFILE;
  density = DENSITY;
  A = a;
  abundance = ABUNDANCE;
  thickness = THICKNESS;
  RI = ri;
  RI_err = ri_ERR;
  ReadXsecFile();
}
Foil::~Foil()
{
}
void Foil::AddCover(Cover cover)
{
  covers.push_back(cover);
}
