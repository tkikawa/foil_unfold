#include <iostream>
#include "Cover.hh"

Cover::Cover(std::string NAME, std::string XSECFILE, double DENSITY, double a, double ABUNDANCE, double THICKNESS)
{
  name = NAME;
  xsecfile = XSECFILE;
  density = DENSITY;
  A = a;
  abundance = ABUNDANCE;
  thickness = THICKNESS;
  ReadXsecFile();
}
Cover::~Cover()
{
}
