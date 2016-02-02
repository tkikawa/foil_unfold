#ifndef Cover_h
#define Cover_h 1

#include "Material.hh"

class Cover : public Material
{
public:
  Cover(std::string XSECNAME, std::string FILE, double DENSITY, double a, double ABUNDANCE, double THICKNESS);
  virtual ~Cover();

public:

};

#endif
