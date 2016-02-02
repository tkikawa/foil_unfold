#ifndef Material_h
#define Material_h 1

#include <TGraph.h>
#include <vector>
#include <string>

class Material
{
public:
  Material();
  virtual ~Material();

public:
  std::string name;
  std::string xsecfile;
  double density;//g/cm3
  double A;//mass number
  double abundance;
  double thickness;//cm
  std::vector<double> energy;
  std::vector<double> xsec;
  int N;
  TGraph* graph;
  void ReadXsecFile();
  double Eval(double E);
};

#endif
