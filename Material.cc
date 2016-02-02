#include <fstream>
#include <iostream>
#include "Material.hh"

Material::Material()
{
}

Material::~Material()
{
}

void Material::ReadXsecFile()
{
  double x,y;
  ifstream data;
  data.open(xsecfile.c_str());
  if(!data){
    std::cout<<"Cross section file "<<xsecfile.c_str()<<" is not found."<<std::endl;
    exit(1);
  }
  while(data>>x>>y){
    x *= 1e6; //MeV->eV
    y *= 1e-24; //barn->cm2
    energy.push_back(x);
    xsec.push_back(y);
  }
  N = energy.size();
  data.close();

  if(energy.size()>0&&xsec.size()>0){
    double* xpointer = &(energy.at(0));
    double* ypointer = &(xsec.at(0));
    graph = new TGraph(energy.size(),xpointer,ypointer);
  }
  else{
    std::cout<<xsecfile<<" is empty."<<std::endl;
    exit(1);
  }

  std::cout<<name<<" "<<xsecfile<<std::endl;
}

double Material::Eval(double E){
  double cross;
  if(E<energy[0])cross = xsec[0]*sqrt(energy[0])/sqrt(E);
  else if(E>energy[N-1])cross = xsec[N-1]*sqrt(energy[N-1])/sqrt(E);
  else cross = graph->Eval(E);
  return cross;
}
