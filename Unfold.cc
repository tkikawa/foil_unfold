/****************************************/
/*                                      */
/*  Unfolding code for the neutron      */
/*  spectrometry with activation foils  */
/*                                      */
/*  Author: Tatsuya Kikawa              */
/*                                      */
/****************************************/

#include <iostream>
#include <iomanip>
#include <vector>
#include <sys/time.h>
#include <TApplication.h>
#include "Fitter.hh"
#include "Material.hh"
#include "Foil.hh"
#include "Cover.hh"
#include "Unfold.hh"

int main(int argc, char *argv[])
{
  if(argc<2) {
    std::cout << "Usage : ./Unfold (input card)" << std::endl;
    exit(1);
  }
  TApplication* app = new TApplication("app", 0, 0);

  struct timeval  start_time, end_time;
  gettimeofday(&start_time, NULL);

  std::string input = argv[1];
  TConfig config;
  ReadInFile(input.c_str(), config);

  int drawspec=2, drawcov=0, drawcor=0;
  std::istringstream(config["GLOBAL"]["Spectrum"])>> drawspec;
  std::istringstream(config["GLOBAL"]["Covariance"])>> drawcov;
  std::istringstream(config["GLOBAL"]["Correlation"])>> drawcor;

  std::string name, xsecfile, cov;
  double density, A, abundance, thickness, area, RI, RI_err;
  std::vector<Cover> covers;
  std::cout<<"Loading cover data"<<std::endl;
  for(std::map<std::string, std::string>::iterator i = config["COVERS"].begin(); i != config["COVERS"].end(); i++){
    name = i->first;
    std::istringstream ss(i->second);
    ss >> xsecfile >> density >> A >> abundance >> thickness;
    Cover *cover = new Cover(name, xsecfile, density, A, abundance, thickness);
    covers.push_back(*cover);
  }

  Fitter *fitter = new Fitter();
  std::cout<<"Loading foil data"<<std::endl;
  for(std::map<std::string, std::string>::iterator i = config["FOILS"].begin(); i != config["FOILS"].end(); i++){
    name = i->first;
    std::istringstream ss(i->second);
    ss >> xsecfile >> density >> A >> abundance >> thickness >> area >> RI >> RI_err;
    Foil *foil = new Foil(name, xsecfile, density, A, abundance, thickness, area, RI, RI_err);
    while(ss >> cov){
      for(unsigned int j=0;j<covers.size();j++){
	if(cov == covers[j].name){
	  foil->AddCover(covers[j]);
	  break;
	}
	else if(j+1 == covers.size()){
	  std::cout<<"Cover "<<cov<<" used in foil "<<name<<" but not defined."<<std::endl;
	  exit(1);
	}
      }
    }
    fitter->AddFoil(*foil);
  }

  std::cout<<"Start unfolding"<<std::endl;
  fitter->Fit();
  if(drawspec>0)fitter->DrawSpectrum(drawspec-1);
  if(drawcov>0)fitter->DrawCovariance();
  if(drawcor>0)fitter->DrawCorrelation();
  fitter->ShowSummary();
  
  gettimeofday(&end_time, NULL);
  double time_diff = end_time.tv_sec - start_time.tv_sec + (double)(end_time.tv_usec - start_time.tv_usec)/1e6;
  std::cout << "CPU time: "<< std::setprecision(4) << time_diff <<" sec."<< std::endl;
  std::cout << "All finished. \\(^o^)/" << std::endl;

  app->Terminate();
  return 0;
}
