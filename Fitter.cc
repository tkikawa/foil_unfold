#include <TF1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "Fitter.hh"

Fitter::Fitter()
{
  min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000); // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
}

Fitter::~Fitter()
{
}

void Fitter::AddFoil(Foil foil)
{
  foils.push_back(foil);
}

void Fitter::Fit()
{
  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(this, &Fitter::chi_sq,npar);

  min->SetFunction(f);

  min->SetLowerLimitedVariable(0,"A",1., 0.01, 0);
  //min->SetFixedVariable(0,"A",1.11883);
  min->SetLowerLimitedVariable(1,"T",75., 0.01, 0);
  //min->SetFixedVariable(1,"T",8.02051e1);
  min->SetLowerLimitedVariable(2,"B",1.e-2, 0.01, 0);
  //min->SetFixedVariable(2,"B",1.27352e-2);
  min->SetLowerLimitedVariable(3,"C",1., 0.01, 0);
  //min->SetFixedVariable(3,"C",1.14278);

  // do the minimization
  min->Minimize();
}

double Fitter::chi_sq(const double *p){
  double chi2 = 0;
  double RI_exp;
  for(unsigned int i=0;i<foils.size();i++){
    RI_exp = 0;
    for(int j=0;j<foils[i].N-1;j++){
      RI_exp += cover((foils[i].energy[j+1]+foils[i].energy[j])/2,foils[i])*spectrum((foils[i].energy[j+1]+foils[i].energy[j])/2,p)*(foils[i].energy[j+1]-foils[i].energy[j])*(foils[i].xsec[j+1]+foils[i].xsec[j])/2*foils[i].density*foils[i].thickness*foils[i].abundance/foils[i].A*NA;
    }
    chi2 += pow(foils[i].RI-RI_exp,2)/pow(foils[i].RI_err,2);
  }
  chi2 /= foils.size();
  return chi2;
}

double Fitter::spectrum(double E, const double *p){
  double flux;
  flux = p[0]*E/pow((kb*p[1]),2)*exp(-E/(kb*p[1]));
  if(E>kb*p[1]*5)flux += p[2]/pow(E,p[3]);
  return flux;
}

double Fitter::spect(double *E, double *p){
  return spectrum(E[0],p);
}

double Fitter::cover(double E, Foil foil){
  double att = 1;
  for(unsigned int i=0;i<foil.covers.size();i++){
    att *= exp(-foil.covers[i].Eval(E)*foil.covers[i].density*foil.covers[i].thickness*foil.covers[i].abundance/foil.covers[i].A*NA);
  }
  return att;
}


void Fitter::ShowSummary()
{
  const double *p = min->X();
  const double chi2 = min->MinValue();
  //min->PrintResults();
}

void Fitter::DrawSpectrum(bool err)
{
  const double *p = min->X();
  TF1 *spec = new TF1("spec", this, &Fitter::spect,1e-4,1e7,npar);
  spec->SetParameters(p[0],p[1],p[2],p[3]);
  spec->SetLineColor(1);

  if(!err){
    TCanvas *c1 = new TCanvas("c1","c1");
    spec->GetXaxis()->SetTitle("Neutron energy (eV)");
    spec->GetYaxis()->SetTitle("Neutron flux (/eV)");
    spec->Draw();
    c1->Draw();
    c1->SetLogx();
    c1->SetLogy();
    c1->WaitPrimitive();
  }
  else{
    double xedge[1110], yedge[1001], ymin, ymax, yint;
    for(int k=0;k<1101;k++){
      xedge[k] = 1e-4*pow(10,0.01*k);
    }
    ymax = spec->GetMaximum(1e-4,1e7)*100.;
    ymin = spec->GetMinimum(1e-4,1e7)/100.;
    yint = (log(ymax)-log(ymin))/1000.;
    for(int k=0;k<1001;k++){
      yedge[k] = ymin*exp(yint*k);
    }
    TH1D *htmp[1100];
    for(int l=0;l<1100;l++){
      htmp[l] = new TH1D(Form("htmp_%d",l),Form("htmp_%d",l),1000,yedge);
    }
    double cov_mat[npar][npar], chol_mat[npar][npar];
    for(int i=0;i<npar;i++){
      for(int j=0;j<npar;j++){
	cov_mat[i][j] = min->CovMatrix(i,j);
      }
    }
    cholcov_conv(cov_mat, chol_mat);
    double nrand[npar], par[npar];
    TRandom3 rand;
    double E;
    std::cout<<"Start toy MC for error estimation"<<std::endl;
    for(int k=0;k<100000;k++){
      if(k%10000==0)std::cout<<k<<"th toy MC event"<<std::endl;
      for(int i=0;i<npar;i++){
	nrand[i]=rand.Gaus();
      }
      for(int i=0;i<npar;i++){
	par[i]=p[i];
	for(int j=0;j<=i;j++){
	  par[i]+=chol_mat[i][j]*nrand[j];
	}
      }
      if(par[0]<0||par[1]<0||par[2]<0||par[3]<0)continue;
      for(int l=0;l<1100;l++){
	E = 1e-4*pow(10,0.005+0.01*l);
	htmp[l]->Fill(spectrum(E,par));
      }
    }
    TH2D *hcont = new TH2D("hcont","hcont",1100,xedge,1000,yedge);
    for(int l=0;l<1100;l++){
      for(int k=0;k<1000;k++){
	hcont->SetBinContent(l+1,k+1,htmp[l]->GetBinContent(k+1)/htmp[l]->GetMaximum());
      }
    }
    
    double param[1]={exp(-0.5)};
    hcont->SetContour(1,param);
    TCanvas *c1 = new TCanvas("c1","c1");
    hcont->GetXaxis()->SetTitle("Neutron energy (eV)");
    hcont->GetYaxis()->SetTitle("Neutron flux (/eV)");
    hcont->SetFillColor(2);
    hcont->SetFillStyle(3001);
    hcont->Draw("cont3");
    spec->Draw("same");
    c1->Draw();
    c1->SetLogx();
    c1->SetLogy();
    c1->WaitPrimitive();
  }
}

void Fitter::DrawCovariance()
{
  TCanvas *c1 = new TCanvas("c1","c1");
  TH2D *mat = new TH2D("mat","mat",npar,0,npar,npar,0,npar);
  for(int i=0;i<npar;i++){
    for(int j=0;j<npar;j++){
      mat->Fill(i,j,min->CovMatrix(i,j));
    }
  }
  mat->Draw("colz");
  c1->Draw();
  c1->WaitPrimitive();
}

void Fitter::DrawCorrelation()
{
  TCanvas *c1 = new TCanvas("c1","c1");
  TH2D *mat = new TH2D("mat","mat",npar,0,npar,npar,0,npar);
  for(int i=0;i<npar;i++){
    for(int j=0;j<npar;j++){
      mat->Fill(i,j,min->Correlation(i,j));
    }
  }
  mat->Draw("colz");
  c1->Draw();
  c1->WaitPrimitive();
}

void Fitter::cholcov_conv(double covmat[npar][npar], double cholcovmat[npar][npar])
{
  memset(cholcovmat,0,sizeof(cholcovmat));
  for ( Int_t j=0; j<npar; j++ ) {
    Double_t s = covmat[j][j] ;
    for ( Int_t k=0; k<j; k++ ) {
      s -= cholcovmat[j][k]*cholcovmat[j][k] ;
    }
    if(s<0){
      std::cout << "Covariance matrix is strange. " << j << " " << s << std::endl ;
      exit(1);
    }
    cholcovmat[j][j] = sqrt(s) ;
    
    for ( Int_t i=j+1; i<npar; i++ ) {
      s = covmat[i][j] ;
      for ( Int_t k=0; k<j; k++ ) {
	s -= cholcovmat[i][k]*cholcovmat[j][k] ;
      }
      if ( TMath::Abs(s)<0.000000000001 )
	cholcovmat[i][j] = 0. ;
      else
	cholcovmat[i][j] = s/cholcovmat[j][j] ;
    }
  }
}
