#include <TF1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TColor.h>

#include "Fitter.hh"

Fitter::Fitter()
{
  min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000); // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  setstyle();
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

  min->SetLowerLimitedVariable(0,"A",4.80268e9, 5.e7, 0);
  //min->SetFixedVariable(0,"A",4.80268e9);
  min->SetLowerLimitedVariable(1,"T",80.2798, 0.8, 0);
  //min->SetFixedVariable(1,"T",80.2798);
  min->SetLowerLimitedVariable(2,"B",5.45585e7, 5.e5, 0);
  //min->SetFixedVariable(2,"B",5.45585e7);
  min->SetLowerLimitedVariable(3,"C",1.1423, 0.01, 0);
  //min->SetFixedVariable(3,"C",1.1423);

  // do the minimization
  bool converge = min->Minimize();
  if(!converge){
    std::cout << "Minimization did not converge." << std::endl;
    int status = min->Status();
    if(status == 1) std::cout << "Covariance was made pos defined." << std::endl;
    if(status == 2) std::cout << "Hesse is invalid." << std::endl;
    if(status == 3) std::cout << "Edm is above max." << std::endl;
    if(status == 4) std::cout << "Reached call limit." << std::endl;
    if(status == 5) std::cout << "Other failure." << std::endl;
    exit(1);
  }
}

double Fitter::chi_sq(const double *p){
  double chi2 = 0;
  double RI_exp;
  for(unsigned int i=0;i<foils.size();i++){
    RI_exp = 0;
    for(int j=0;j<foils[i].N-1;j++){
      RI_exp += cover((foils[i].energy[j+1]+foils[i].energy[j])/2,foils[i])*spectrum((foils[i].energy[j+1]+foils[i].energy[j])/2,p)*(foils[i].energy[j+1]-foils[i].energy[j])*(foils[i].xsec[j+1]+foils[i].xsec[j])/2*foils[i].density*foils[i].thickness*foils[i].area*foils[i].abundance/foils[i].A*NA;
    }
    chi2 += pow(foils[i].RI-RI_exp,2)/pow(foils[i].RI*foils[i].RI_err,2);
  }
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
  const double *p   = min->X();
  const double *err = min->Errors();
  const double chi2 = min->MinValue();
  std::cout<<"A = "<<p[0]<<" +/- "<<err[0]<<" ("<<err[0]/p[0]*100<<"% error)"<<std::endl;
  std::cout<<"T = "<<p[1]<<" +/- "<<err[1]<<" ("<<err[1]/p[1]*100<<"% error)"<<std::endl;
  std::cout<<"B = "<<p[2]<<" +/- "<<err[2]<<" ("<<err[2]/p[2]*100<<"% error)"<<std::endl;
  std::cout<<"C = "<<p[3]<<" +/- "<<err[3]<<" ("<<err[3]/p[3]*100<<"% error)"<<std::endl;
  std::cout<<"Minimum chi2 = "<<chi2<<std::endl;
  //min->PrintResults();
}

void Fitter::DrawSpectrum(bool err)
{
  const double *p = min->X();
  TF1 *spec = new TF1("spec", this, &Fitter::spect,1e-4,1e7,npar);
  spec->SetParameters(p[0],p[1],p[2],p[3]);
  spec->SetLineColor(1);
  spec->SetLineWidth(2);

  if(!err){
    TCanvas *c1 = new TCanvas("c1","c1");
    spec->GetXaxis()->SetTitle("Neutron energy (eV)");
    spec->GetYaxis()->SetTitle("Neutron flux (/cm^{2}/eV/sec.)");
    spec->Draw();
    c1->Draw();
    c1->SetLogx();
    c1->SetLogy();
    c1->WaitPrimitive();
  }
  else{
    double xedge[1110], yedge[1001], ymin, ymax, yint, cedge[111];
    for(int k=0;k<111;k++){
      cedge[k] = 1e-4*pow(10,0.1*k);
    }
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
    double ene[1100], cene[110], nom[110], var[110];
    TH2D *hcov = new TH2D("hcov","hcov",110,cedge,110,cedge);
    TH2D *hcor = new TH2D("hcor","hcor",110,cedge,110,cedge);
    
    for(int l=0;l<1100;l++){
      ene[l] = 1e-4*pow(10,0.005+0.01*l);
    }
    for(int l=0;l<110;l++){
      cene[l] = 1e-4*pow(10,0.05+0.1*l);
      nom[l] = spectrum(cene[l],p);
    }

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
	htmp[l]->Fill(spectrum(ene[l],par));
      }
     
      for(int l=0;l<110;l++){
	var[l] = spectrum(cene[l],par);
      }

      for(int l=0;l<110;l++){
	for(int m=0;m<110;m++){
	  hcov->Fill(cene[l],cene[m],(var[l]-nom[l])*(var[m]-nom[m])/nom[l]/nom[m]/100000.);
	}
      }

    }

    TH2D *hcont = new TH2D("hcont","hcont",1100,xedge,1000,yedge);
    TH1D *herr = new TH1D("herr","herr",1100,xedge);
    for(int l=0;l<1100;l++){
      for(int k=0;k<1000;k++){
	hcont->SetBinContent(l+1,k+1,htmp[l]->GetBinContent(k+1)/htmp[l]->GetMaximum());
      }
      herr->SetBinContent(l+1,htmp[l]->GetRMS()/htmp[l]->GetMean()*100);
    }

    double max = 0;
    for(int l=0;l<110;l++){
      for(int m=0;m<110;m++){
	hcor->SetBinContent(l+1,m+1,hcov->GetBinContent(l+1,m+1)/sqrt(hcov->GetBinContent(l+1,l+1)*hcov->GetBinContent(m+1,m+1)));
      if(fabs(hcov->GetBinContent(l+1,l+1))>max)max = fabs(hcov->GetBinContent(l+1,l+1));
      }
    }     

    double param[1]={exp(-0.5)};
    hcont->SetContour(1,param);
    TCanvas *c1 = new TCanvas("c1","c1");
    hcont->GetXaxis()->SetTitle("Neutron energy (eV)");
    hcont->GetYaxis()->SetTitle("Neutron flux (/cm^{2}/eV/sec.)");
    hcont->Draw("cont3");
    spec->Draw("same");
    c1->Draw();
    c1->SetLogx();
    c1->SetLogy();
    c1->WaitPrimitive();

    TCanvas *c2 = new TCanvas("c2","c2");
    herr->GetYaxis()->SetRangeUser(0,100);
    herr->SetLineColor(2);
    herr->SetLineWidth(2);
    herr->GetXaxis()->SetTitle("Neutron energy (eV)");
    herr->GetYaxis()->SetTitle("Neutron flux error (%)");
    herr->Draw("C");
    c2->Draw();
    c2->SetLogx();
    c2->WaitPrimitive();

    TCanvas *c3 = new TCanvas("c3","c3",0,0,600,600);
    hcov->GetZaxis()->SetRangeUser(-max,max);
    hcov->GetXaxis()->SetTitle("Neutron energy (eV)");
    hcov->GetYaxis()->SetTitle("Neutron energy (eV)");
    hcov->Draw("colz");
    c3->Draw();
    c3->SetLogx();
    c3->SetLogy();
    c3->WaitPrimitive();

    TCanvas *c4 = new TCanvas("c4","c4",0,0,600,600);
    hcor->GetZaxis()->SetRangeUser(-1,1);
    hcor->GetXaxis()->SetTitle("Neutron energy (eV)");
    hcor->GetYaxis()->SetTitle("Neutron energy (eV)");
    hcor->Draw("colz");
    c4->Draw();
    c4->SetLogx();
    c4->SetLogy();
    c4->WaitPrimitive();

  }
}

void Fitter::DrawCovariance()
{
  const double *p = min->X();
  TH2D *cov_mat = new TH2D("cov_mat","cov_mat",npar,0,npar,npar,0,npar);
  double max = 0;
  for(int i=0;i<npar;i++){
    for(int j=0;j<npar;j++){
      cov_mat->SetBinContent(i+1,j+1,min->CovMatrix(i,j)/p[i]/p[j]);
      if(fabs(cov_mat->GetBinContent(i+1,j+1))>max)max = fabs(cov_mat->GetBinContent(i+1,j+1));
    }
  }

  std::cout<<"Covariance matrix"<<std::endl;
  std::cout<<cov_mat->GetBinContent(1,1)<<" "<<cov_mat->GetBinContent(2,1)<<" "<<cov_mat->GetBinContent(3,1)<<" "<<cov_mat->GetBinContent(4,1)<<std::endl;
  std::cout<<cov_mat->GetBinContent(1,2)<<" "<<cov_mat->GetBinContent(2,2)<<" "<<cov_mat->GetBinContent(3,2)<<" "<<cov_mat->GetBinContent(4,2)<<std::endl;
  std::cout<<cov_mat->GetBinContent(1,3)<<" "<<cov_mat->GetBinContent(2,3)<<" "<<cov_mat->GetBinContent(3,3)<<" "<<cov_mat->GetBinContent(4,3)<<std::endl;
  std::cout<<cov_mat->GetBinContent(1,4)<<" "<<cov_mat->GetBinContent(2,4)<<" "<<cov_mat->GetBinContent(3,4)<<" "<<cov_mat->GetBinContent(4,4)<<std::endl;

  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
  cov_mat->GetZaxis()->SetRangeUser(-max,max);
  cov_mat->SetMarkerSize(2);
  cov_mat->Draw("text colz");
  c1->Draw();
  c1->WaitPrimitive();
}

void Fitter::DrawCorrelation()
{
  TH2D *cor_mat = new TH2D("cor_mat","cor_mat",npar,0,npar,npar,0,npar);
  for(int i=0;i<npar;i++){
    for(int j=0;j<npar;j++){
      cor_mat->Fill(i,j,min->Correlation(i,j));
    }
  }

  std::cout<<"Correlation matrix"<<std::endl;
  std::cout<<cor_mat->GetBinContent(1,1)<<" "<<cor_mat->GetBinContent(2,1)<<" "<<cor_mat->GetBinContent(3,1)<<" "<<cor_mat->GetBinContent(4,1)<<std::endl;
  std::cout<<cor_mat->GetBinContent(1,2)<<" "<<cor_mat->GetBinContent(2,2)<<" "<<cor_mat->GetBinContent(3,2)<<" "<<cor_mat->GetBinContent(4,2)<<std::endl;
  std::cout<<cor_mat->GetBinContent(1,3)<<" "<<cor_mat->GetBinContent(2,3)<<" "<<cor_mat->GetBinContent(3,3)<<" "<<cor_mat->GetBinContent(4,3)<<std::endl;
  std::cout<<cor_mat->GetBinContent(1,4)<<" "<<cor_mat->GetBinContent(2,4)<<" "<<cor_mat->GetBinContent(3,4)<<" "<<cor_mat->GetBinContent(4,4)<<std::endl;

  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
  cor_mat->GetZaxis()->SetRangeUser(-1,1);
  cor_mat->SetMarkerSize(2);
  cor_mat->Draw("text colz");
  c1->Draw();
  c1->WaitPrimitive();
}

void Fitter::cholcov_conv(double covmat[npar][npar], double cholcovmat[npar][npar])
{
  memset(cholcovmat,0,sizeof(cholcovmat));
  for ( int j=0; j<npar; j++ ) {
    double s = covmat[j][j] ;
    for ( int k=0; k<j; k++ ) {
      s -= cholcovmat[j][k]*cholcovmat[j][k] ;
    }
    if(s<0){
      std::cout << "Covariance matrix is strange. " << j << " " << s << std::endl ;
      exit(1);
    }
    cholcovmat[j][j] = sqrt(s) ;
    
    for ( int i=j+1; i<npar; i++ ) {
      s = covmat[i][j] ;
      for ( int k=0; k<j; k++ ) {
	s -= cholcovmat[i][k]*cholcovmat[j][k] ;
      }
      if ( TMath::Abs(s)<0.000000000001 )
	cholcovmat[i][j] = 0. ;
      else
	cholcovmat[i][j] = s/cholcovmat[j][j] ;
    }
  }
}

void Fitter::setstyle(){
  gStyle->SetLabelFont(132,"axis1");
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameFillColor (0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTextFont(132);
  gStyle->SetTextSize(0.08);
  gStyle->SetLabelFont(132,"x");
  gStyle->SetLabelFont(132,"y");
  gStyle->SetLabelFont(132,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.06,"z");
  gStyle->SetLabelFont(132,"t");
  gStyle->SetTitleFont(132,"x");
  gStyle->SetTitleFont(132,"y");
  gStyle->SetTitleFont(132,"z");
  gStyle->SetTitleFont(132,"t");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleX(0.25);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(132,"pad");
  gStyle->SetMarkerStyle(20);
  gStyle->SetHistLineWidth(1.85);
  gStyle->SetLineStyleString(2,"[12 12]");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPalette(1,0);

  const int NRGBs = 3;
  const int NCont = 255;
  double stops[NRGBs] = { 0.00, 0.50, 1.00 };
  double red[NRGBs]   = { 0.00, 1.00, 1.00 };
  double green[NRGBs] = { 0.00, 1.00, 0.00 };
  double blue[NRGBs]  = { 1.00, 1.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}
