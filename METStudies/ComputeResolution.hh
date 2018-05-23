#ifndef _LogL_H_
#define _LogL_H_

//ROOT libs
#include <TH1.h>
#include <TF1.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TMath.h>
#include <Minuit2/FCNBase.h>

//C++ libs
#include <fstream>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

bool Dat=false;
bool CorResp=true;

double Chi2(double di, double mci);

double respfuncCor(double pt) {
    float r =1;
  if(Dat) {
    r = 1.00678+ 1.17817/(15.7148 +pt) - 38.5088/(108.963+pt*pt);
 }
  else {
    r = 1.00578+ 1.92346/(17.2487+pt)  -33.5747/(91.5004+pt*pt);
  }
  return r;
}

double fitfunc2(double *x, double *par) {

  double response = 1./respfuncCor(x[0]); 
  if(!CorResp) response=1.;
  int nIT = (par[4])/0.7;

  double f = sqrt( pow( sqrt(x[0])*par[0] + par[1] ,2 )   +par[2]*par[2]*pow(response, 2 ) + (par[3]*par[3]*nIT) *pow(response, 2 ));
  //double f = sqrt( (x[0])*par[0] + par[1]  +par[2]*par[2]*pow(response, 2 ) + (par[3]*par[3]*nIT) *pow(response, 2 ));

  return f;

}



class UseSimultaneousFit: public ROOT::Minuit2::FCNBase {
   
  public: 
   
  UseSimultaneousFit();
  ~UseSimultaneousFit();
  
  double operator()(const vector<double>& param) const;   
  double Up() const { return _fUp; }
  void Parametrize( map<int,TGraphErrors*> graphs, vector<double> fixPar, vector<double> fixParError, bool DSD, bool CResp);
  //double fitfunc(double *x, double *par);


  void LoadLaw();
  void SetFitParameters();

private: 
  vector<pair<int,TGraphErrors*> > _Histos;
  //map<int,TGraphErrors*>::iterator iter;

  vector<double> _FixPar;
  vector<double> _DataFixPar;
  vector<double> _FixParError;
  vector<double> _DataFixParError;

  TF1* _law;
  double _fUp;

};


UseSimultaneousFit::UseSimultaneousFit():
_law(0)
{

  cout << "SimFit::Preparing the FCNBase object to give to Minuit" << endl;

 
  _fUp = 0.5;
  
} 

UseSimultaneousFit::~UseSimultaneousFit() {}

void UseSimultaneousFit::SetFitParameters() {

  cout<<" Parameters filled "<<endl;

  _law->SetParNames("sqrt factor", "sqrt biais","linear factor", "bias");
 
  for(int unsigned i=0;i<_FixPar.size();i++) {
    cout<<_FixPar[i]<<endl;
    _law->FixParameter(i,_FixPar[i]);
  }
  
}


void UseSimultaneousFit::Parametrize(map<int, TGraphErrors*> histos, vector<double> fixPar, vector<double> fixParError, bool DSD, bool CResp) {

  Dat = DSD;
  CorResp = CResp;

  //_Histos = histos;

  //chier...
  map<int,TGraphErrors*>::iterator iT;
  //iter=_Histos.begin();
  // iter++;
  for(iT=histos.begin();//_Histos.begin();
      iT!=histos.end();
      iT++) {
    pair<int,TGraphErrors*> p( iT->first, iT->second );
    _Histos.push_back(p);
  }

  _FixPar = fixPar;
  _FixParError = fixParError;

  cout<<" Parameters loaded "<<endl;

}


void UseSimultaneousFit::LoadLaw() {
  
  _law = new TF1("f2fit",fitfunc2,0.,200.,5); //7
   
}

double UseSimultaneousFit::operator()( const std::vector<double>& param ) const { //machin à minimiser 


  // _law->FixParameter(3,param[0]);
 
  //==================== simul fit
  _law->FixParameter(0,param[0]); //fact sqrt
  _law->FixParameter(1,param[1]); //bias 
  _law->FixParameter(2,param[2]); //noise
  _law->FixParameter(3,param[3]); //euh.. PU
  /* _law->FixParameter(5,param[4]); //OOT PU */
  /* _law->FixParameter(6,param[5]); //fact lin */

  // _law->FixParameter(2,0);
  
  //simul fit =======================
  
  double chi2=0;
  
  for(int i=0;i<_Histos.size();i++) {
    
    

    _law->FixParameter(4, _Histos[i].first );
    
    for(int ib=0;ib<_Histos[i].second->GetN();ib++) {

      double error = _Histos[i].second->GetErrorY(ib);

      double valx, valy;
      _Histos[i].second->GetPoint(ib,valx,valy);
      
      double vallaw = _law->Eval(valx);
   
      if( valx != -1000 ) {
	chi2 += Chi2(valy,vallaw)/pow(error,2);
      }
      
    }
  }
  cout<< chi2<<endl;
  return chi2;
  
}

double Chi2(double di, double mci ){

  double chi2=0;

  chi2 = pow( di-mci , 2);
 
  return chi2;
}



#endif




