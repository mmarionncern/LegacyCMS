#ifndef _B_H_
#define _B_H_

//ROOT libs
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <Minuit2/FCNBase.h>

//C++ libs
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

class ProbaFit: public ROOT::Minuit2::FCNBase {
   
  public: 
   
  ProbaFit();
  ~ProbaFit();
  
  double operator()(const vector<double>& param) const;   
  double Up() const { return _fUp; }
  // void setPoints(vector<std::pair<vector<float>, vector<int> > >);
  void setPoint(float val, float eval, float p1, float p2);

private: 

  vector<std::pair<vector<float>, vector<int> > > _vals;
  double _fUp;

};


ProbaFit::ProbaFit()
{

  cout << "SimFit::Preparing the FCNBase object to give to Minuit" << endl;

  _fUp = 0.5;
  
} 

ProbaFit::~ProbaFit() {
}

void ProbaFit::setPoint(float val, float eval, float p1, float p2) {

  vector<float> vals(2,0);
  vector<int> bins(2,0);
  
  vals[0]=val;
  vals[1]=eval;
  bins[0]=p1;
  bins[1]=p2;

  std::pair<vector<float>, vector<int> > p(vals,bins);
  _vals.push_back(p);

}

double ProbaFit::operator()( const std::vector<double>& param ) const { //machin à minimiser 

  double chi2=0;

  float val,eval;
  int p1,p2;
  for(unsigned int ip=0;ip<_vals.size();ip++) {
    val=_vals[ip].first[0];
    eval=_vals[ip].first[1];
    p1=_vals[ip].second[0];
    p2=_vals[ip].second[1];
    
    chi2 += pow( val-(param[p1]+param[p2]), 2)/pow(eval,2);
  }

  return chi2;
}


#endif



