// My libs
#include "ComputeResolution.hh"
//#include "process.C"
// ROOT libs
#include "TTree.h"
#include "TH2.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFitterMinuit.h"
#include "TChain.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMinuit.h"
#include <TVector2.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TPolyLine.h>

// C++ libs
#include <map>
//#include <multimap>
#include <vector>

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TColor.h>

#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <TStyle.h>
#include <RooVoigtian.h>

#include <RooMsgService.h>

#include "RecoilCorrector.cc"

using namespace RooFit;
using namespace std;

vector<int> colors(10,1);
vector<int> mtypes(10,1);

TCanvas* c35;

bool fillH=false;
bool _resp=false;
bool loaded = false;
TTree* Data;

bool apf=false;
bool isLoaded=false;
bool sumUnc=false;

bool noPU=false;

TH1F* puweights;
TH1F* vtxW;
TH1F* qTW;

bool fit;
bool dat;
bool doQtW;
bool Run2012A;
bool Run2012B;
bool ttb;
bool _sumEt;
bool computeUnc=false;


vector<vector<float> > Systematics;

vector<string> metNames;
vector<string> metUncNames;

map<string,string> nameMap;
bool initM=false;

typedef map<string, vector<vector<float> > > recMap;
typedef map<string, vector<vector<float> > >::iterator recMapIter;



//Function declaration
void loadMETTypes(string metMC="", string metData="");

void InitMapName() {

  if(initM) return;

  nameMap[ "mva___mvaIdMedium" ]="mva";
  nameMap[ "mva_RedUnc_atlas_mvaIdMedium" ]="mvaRU";
  nameMap[ "mva_RedAll_atlas_mvaIdMedium" ]="mvaRA";
  nameMap[ "mva_RedAll_vtx_mvaIdMedium" ]="mvaRAn";
  nameMap[ "basic__" ]="jvs";
  nameMap[ "basic_RedAll_atlas" ]="jvsRA";
  nameMap[ "atlasAssoc__" ]="jvf";
  nameMap[ "atlasAssoc_RedAll_atlas" ]="jvfRA";
  nameMap[ "atlasred_mva_atlas_10" ]="mvaRA10";
  nameMap[ "atlasred_mva_atlas_20" ]="mvaRA20";
  nameMap[ "atlasred_mva_atlas_30" ]="mvaRA30";
  nameMap[ "pf_pfType1p0PhiCorrectedMet_" ]="pfType1";
  nameMap[ "pf_pfType1p0PhiCorrectedMet_Smear" ]="pfType1";


  nameMap[ "pf_pfMEtMVA_" ]="mva";
  nameMap[ "pf_pfMEtMVAPhi_" ]="mva";
  nameMap[ "pat_patPFMetMVA_" ]="mva";
  nameMap[ "pat_patPFMetMVAPhi_" ]="mva";
  nameMap[ "pf_noPileUpPFMEt_" ]="noPU";
  nameMap[ "pf_noPileUpPFMEtPhi_" ]="noPU";

  nameMap[ "pat_patPFMetMVASmeared" ]="mva";
  nameMap[ "pat_patPFMetMVAPhiSmeared" ]="mva";
  nameMap[ "pat_patPFMetMVANoSmear" ]="mva";
  nameMap[ "pat_patPFMetMVAPhiNoSmear" ]="mva";
  
  nameMap[ "pat_patPFMetNoPileUpSmeared" ]="mva";
  nameMap[ "pat_patPFMetNoPileUpPhiSmeared" ]="mva";
  nameMap[ "pat_patPFMetNoPileUpNoSmear" ]="mva";
  nameMap[ "pat_patPFMetNoPileUpPhiNoSmear" ]="mva";
  
  nameMap[ "pat_patType1CorrectedPFMet" ]="pfType1";
  nameMap[ "pat_patType1PhiCorrectedPFMet" ]="pfType1";

  nameMap[ "pat_patType1PhiCorrectedPFMetSmeared" ]="pfType1";
  nameMap[ "pat_patType1PhiCorrectedPFMetNoSmear" ]="pfType1";

  nameMap[ "Calo_caloType1CorrectedMet_" ]="caloType1";
  nameMap[ "Calo_caloType1PhiCorrectedMet_" ]="caloType1";
  
  nameMap[ "pf_pfMet_" ]="rawPF";
  nameMap[ "pf_pfMet_" ]="rawPF";
  nameMap[ "pf__" ]="typeIPF";
  nameMap[ "pf__" ]="typeIPF";


  initM=true;

  TGaxis::SetMaxDigits(3);
  //  gStyle->Set

}

float GetVtwW(int nvtx);
void LoadVtxReweight();

void LoadQtReweight();
float GetQtW(float qT);

vector<TPad*> PreparePads(TCanvas *& c2);


float ComputeUncertainties( map<string, TH1*> uncs );

TGraphAsymmErrors* ComputeChristianUncertainties(string type, string comp,TGraph* mc);

TVector2 phiCorrection(TVector2 met, int Nvtx, bool isD, string type);
vector<TGraphAsymmErrors*> getRatios(vector<TGraphAsymmErrors*> graphs);

void  
LoadPUWeights() {
 
  TFile * file = new TFile("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/puWeightsS10.root","READ");
  puweights = (TH1F*)file->Get("pileup");
  
}

float 
SearchWeight(float trueNint) {
  return puweights->GetBinContent( puweights->GetXaxis()->FindBin( trueNint)  );
}



double FWHM(double sigma, double gamma) {

  double f_g = 2*sigma*sqrt(2*log(2));
  double f_l = 2* gamma;

  return 0.5346*2*gamma+sqrt(0.2166*f_l*f_l + f_g*f_g);
}

double FWHMError(/*double sigma, double gamma,*/ double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg) {
  double ef_g =  2*esigma*sqrt(2*log(2));
  double ef_l = 2* egamma;

  double p1 = ef_l*ef_l*Vgg;
  double p2 = ef_g*ef_l*Vsg; //identical (should be)
  double p3 = ef_g*ef_l*Vgs;
  double p4 = ef_g*ef_g*Vss;

  return sqrt(p1 + p2 + p3 + p4);
}


//=== 2D parametrization
bool loadDB=false;
bool corResp=false;

vector<double> dataFixPar(6,0);
vector<double> dataFixParError(6,0);
double FitPar[2][3][4];
double FitParError[2][3][4];

recMap LoadDB(string metname, string comp, string ds, bool FIT);
void GetGeneralPar(string comp, map<int, TGraphErrors*> oMap);
vector<TObject*> FitSigmaP0(string metname, string comp, map<int, TGraphErrors*> );
double fitfunc(double *x, double *par);
double respfunc(double pt);

//For uncertainties
double matrixMC[6][6];
double matrixData[6][6];
double valpData[6];
double valpMC[6];
double errpData[6];
double errpMC[6];
//======================



void ConfigAnalysis(bool isData, bool Fit, bool doqtw, bool unc, bool ttbar ) {

  InitMapName();

  dat=isData;
  //Run2012A = run2012A;
  doQtW = doqtw;
  computeUnc = unc;
  //Run2012B = run2012B;
  ttb = ttbar;
  loadMETTypes();
  fit = Fit;

  int Colors[10]={kMagenta+2, kViolet-6, kBlue+1,
		  kAzure+7, kTeal+5, kOrange,
		  kOrange+7, kRed+1, kPink+5, kPink-7};
  int Mtypes[10]={20,24,21,25,22,26,23,32,34,28};

  int Colors2[4]={kBlue+1,kRed+1,kGreen+1,kViolet-6};

  for(int i=0;i<10;i++)
    {

      colors[i] = Colors[i];
      mtypes[i] = Mtypes[i];
      if(i<4) colors[i] = Colors2[i];
    }

  
}


void loadMETTypes(string metMC, string metData) {

  noPU=false;

  if(isLoaded) return;

  
  if(metMC!="" && metData!="") {
    metNames.push_back(metMC);
    metNames.push_back(metData);

    //FIXME noPU MET very ugly....
    noPU = (metMC.find("NoPileUp")!=(size_t)-1 );

  } 
  else {

  ifstream ifile("metTypes.dat", ios::in);

  if(ifile) {

    while( !ifile.eof() ) {

      string tmp;
      //bool find=false;
      ifile >> tmp;

      if(tmp.substr(0,2) == "//") continue;

      // for(size_t i=0;i<metNames.size();i++) {
      // 	if(metNames[i] == tmp )
      // 	  {find=true; break;}
      // }

      //FIXME, for data
      //      if(!find)
	metNames.push_back(tmp);
    
	//FIXME noPU MET very ugly....
	noPU = (tmp.find("NoPileUp")!=(size_t)-1 );

    }
    
  }
  else { cout<<" No Such file! "<<endl; abort(); }
  
  }

//  "pat_patType1CorrectedPFMet
  metUncNames.push_back("JetEnDown");
  metUncNames.push_back("JetEnUp");
  metUncNames.push_back("JetResDown");
  //metUncNames.push_back("JetResUp");
  metUncNames.push_back("UnclusteredEnDown");
  metUncNames.push_back("UnclusteredEnUp");
  metUncNames.push_back("MuonEnDown");
  metUncNames.push_back("MuonEnUp");
  
  isLoaded=true;

  cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<noPU<<endl;

}

float XSectWeight(int i, vector<int> mcevts, int nZ, int nTT) {

  if( i<= mcevts[0] ) return 1./(nZ/3503.71)*10259;
  if(i> mcevts[0] && i<= mcevts[1]) return 1./(nTT/234.)*10259;
  
  return 1.;
}

bool isMCZ(int i, vector<int> mcevts) {
  if( i<= mcevts[0] ) return 1;
  else return 0;
}

void quadSumUnc(TGraphAsymmErrors*& gr, TGraphAsymmErrors* uncgr) {

  cout<<" quadratic sum "<<endl;

  double x,y, esyh, esyl, euyh, euyl;
  for(int i=0;i<gr->GetN();i++) {

    cout<<i<<endl;

    esyh = gr->GetErrorYhigh(i);
    esyl = gr->GetErrorYlow(i);

    cout<<" glou "<<endl;

    euyh = uncgr->GetErrorYhigh(i);
    euyl = uncgr->GetErrorYlow(i);

    cout<<" gliu "<<endl;

    gr->SetPointEYhigh(i, sqrt( esyh*esyh + euyh*euyh ) );
    gr->SetPointEYlow(i, sqrt( esyl*esyl + euyl*euyl ) );

  }

}


TTree* LoadTree(bool b=false) {

  //string rep = "vPt20/";
  string rep = "/afs/cern.ch/user/m/mmarionn/workspace/private/ANSNtuples53Xv14/";

  if(!(Run2012A && Run2012B) ) {
    if(Run2012A) Run2012B=false;
    if(Run2012B) Run2012A=false;
  }
  LoadQtReweight();
  vector<int> MCNevt;
  TChain* ITreeD=new TChain("metTree");
 
  double ZNent=0; double TTBNent=0;
  TH1F* ZEntries=NULL;
  TH1F* ttbarEntries=NULL;
  if(dat) {
    string pathdata=rep;//"~/workspace/private/ANSNtuples/";
    // ITreeD->Add( (pathdata+"DoubleMu_190450_193686_1.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_B_193752_197044_1.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_B_193752_197044_2.root").c_str() );
   
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_190456_193621_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_A_r_190782_190949_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_B_190456_196531_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_B_190456_196531_2.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_C_v1_197770_198913_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_C_v2_198934_201229_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_C_v2_198934_201229_2.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_C_v2_201230_201678_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv4/DoubleMu_C_v2_201679_202016_1.root").c_str() );

    // ITreeD->Add( (pathdata+"../53Xv8/DoubleMu_190456_193621_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv8/DoubleMu_A_r_190782_190949_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv8/DoubleMu_B_190456_196531_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv8/DoubleMu_B_190456_196531_2.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv8/DoubleMu_C_v1_198022_198523_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv8/DoubleMu_C_v2_190456_204567_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv8/DoubleMu_C_v2_190456_204567_2.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv8/DoubleMu_C_v2_190456_204567_3.root").c_str() );
    
    // ITreeD->Add( (pathdata+"DoubleMu_190456_193621_1.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_A_r_190782_190949_1.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_B_190456_196531_1.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_B_190456_196531_2.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_C_v1_198022_198523_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv12/DoubleMu_C_v2_190456_206539_1.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv12/DoubleMu_C_v2_190456_206539_2.root").c_str() );
    // ITreeD->Add( (pathdata+"../53Xv12/DoubleMu_C_v2_190456_206539_3.root").c_str() );

    
    // ITreeD->Add( (pathdata+"DoubleMu_C_v2_190456_206940_1.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_C_v2_190456_206940_2.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_C_v2_190456_206940_3.root").c_str() );
    
    ITreeD->Add( (pathdata+"DoubleMu_190456_193621_1.root").c_str() );
    ITreeD->Add( (pathdata+"DoubleMu_A_r_190782_190949_1.root").c_str() );
    ITreeD->Add( (pathdata+"DoubleMu_B_190456_196531_1.root").c_str() );
    //ITreeD->Add( (pathdata+"DoubleMu_B_190456_196531.root").c_str() );
    ITreeD->Add( (pathdata+"DoubleMu_C_v1_198022_198523_1.root").c_str() );
    ITreeD->Add( (pathdata+"DoubleMu_C_v2_190456_208686_1.root").c_str() );
  

    // ITreeD->Add( (pathdata+"DoubleMu_C_v2_190456_208686_1.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_C_v2_190456_208686_2.root").c_str() );
    // ITreeD->Add( (pathdata+"DoubleMu_C_v2_190456_208686_3.root").c_str() );

    cout<<" ouééééééééééééééééééé data "<<endl;
    
  }
  else {
    // if(Run201A) rep = "_A";
//     if(Run2011B) rep = "_B";
    //rep =  "v1_2011"; //FIXME

    //mc_v6_2011_testJEC
    //Z_2m_V13_v2
    rep = "/afs/cern.ch/user/m/mmarionn/workspace/private/ANSNtuples53Xv16/"; //53Xv6
    LoadPUWeights();
    LoadVtxReweight();
   
    ITreeD->Add( (rep+"ZJets_2l.root").c_str() );
    //ITreeD->Add( (rep+"ZJets_2l.root").c_str() );
    MCNevt.push_back( ITreeD->GetEntries() );
    //TTbar
     if(ttb) {
       ITreeD->Add( (rep+"TTbarJets.root").c_str() );
       MCNevt.push_back( ITreeD->GetEntries() );
     }
     //ITreeD->Add( ("Ntuple.root") );
    

    // cout<<" MC event : "<<MCNevt[0];
    if(MCNevt.size()==2) cout<<" / "<<MCNevt[1]<<endl;
    else cout<<endl;

    TFile* ZFile=new TFile((rep+"ZJets_2l.root").c_str(),"READ");
    // TFile* ZFile=new TFile(("Ntuple.root"),"READ"); cout<<ZFile<<endl;
    ZEntries = (TH1F*)ZFile->Get("all/data/nEvtProc__run"); cout<<ZEntries<<endl;
    ZNent = ZEntries->Integral(0,10001);
    cout<<ZNent<<endl;
    ZFile->Close();
     if(ttb) {
      
       TFile* TTBFile=new TFile((rep+"TTbarJets.root").c_str(),"READ");
       ttbarEntries= (TH1F*)TTBFile->Get("all/data/nEvtProc__run");
       TTBNent = ttbarEntries->Integral(0,10001);
       TTBFile->Close();
       //cout<<" Adding ttbar : Z-> "<<ZNent<<" -> "<<TTBNent<<endl;
     }
    cout<<" ouééééééééééééééééééé MC "<<endl;
  }

  if(b)
    return ITreeD;

  cout<<" Trees Merged   "<<endl;
  float ZPt;
  float ZMass;
  float ZPhi;

  int njet;

  int run;
  int event;

  //general variables

  // ITreeD->SetBranchStatus("add*",0);
  // ITreeD->SetBranchStatus("miss*",0);
  //cout<<" coin "<<endl;
  ITreeD->SetBranchAddress("phiZ",&ZPhi);
  ITreeD->SetBranchAddress("massZ",&ZMass);
  ITreeD->SetBranchAddress("ptZ",&ZPt);

  ITreeD->SetBranchAddress("run",&run);
  ITreeD->SetBranchAddress("event",&event);
  ITreeD->SetBranchAddress("njet",&njet);
  //ITreeD->SetBranchAddress("LS",&LS);

  //muon needed for additional smearing
  float pt1,eta1,phi1;
  float pt2,eta2,phi2;
  ITreeD->SetBranchAddress("pt_m1",&pt1);
  ITreeD->SetBranchAddress("eta_m1",&eta1);
  ITreeD->SetBranchAddress("phi_m1",&phi1);
  ITreeD->SetBranchAddress("pt_m2",&pt2);
  ITreeD->SetBranchAddress("eta_m2",&eta2);
  ITreeD->SetBranchAddress("phi_m2",&phi2);
  
  //MissingETs

  map<string, float > metPts; 
  map<string, float > metPhis; 
  map<string, float > metSumEts;
  map<string, float > metXs;
  map<string, float > metYs;

  map<string, float > metPtsUnc; 
  map<string, float > metPhisUnc; 

  map<string, float >::iterator iter;
  map<string, float >::iterator iter2;
  map<string, float >::iterator itUM;
  //map<string, float >::iterator iter3;

  //cout<<" coin "<<endl;
  //Factors & others
  // float fAtlas;
  // float fVtx;
  int nVtx;
  float trueNint;

  // ITreeD->SetBranchAddress("fAtlas",&fAtlas);
  // ITreeD->SetBranchAddress("fVtx",&fVtx);
  ITreeD->SetBranchAddress("nVtx",&nVtx);
  ITreeD->SetBranchAddress("trueNint",&trueNint);

  //Weights
  //cout<<" coin "<<endl;
  //internal objects
  TVector2 Z(0,0);
  map<string, TVector2> mets;
  map<string, TVector2> recoils;
  map<string, float> uParas;
  map<string, float> uPerps;
  map<string, float> redUParas;
  map<string, float>::iterator iterF;
  map<string, float>::iterator iterFU;
  
  map<string, float> uParasUnc;
  map<string, float> uPerpsUnc;
  map<string, float> redUParasUnc;
  map<string, float> sumEt;

  //bool phi=false;

  map<string,string> rNames;
  map<string,string> rNamesUnc;

 
  map<string, TH2F*> histoPara;
  map<string, TH2F*> histoPerp;


  //Recoil corrector
  // RecoilCorrector rcCor("/home/mmarionn/Documents/CMS/LAURE/modules/RecoilCorrector_v6/recoilfits/recoilfit_zmm53X_20pv_njet.root");

  // // Works perfectly with the mvaMET (no jet smearing)
  // rcCor.addMCFile("/home/mmarionn/Documents/CMS/LAURE/modules/RecoilCorrector_v6/recoilfits/recoilfit_zmm53X_20pv_njet.root");
  // rcCor.addDataFile("/home/mmarionn/Documents/CMS/LAURE/modules/RecoilCorrector_v6/recoilfits/recoilfit_datamm53X_20pv_njet.root");

  RecoilCorrector rcCorMva("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/RecoilCorrector_v8/recoilfit_mc_53X_cv_v4_njet.root");//recoilfit_data_53X_cv_np_inc_v2_njet.root
  rcCorMva.addMCFile("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/RecoilCorrector_v8/recoilfit_mc_53X_cv_v4_njet.root");
  rcCorMva.addDataFile("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/RecoilCorrector_v8/recoilfit_data_53X_cv_v4_njet.root");

  RecoilCorrector rcCor("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/RecoilCorrector_v8/recoilfit_mc_53X_cv_np_v3_njet.root");//recoilfit_data_53X_cv_np_inc_v2_njet.root
  rcCor.addMCFile("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/RecoilCorrector_v8/recoilfit_mc_53X_cv_np_v3_njet.root");
  rcCor.addDataFile("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/RecoilCorrector_v8/recoilfit_data_53X_cv_np_v3_njet.root");

  string postfix="Smeared";
  if(dat) postfix = "NoSmear" ;

  if(!apf) postfix="";

  ////cout<<" coin "<<endl;
  for(size_t i=0;i<metNames.size();i++) {
   
    if(metNames[i].find("Phi")!=(size_t)-1) {
      //  phi =true;
      uParas[ metNames[i] ] = 0;
      uPerps[ metNames[i] ] = 0;
      redUParas[ metNames[i] ] = 0;
      sumEt[ metNames[i] ] = 0;
      string s = metNames[i];
      rNames[ metNames[i] ] = s.replace( metNames[i].find("Phi"), 3, "") ;

      if(apf) rNames[ metNames[i] ]+=postfix;

      cout<<metNames[i]<<"   "<<rNames[ metNames[i] ]<<endl;
      
      if(fillH) {
	//cout<<metNames[i]<<endl;
	histoPara[ metNames[i] ] = new TH2F( (metNames[i]+"para").c_str(), (metNames[i]+"para").c_str(), 51,-0.5,50.5, 200,-100.,100.);
	histoPerp[ metNames[i] ] = new TH2F( (metNames[i]+"perp").c_str(), (metNames[i]+"perp").c_str(), 51,-0.5,50.5, 200,-100.,100.);
      }

    }
    else {
      uParas[ metNames[i] ] = 0;
      uPerps[ metNames[i] ] = 0;
      redUParas[ metNames[i] ] = 0;
      sumEt[ metNames[i] ] = 0;
      rNames[ metNames[i] ] = metNames[i] ;

      if(apf) rNames[ metNames[i] ]+=postfix;

      if(fillH) {
	//cout<<metNames[i]<<endl;
	histoPara[ metNames[i] ] = new TH2F( (metNames[i]+"para").c_str(), (metNames[i]+"para").c_str(), 51,-0.5,50.5, 200,-100.,100.);
	histoPerp[ metNames[i] ] = new TH2F( (metNames[i]+"perp").c_str(), (metNames[i]+"perp").c_str(), 51,-0.5,50.5, 200,-100.,100.);
      }
    }
     //cout<<" ddcoin "<<endl;
    //uncertainties
    if(computeUnc && !dat && !noPU) {
      for(size_t j=0;j<metUncNames.size();j++) {
	//cout<<metUncNames[j]+"_"+metNames[i]<<endl;
	if(metNames[i].find("Phi")!=(size_t)-1) {
	  //phi =true;
	  uParasUnc[ metUncNames[j]+"_"+metNames[i] ] = 0;
	  uPerpsUnc[ metUncNames[j]+"_"+metNames[i] ] = 0;
	  redUParasUnc[ metUncNames[j]+"_"+metNames[i] ] = 0;
	  if(i==0) {
	    string s = metUncNames[j];
	    rNamesUnc[ metUncNames[j] ] = s;//s.replace( (metUncNames[j]+"_"+metNames[i]).find("Phi"), 3, "") ;
	    //if(apf) rNamesUnc[ metUncNames[i] ] = postfix + rNamesUnc[ metUncNames[i] ];
	  }
	}
	else {
	  uParasUnc[ metUncNames[j]+"_"+metNames[i] ] = 0;
	  uPerpsUnc[ metUncNames[j]+"_"+metNames[i] ] = 0;
	  redUParasUnc[ metUncNames[j]+"_"+metNames[i] ] = 0;
	  if(i==0) {
	    rNamesUnc[ metUncNames[j] ] = metUncNames[j] ;
	    //if(apf) rNamesUnc[ metUncNames[i] ] = postfix + rNamesUnc[ metUncNames[i] ];
	  }
	}
      }
    }

  }
  
  //cout<<" ddcoddddin "<<endl;
  //initialisation of MET branches
  for(size_t i=0;i<metNames.size();i++) {
    
    metPts[ metNames[i] ] = 0;
    metPhis[ metNames[i] ] = 0;
    metSumEts[ metNames[i] ] = 0;
    ITreeD->SetBranchAddress( (rNames[metNames[i]]+"__pt").c_str(), &metPts[ metNames[i] ] );
    ITreeD->SetBranchAddress( (rNames[metNames[i]]+"__phi").c_str(), &metPhis[ metNames[i] ] );
    ITreeD->SetBranchAddress( (rNames[metNames[i]]+"__sumEt").c_str(), &metSumEts[ metNames[i] ] );

    cout<<metNames[i]<<" --> "<<rNames[metNames[i]]<<endl;
  }
  // metPts[ "pat_patType1CorrectedPFMet_" ] = 0;
  // metPhis[ "pat_patType1CorrectedPFMet_" ] = 0;
  // metSumEts[ "pat_patType1CorrectedPFMet_" ] = 0;
  // ITreeD->SetBranchAddress( "pat_patType1CorrectedPFMet__pt", &metPts[ "pat_patType1CorrectedPFMet_" ] );
  // ITreeD->SetBranchAddress( "pat_patType1CorrectedPFMet__phi", &metPhis[ "pat_patType1CorrectedPFMet_" ] );
  // ITreeD->SetBranchAddress( "pat_patType1CorrectedPFMet__sumEt", &metSumEts[ "pat_patType1CorrectedPFMet_" ] );

  if(computeUnc && !dat && !noPU) {
    cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!! Uncertainties "<<endl;
    for(size_t i=0;i<metUncNames.size();i++) {
      //cout<<rNamesUnc[metUncNames[i]]<<"   "<<metUncNames[i]<<endl;
      string base=rNames[metNames[0]].substr(0,rNames[metNames[0]].size()-postfix.size());
      cout<<base+rNamesUnc[metUncNames[i]]<<"   "<<postfix<<endl;
      metPtsUnc[ metUncNames[i] ] = 0;
      metPhisUnc[ metUncNames[i] ] = 0;
      ITreeD->SetBranchAddress( (base+rNamesUnc[metUncNames[i]]+postfix+"__pt").c_str(), &metPtsUnc[ metUncNames[i] ] );
      ITreeD->SetBranchAddress( (base+rNamesUnc[metUncNames[i]]+postfix+"__phi").c_str(), &metPhisUnc[ metUncNames[i] ] );
    }
  }

 //cout<<" ddcccccccccoddddin "<<endl;

  //output Tree
  TTree* ZData=new TTree("ZData","data");
  
  float oZpt;
  int onVtx;
  float weight;

  ZData->Branch("Zpt",&oZpt,"oZpt/F");
  ZData->Branch("nVtx",&onVtx,"onVtx/I");
  ZData->Branch("weight",&weight,"weight/F");
  
  //cout<<" flibli "<<endl;

  for(iterF=uParas.begin();iterF!=uParas.end();iterF++) {
    //cout<<iterF->first<<endl;
    ZData->Branch( (iterF->first + "_upara").c_str(), &uParas[  iterF->first ]  );
    ZData->Branch( (iterF->first + "_uperp").c_str(), &uPerps[  iterF->first ] );
    ZData->Branch( (iterF->first + "_redupara").c_str(), &redUParas[  iterF->first ] );
    ZData->Branch( (iterF->first + "_sumEt").c_str(), &sumEt[  iterF->first ] );
    
    //uncertainties
    for(iterFU=uParasUnc.begin();iterFU!=uParasUnc.end();iterFU++) {
      //cout<<" ----> "<<iterFU->first<<endl;
      ZData->Branch( (iterFU->first + "_upara").c_str(), &uParasUnc[  iterFU->first ]  );
      ZData->Branch( (iterFU->first + "_uperp").c_str(), &uPerpsUnc[  iterFU->first ] );
      ZData->Branch( (iterFU->first + "_redupara").c_str(), &redUParasUnc[  iterFU->first ] );
      //cout<<" booking uncertainties "<<(iterFU->first + "_uperp")<<endl;
    }

  }

  //cout<<" flibddddli "<<endl;

  TLorentzVector l1,l2;
  float smearB = 0.011; //0.011
  float smearE = 0.011; //0.015
  TRandom3 rand,rand2,rnd;
  rand.SetSeed(18184);
  rand2.SetSeed(18630);
  
  TVector2 tmpunc;

  float metUnc,angleUnc;
  
  cout<<" ============> "<<_sumEt<<endl;

  for(int id=0;id< ITreeD->GetEntries();id++) {

    ITreeD->GetEntry(id);
 
    //Z.SetMagPhi(ZPt,ZPhi);
    //cout<<" event "<<id<<endl;
    //if(id>20000) continue;
   
    onVtx = nVtx;

    if(ZMass < 60 || ZMass > 120) continue;
       //if(njet!=0) continue;
 
    l1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
    l2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);
    if(dat) {
      l1 *= (fabs(eta1)>1.3)?1.000824:1.001305;
      l2 *= (fabs(eta2)>1.3)?1.000824:1.001305;
      weight =1.*(doQtW?GetQtW( (l1+l2).Pt() ):1);//*GetQtW( (l1+l2).Pt() ); //FIXME
    }
    else {
      //cout<<id<<"  ";
   
      l1 *= rnd.Gaus(1,(fabs(eta1)<1.3)?smearB:smearE); 
      l2 *= rnd.Gaus(1,(fabs(eta2)<1.3)?smearE:smearB); 
      weight = SearchWeight(trueNint)*XSectWeight(id, MCNevt, ZNent , TTBNent )*GetVtwW(nVtx)*(doQtW?GetQtW( (l1+l2).Pt() ):1); //FIXME
      //  cout<<weight<<"   "<<l1.Pt()<<"   "<<l2.Pt()<<endl;
    }

    //Fixme
    if(doQtW && (l1+l2).Pt()<100) continue;

    //FIXME pt lept > 30 GeV
    if(l1.Pt()<20 || l2.Pt()<20) continue;

    Z = (l1+l2).Vect().XYvector();
    oZpt = Z.Mod();

    for(iter=metPts.begin();iter!=metPts.end();iter++) {
     
      //if(iter->first=="pat_patType1CorrectedPFMet_") continue;

      TVector2 tmp(0,0); tmp.SetMagPhi( iter->second , metPhis[ iter->first ] );

      // if(id<100)
      // 	cout<<iter->first<< "== "<<tmp.X()<<"  "<<tmp.Y()<<endl;
      if(iter->first.find("Phi")!=(size_t)-1)
	tmp = phiCorrection(tmp, nVtx, dat, iter->first);
      // if(id<100)
      // 	cout<<iter->first<< ">> "<<tmp.X()<<"  "<<tmp.Y()<<endl;
      mets[ iter->first ] = tmp;
      TVector2 tmp2(0,0); tmp2 -= (Z + tmp);

      if(corResp)
	tmp2 /=respfunc(ZPt);

      recoils[ iter->first ] = tmp2;
	
      uParas[ iter->first ] = _sumEt?tmp.X():(( Z*tmp2)/Z.Mod());
      tmp2 = tmp2.Rotate( TMath::Pi()/2);
      uPerps[ iter->first ] = _sumEt?tmp.Y():(( Z*tmp2)/Z.Mod());

      double tmet = iter->second;
      double tmetphi = metPhis[ iter->first ];
      double tupara = uParas[ iter->first ];
      double tuperp = uPerps[ iter->first ];

      //recoil correction
      if(!dat && isMCZ(id, MCNevt) && metNames[0].find("NoPileUp")!=(size_t)-1 ) {
	if(id<5) cout<<" correcting recoil -!-> "<<tupara<<endl;
	rcCor.CorrectType1( tmet, tmetphi,Z.Mod(),Z.Phi(),
			    Z.Mod(),Z.Phi(),tupara,tuperp,0.,0.,njet);
	if(id<5) cout<<" correcting recoil -!-> "<<tupara<<endl;
      }
      if(!dat && isMCZ(id, MCNevt) && metNames[0].find("MVA")!=(size_t)-1 ) {
      	if(id<5) cout<<" correcting recoil -!-> "<<tupara<<endl;
      	rcCorMva.CorrectType1( tmet, tmetphi,Z.Mod(),Z.Phi(),
      			    Z.Mod(),Z.Phi(),tupara,tuperp,0.,0.,njet);
	if(id<5) cout<<" correcting recoil -!-> "<<tupara<<endl;
      }

      tmp2.SetMagPhi( tmet, tmetphi );

      uParas[ iter->first ] = _sumEt?tmp.X():tupara;
      uPerps[ iter->first ] = _sumEt?tmp.Y():tuperp;

      redUParas[ iter->first ] = _sumEt?tmp.X():(uParas[ iter->first ]+oZpt);

      sumEt[ iter->first ] = metSumEts[ iter->first ]-pt1-pt2;
	

	// if(id<100)
	//   cout<<_sumEt<<"   "<<tmp.X()<<"   "<<uParas[ iter->first ]<<" ---- "<<tmp.Y()<<"   "<<uPerps[ iter->first ]<<" // "<<weight<<endl;

	// if(Z.Mod() >300 && Z.Mod() <400 ) {
	//   cout<<Z.Mod()<<"   "<<-uParas[ iter->first ]<<"   "<<-uParas[ iter->first ]/Z.Mod()<<endl;
	// }

	//Fill the 2D RMS/Vtx histos 
	if(fillH) {
	  if(Z.Mod() <300 ) {
	    if( uParas[ iter->first ]<75 && uParas[ iter->first]>-500 )
	      histoPara[ iter->first ]->Fill( onVtx, redUParas[ iter->first ], weight );
	    if( uPerps[ iter->first ]<75 && uPerps[ iter->first]>-75 )
	      histoPerp[ iter->first ]->Fill( onVtx, uPerps[ iter->first ], weight );
	  }
	}
	
	if(computeUnc && !dat && !noPU) {
	  
	  // tmp.SetMagPhi( metPts[ "pat_patType1CorrectedPFMet_" ] , metPhis[ "pat_patType1CorrectedPFMet_" ] );
	  // if(iter->first.find("Phi")!=(size_t)-1)
	  //   tmp = phiCorrection(tmp, nVtx, dat, iter->first);
	  
	  for(itUM=metPtsUnc.begin();itUM!=metPtsUnc.end();itUM++) {
	    tmpunc.SetMagPhi( itUM->second, metPhisUnc[ itUM->first ]);
	    if(iter->first.find("Phi")!=(size_t)-1)
	      tmpunc = phiCorrection(tmpunc, nVtx, dat, iter->first);  
	    // if(id<100)
	    //   cout<<itUM->first<<" --> "<<tmpunc.Mod()<<"  "<<tmpunc.Phi()<<endl;

	    metUnc = tmpunc.Mod();
	    angleUnc = tmpunc.Phi(); //works only on 0 to 2 pi
	    
	    //	    cout<<iter->first<<"   "<<itUM->first<<"   "<<metUnc<<"  // "<<mets[ iter->first ].Mod()<<"   "<<tmpunc.Mod()<<"   "<<tmp.Mod()<<"  "<<itUM->second<<"   "<<endl;

	    tmpunc.SetMagPhi(metUnc,angleUnc);
	   
	    tmp2.SetMagPhi(0,0);
	    tmp2 -= (Z + tmpunc);

	    if(corResp)
	      tmp2 /=respfunc(ZPt);
	    
	    uParasUnc[ itUM->first+"_"+iter->first ] = ( Z*tmp2)/Z.Mod();
	    tmp2 = tmp2.Rotate( TMath::Pi()/2);
	    uPerpsUnc[ itUM->first+"_"+iter->first ] = ( Z*tmp2)/Z.Mod();
	    redUParasUnc[ itUM->first+"_"+iter->first ] = uParasUnc[ itUM->first+"_"+iter->first ]+ZPt;	
	    
	    double tmetu = tmpunc.Mod();
	    double tmetphiu = tmpunc.Phi();
	    double tuparau = uParasUnc[ itUM->first+"_"+iter->first ];
	    double tuperpu = uPerpsUnc[ itUM->first+"_"+iter->first ];
	    
	    //recoil correction
	    if(!dat && isMCZ(id, MCNevt) && metNames[0].find("NoPileUp")!=(size_t)-1 ) {
	      if(id<5) cout<<" correcting recoil uncertainty "<<itUM->first+"_"+iter->first<<"  "<<tuparau<<endl;
	      rcCor.CorrectType1( tmetu, tmetphiu,Z.Mod(),Z.Phi(),
				  Z.Mod(),Z.Phi(),tuparau,tuperpu,0.,0.,njet);
	      if(id<5) cout<<" --> "<<tuparau<<endl;
	    }
	     if(!dat && isMCZ(id, MCNevt) && metNames[0].find("MVA")!=(size_t)-1 ) {
	       if(id<5) cout<<" correcting recoil uncertainty "<<tuparau<<endl;
	       rcCorMva.CorrectType1( tmetu, tmetphiu,Z.Mod(),Z.Phi(),
	     			     Z.Mod(),Z.Phi(),tuparau,tuperpu,0.,0.,njet);
	       if(id<5) cout<<" --> "<<tuparau<<endl;
	     }
	    
	    tmp2.SetMagPhi( tmetu, tmetphiu ); 

	    uParasUnc[ itUM->first+"_"+iter->first ] = _sumEt?tmp.X():tuparau;
	    uPerpsUnc[ itUM->first+"_"+iter->first ] = _sumEt?tmp.Y():tuperpu;

	    redUParasUnc[ itUM->first+"_"+iter->first ] = _sumEt?tmp.X():(uParasUnc[ itUM->first+"_"+iter->first ]+oZpt);
	  
	    // if(id<100)
	    //   cout<<itUM->first<<" --> "<<tmpunc.Mod()<<"  "<<tmpunc.Phi()<<"   "<<Z.Mod()<<"   "<<Z.Phi()<<" ---> "<<itUM->first+"_"+iter->first<<"  "<<uParasUnc[ itUM->first+"_"+iter->first ]+ZPt<<endl;
	  }

	}

    }
    
    ZData->Fill();
  }
  cout<<" done "<<ZData->GetEntries()<<endl;
   
  if(fillH) {
    TFile* oFile=new TFile( (dat?"data.root":"MC.root"), "RECREATE");
    for(iter=metPts.begin();iter!=metPts.end();iter++) {
      cout<<iter->first<<endl;
      histoPara[ iter->first ]->Write();
      histoPerp[ iter->first ]->Write();
    }
    oFile->Close();
  }

  loaded = true;
  return ZData;
  
}



//key=para/perp w,w/o error, object = vtx, pt, value 
recMap GetSampleStat(string comp, string metType, bool AllVtx, bool AllPts, bool resp=false) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  _resp=resp;
  bool debug_=false;
  _sumEt = AllVtx && AllPts;
  if(!loaded)
    Data = LoadTree();

  cout<<" Getting sample :: component = "<<comp<<" for "<<metType
      <<"  allvtx="<<AllVtx<<" allpts="<<AllPts<<"  "<<dat<<endl;

  cout<<" number of entries "<<Data->GetEntries()<<endl;

  TTree* ZData = (TTree*)Data->Clone();
  cout<<" Tree cloned "<<endl;
  //int nBin = 13;
  //double BinVar[17] = {0,5,10,20,30,40,50,60,70,80,90,100,110,120,140,160,200};
  //double BinVar[14] = {0.,10.,20.,30.,40.,60.,80.,100.,120.,140.,170.,200.,250.,300.};
  //double BinVar[35] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,230,260,300};
  // double BinVar[24] = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,230,260,300};
  
  int nBin=17;
  double bVar[18] = {0,5,10,15,20,25,30,40,60,85,100,120,145,165,200,250,300,400};//,600};
  double *BinVar(0);
  //  int nBin=19;
  //double BinVar[20] = {0,5,10,15,20,25,30,40,60,85,100,120,145,165,200,250,300,350,400,600};

  //int mType[16] = {100,20,24,25,26,27,21,2,2,2,2,2,2,2,2,2};
  //int colors[16] = {1000,kOrange+9,1,38,896,814,kMagenta+3,2,3,4,5,6,7,8,9,11};
  
  // nBin=1;
  // double BinVar[2] = {10,20};

  int nVtxMax=35;

  int offset=1;
  int NVMax = nVtxMax;
  
  //bool sumEt=false;

  if(AllVtx) {
    offset=0;
    NVMax=1;
  }

  if(AllPts) {
    nBin=1;
  }

  if(AllVtx && AllPts) {
    offset=0;
    NVMax=1;
    nBin=23;
    //sumEt=true;
    BinVar=new double[24];
    for(int i=0;i<24;i++)
      BinVar[i]= 100*i;
  }
  else {
    BinVar=new double[19];
    for(int i=0;i<19;i++) BinVar[i]=bVar[i];
  }

  //prepare output ====================================
  map<string, vector<vector<float> > > oMap;
  vector<vector<float> > cqt(NVMax,vector<float>(nBin,0));
  vector<vector<float> > cqte(NVMax,vector<float>(nBin,0));
  vector<vector<float> > upx(NVMax,vector<float>(nBin,0));
  vector<vector<float> > upxe(NVMax,vector<float>(nBin,0));
  vector<vector<float> > upxu(NVMax,vector<float>(nBin,0));

  // ====================================================
  map<string, TH1* > uncHistosU;
 

  for(int Nvtx=offset;Nvtx<NVMax;Nvtx++) {
      RooRealVar weight("weight","weight",-10,100);
      RooRealVar Vtxcat("nVtx","nVtx", AllVtx?-0.1:(Nvtx - 0.1) , AllVtx?100.:(Nvtx+0.1) ); //NVertex have to be changed 
      //with the name of the
      //corresponding tree branch

      for(int iqt=0;iqt<nBin;iqt++) {
	cout<<" Start bin "<<iqt<<"   "<<BinVar[iqt]<<"    "<<BinVar[iqt+1]<<"   "<<Nvtx<<"   with sumET?"<<_sumEt<<"   "<<comp<<endl;

	string n=_sumEt?((metType+"_sumEt").c_str()):"Zpt"; //cout<<n<<"  "<<AllVtx<<"  "<<AllPts<<"  "<<sumEt<<endl;
	RooRealVar qTcat(n.c_str(),n.c_str(),(AllPts && !_sumEt)?0.:BinVar[iqt],(AllPts && !_sumEt)?1000.:BinVar[iqt+1]);
	//RooRealVar qTcat("Zpt","Zpt",(AllPts)?0.:BinVar[iqt],(AllPts)?1000.:BinVar[iqt+1]);
	double m,M,c;
	if(comp=="para")
	  {
	    c = (BinVar[iqt]+BinVar[iqt+1])/2.;
	    m = AllPts? -400 : -500;//-400+(-10*iqt)+5; // -c-100;
	    M = AllPts? 400  : 75;//(-10*iqt)+400+5; //-c+100
	  }
	else {
	  //	c = 0;
	  m = -75; //FIXME 100
	  M = 75;
	}
	//cout<<iqt<<"   "<<BinVar[iqt]<<"    "<<BinVar[iqt+1]<<"   "<<m<<"    "<<M<<endl;
	RooRealVar upara((metType+"_u"+comp).c_str(),"upara",m,M);
	RooRealVar resoP( (metType+"_redupara").c_str(),"resoPar", -75,75 ); //100

	//by pass for debugging
	if(debug_) {
	  cqt[Nvtx-offset][iqt] = (BinVar[iqt]+BinVar[iqt+1])/2.;
	  cqte[Nvtx-offset][iqt] = 0;

	  upx[Nvtx-offset][iqt] = Nvtx;
	  upxe[Nvtx-offset][iqt] = 0;
	  continue;
	}

	RooDataSet* Dpara;
	if(dat)
	  { 
	    if(!_sumEt)
	      Dpara = new RooDataSet("data","Data",ZData, RooArgSet(Vtxcat,qTcat,upara,resoP) );
	    else
	      Dpara = new RooDataSet("data","Data",ZData, RooArgSet(Vtxcat,qTcat,upara,weight,resoP),"","weight" );
	  }
	else
	  Dpara = new RooDataSet("data","Data",ZData, RooArgSet(Vtxcat,qTcat,upara,weight,resoP),"","weight" );
	
	if(computeUnc && !dat && !noPU) {
	  for(size_t iu=0;iu<metUncNames.size();iu++) {
	    delete uncHistosU[ metUncNames[iu] ];
	    delete uncHistosU[ "central" ];
	    RooRealVar uncU((metUncNames[iu]+"_"+metType+"_u"+comp).c_str(),"upara",m,M);
	    RooRealVar uncRU( (metUncNames[iu]+"_"+metType+"_redupara").c_str(),"resoPar", -75,75 );
	    RooDataSet* uncDS = new RooDataSet("data","Data",ZData, RooArgSet(Vtxcat,qTcat,uncU,weight,uncRU),"","weight" );
	    // cout<<uncDS->sumEntries()<<"   "<<uncDS->numEntries()<<endl;
	    if(comp=="para" && !resp) {
	      uncHistosU[ metUncNames[iu] ] = uncDS->createHistogram( (metUncNames[iu]+"_"+metType+"_redupara").c_str(), 75 );
	      uncHistosU[ "central" ] = Dpara->createHistogram( (metType+"_redupara").c_str(), 75 );
	    }
	    else if(resp) {
	      uncHistosU[ metUncNames[iu] ] = uncDS->createHistogram( (metUncNames[iu]+"_"+metType+"_upara").c_str(), 75 );
	      uncHistosU[ "central" ] = Dpara->createHistogram( (metType+"_upara").c_str(), 75 );
	    }
	    else {
	      uncHistosU[ metUncNames[iu] ] = uncDS->createHistogram( (metUncNames[iu]+"_"+metType+"_u"+comp).c_str(), 75 );
	      uncHistosU[ "central" ] = Dpara->createHistogram( (metType+"_u"+comp).c_str(), 75 );
	    }
	   delete uncDS;
	  }//unc loop
	}

	//cout<<" gloubi "<<endl;
	RooRealVar* meanqT = (*Dpara).meanVar(qTcat);
	RooRealVar* rms = (*Dpara).rmsVar(upara);
	RooRealVar* mean = (*Dpara).meanVar(upara);
	// if(comp=="para" && !resp)
	//   rms = (*Dpara).rmsVar(resoP);

	
	double fwhm=-1, efwhm=-1;
	double peak=-1, epeak=-1;

	if(fit) {
	  double qt= 0;//-1* ( BinVar[iqt+1] + BinVar[iqt] )/2.;
	  double qtm = -30;//qt -100;
	  double qtM = 30;//qt +100;
	  // double qt= -1* ( BinVar[iqt+1] + BinVar[iqt] )/2.;
	  // double qtm = qt -100;
	  // double qtM = qt +100;

	  if(comp=="perp") {
	    qt =0;
	    qtm = -20; //-100
	    qtM = 20; //100
	  }

	  // if(Dpara->sumEntries() < 0.001 || Dpara->numEntries() ==0 ) 
	  //   cout<< " No Entry for bin ["<<m<<":"<<M<<"] / ["<<(AllVtx?0:Nvtx)<<":"<<(AllVtx?100.:Nvtx)<<"]"<<endl;
      
	  RooRealVar g_w("g_w","width Gaus", 10., 0., 40., "GeV");
	  RooRealVar gamma_Z0( "gamma_Z0", "Z0 width",2.3, 0, 10, "GeV" );
	  RooRealVar v_m("v_m","v_m",qt, qtm , qtM , "GeV" );
	  RooVoigtian* voigt= new RooVoigtian("voigt","Voightian",upara, v_m, gamma_Z0, g_w);
	  if(comp=="para")
	    voigt= new RooVoigtian("voigt","Voightian",resoP, v_m, gamma_Z0, g_w);

	  RooFitResult* result = voigt->fitTo( (*Dpara) ,RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );
	  //  cout<<qt<<"    "<<Nvtx<<" ==================> The status "<<result->status()<<endl;
	  
	  // ostringstream os;
	  // os<<iqt;
	  // TCanvas* c1 =new TCanvas("c1","ZPeak",600,600);
	  // RooPlot* frame = upara.frame() ;
	  // if(comp=="para")
	  //   frame = resoP.frame() ;
	  // Dpara->plotOn(frame) ;
	  // voigt->plotOn(frame) ;
	  // frame->Draw();
	  // c1->SaveAs( ("projs/test_"+os.str()+(dat?"data":"MC")+".root").c_str() );

	  //Get the FWHM
	  double sigma = g_w.getVal();
	  double gamma = gamma_Z0.getVal();
	  double esigma = g_w.getError();
	  double egamma = gamma_Z0.getError();

	  double Vsg = result->correlation(g_w,gamma_Z0) ;
	  double Vgs = result->correlation(gamma_Z0,g_w) ;
	  double Vss = result->correlation(g_w,g_w) ;
	  double Vgg = result->correlation(gamma_Z0,gamma_Z0) ;

	  fwhm = FWHM( sigma, gamma );
	  efwhm = FWHMError(/*sigma, gamma,*/ esigma, egamma,
			    Vss, Vsg, Vgs, Vgg);

	  peak = v_m.getVal();
	  epeak = v_m.getError();
	  delete voigt;
	}
	
	if( (*Dpara).sumEntries()>20) {
	
	  cqt[Nvtx-offset][iqt] = meanqT->getVal();
	  cqte[Nvtx-offset][iqt] = meanqT->getError();

	  if(!fit) {
	    if(dat)
	      upx[Nvtx-offset][iqt] = resp?mean->getVal():rms->getVal();// /(0.985 + 0.000387* meanqT->getVal() ); //FIXME
	    else
	      upx[Nvtx-offset][iqt] = resp?mean->getVal():rms->getVal();
	    upxe[Nvtx-offset][iqt] = resp?mean->getError():rms->getError ();

	    if(isnan(rms->getVal())) {
	      TH1* ht(0);
	      ht = Dpara->createHistogram( (metType+"_u"+comp).c_str(),100);
	      if(comp=="para" && !resp)
	      	ht = Dpara->createHistogram( (metType+"_redupara").c_str(),100);
	      TH1* ht2 = Dpara->createHistogram("Zpt",100);
	      TFile* tr=new TFile("fr.root","RECREATE");
	      ht->Write();
	      ht2->Write();
	      tr->Close();
	      cout<<rms->getVal()<<"   "<<rms->getError()<<"   "<<mean->getVal()<<"  "<<ht->GetRMS()<<"   "<<ht2->GetMean()<<endl;

	      cqt[Nvtx-offset][iqt] = ht2->GetMean();
	      cqte[Nvtx-offset][iqt] =ht2->GetMeanError();

	      upx[Nvtx-offset][iqt] =resp?ht->GetMean():ht->GetRMS();
	      upxe[Nvtx-offset][iqt] =resp?ht->GetMeanError():ht->GetRMSError();

	    }

	    cout<<" ---> "<<mean->getVal()<<"  "<<mean->getError()<<"  "<<rms->getVal()<<"  "<<Dpara->numEntries()<<"  "<<(*Dpara).sumEntries()<<"   "<<meanqT->getVal()<<endl;
	  }
	  else {
	    upx[Nvtx-offset][iqt] = resp?peak:fwhm/2.3548;
	    upxe[Nvtx-offset][iqt] = resp?epeak:efwhm/2.3548;
	  }
	  
	  if(computeUnc && !dat && !noPU) {
	    float unc = ComputeUncertainties(uncHistosU);
	    upxu[Nvtx-offset][iqt] = resp?unc:unc;
	    //upx[Nvtx-offset][iqt] = upx[Nvtx-offset][iqt]-unc;
	  }	  
	  

	  }// number of event > 20
	//cout<<Dpara<<endl;
	delete Dpara;
      }//qt

  }//Nvtx
  
  
  //filling output ==============
  oMap [ "centerQt" ]     = cqt;
  oMap [ "centerQtErr" ]  = cqte;
  oMap [ "u"+comp ]       = upx;
  oMap [ "u"+comp+"Err" ] = upxe;

  if(computeUnc && !dat && !noPU) {
    oMap [ "u"+comp+"Unc" ] = upxu;
  }

  
  delete ZData;
  return oMap;
  
}


vector<TGraph*> computeResponse(string metType) {

 
  recMap oMap = GetSampleStat("para",metType, true, false, true );

  int NptBin = (oMap[ "centerQt" ])[0].size();

  vector<TGraph*> vg;

  TGraphAsymmErrors* og=new TGraphAsymmErrors();
  TGraphAsymmErrors* ogunc=new TGraphAsymmErrors();

  double bVar[18] = {0,5,10,15,20,25,30,40,60,85,100,120,145,165,200,250,300,400};//,600};

  int n=0;
  for(int i=0;i<NptBin;i++) {

    float x = (oMap[ "centerQt" ])[0][i];
    float y = -(oMap[ "upara" ])[0][i]/x;


   
    
    float ex = (oMap[ "centerQtErr" ])[0][i];
    float ey = (oMap[ "uparaErr" ])[0][i]/x;
    
    cout<<x<<"  "<<(oMap[ "centerQtErr" ])[0][i]<<"  "
 	<<(oMap[ "upara" ])[0][i]<<"  "
 	<<(oMap[ "uparaErr" ])[0][i]<<" ---> "<<y<<"   "<<ey<<endl;


    if(x!=0) {
      og->SetPoint(n, x, y );
      og->SetPointError(n, ex,ex,ey, ey );
    
      if(computeUnc && !dat && !noPU) {
	float uy = (oMap[ "uparaUnc" ])[0][i];
	ogunc->SetPoint(n,x,y);

	float exl=(bVar[i+1]-bVar[i])/2.;
	float exh=((i+1)!=NptBin)?((bVar[i+2]-bVar[i+1])/2.):(400);
	//ogunc->SetPointError(n,0,uy);
	ogunc->SetPointError(n,exl,exh,uy/x,uy/x);
	// if((i+1)==NptBin) {
	//   ogunc->SetPoint(n+1, 400, y );
	//   ogunc->SetPointError(n+1, 350/2.,0.,uy/x,uy/x );
	// }
	cout<<x<<" --> "<<exl<<"   "<<exh<<"  "<<uy<<endl;
      }

      n++;
    }
  }

  vg.push_back((TGraph*)og);
  if(computeUnc && !dat && !noPU)
    vg.push_back( (TGraph*)ogunc);

 if(computeUnc && !dat && !noPU)
    vg[1]->SetName("unc");

  return vg;

}


vector<TGraph*> computeResolution(string comp, string metType, bool AllVtx, bool AllQts ) {

  
  recMap oMap = GetSampleStat(comp,metType, AllVtx, AllQts );

  int NptBin = (oMap[ "centerQt" ])[0].size();
  int NvtxBin = (oMap[ "centerQt" ]).size();

  if(AllVtx && AllQts) { cout<<" Don't use AllVtx and AllQts options in the same time except if you want to use sumEt "<<endl; }

  // if(AllVtx) AllQts=false;
  // if(AllQts) AllVtx=false;

  vector<TGraph*> vg;
  double bVar[18] = {0,5,10,15,20,25,30,40,60,85,100,120,145,165,200,250,300,400};//,600
  //double Err[18]={2.5,2.5,2.5,2.5,2.5,2.5,3.75,7.5,8.75,7.5,8.75,
  if(AllVtx) {
    TGraphErrors* og=new TGraphErrors();
    TGraphAsymmErrors* ogunc=new TGraphAsymmErrors();
    
    //og->SetName();
    
    int n=0;
    for(int i=0;i<NptBin;i++) {
      
      float x = (oMap[ "centerQt" ])[0][i];
      float y = (oMap[ "u"+comp ])[0][i];
      
      float ex = (oMap[ "centerQtErr" ])[0][i];
      float ey = (oMap[ "u"+comp+"Err" ])[0][i];
    
      
      if(x!=0) {
	og->SetPoint(n, x, y );
	og->SetPointError(n, ex, ey );
	
	if(computeUnc && !dat && !noPU) {
	  float uy = (oMap[ "u"+comp+"Unc" ])[0][i];
	  ogunc->SetPoint(n,x,y);

	  float exl=(bVar[i]+bVar[i+1])/2.;
	  float exh=((i+1)!=NptBin)?((bVar[i+1]+bVar[i+2])/2.):(400);
	  //ogunc->SetPointError(n,0,uy);
	  ogunc->SetPointError(n,exl,exh,uy,uy);
	  // if((i+1)==NptBin) {
	  //   ogunc->SetPoint(n+1, 400, y );
	  //   ogunc->SetPointError(n+1, 350/2.,0.,uy/x,uy/x );
	  // }
	}

	n++;
      }
    }//pt bins

    vg.push_back((TGraph*)og);
    if(computeUnc && !dat && !noPU)
      vg.push_back( (TGraph*)ogunc);

  } //if allVtxs
  else if(AllQts) {
    TGraphAsymmErrors* og=new TGraphAsymmErrors();
    TGraphAsymmErrors* ogunc=new TGraphAsymmErrors();
    //og->SetName();
    
    int n=0;
    for(int i=0;i<NvtxBin;i++) {
      
      float x = i;//(oMap[ "centerQt" ])[i][0];
      float y = (oMap[ "u"+comp ])[i][0];
      
      float ex = 0;//(oMap[ "centerQtErr" ])[i][0];
      float ey = (oMap[ "u"+comp+"Err" ])[i][0];
      
      if(x!=0 && fabs(y)>0.0001) {
	og->SetPoint(n, x, y );
	og->SetPointError(n, ex,ex,ey, ey );

	if(computeUnc && !dat && !noPU) {
	  float uy = (oMap[ "u"+comp+"Unc" ])[i][0];
	  ogunc->SetPoint(n,x,y);
	  //ogunc->SetPointError(n,0.5,uy);
	  ogunc->SetPointError(n,0.5,0.5,uy,uy);
	}	  

	n++;
      }
    }//nvtx bins

    vg.push_back((TGraph*)og);
    if(computeUnc && !dat && !noPU) 
      vg.push_back((TGraph*)ogunc);

  } //if allqts
  else {

    for(int i=0;i<NvtxBin;i++) {
      TGraphErrors* og=new TGraphErrors();
      TGraphAsymmErrors* ogunc=new TGraphAsymmErrors();
      
      int n=0;
      for(int j=0;j<NptBin;j++) {
	
	float x = (oMap[ "centerQt" ])[i][j];
	float y = (oMap[ "u"+comp ])[i][j];
      
	float ex = (oMap[ "centerQtErr" ])[i][j];
	float ey = (oMap[ "u"+comp+"Err" ])[i][j];
      
	if(x!=0 && fabs(y)>0.0001) {
	  og->SetPoint(n, x, y );
	  og->SetPointError(n, ex, ey );
	
	  if(computeUnc && !dat && !noPU) {
	    float uy = (oMap[ "u"+comp+"Unc" ])[i][j];

	    ogunc->SetPoint(n,x,y);
	    ogunc->SetPointError(n,0,0,uy,uy);
	  }

	  n++;
	}
      }//pt bins
      
      vg.push_back((TGraph*)og);
      if(computeUnc && !dat && !noPU)
	vg.push_back((TGraph*)ogunc);
      
    }//vtx bins
  
  }//condition
  
  //uncertainties
  if(computeUnc && !dat && !noPU)
    vg[1]->SetName("unc");


  //by pass for sumET
  if(AllQts && AllVtx) {

  }

  return vg;

}

void CompaRespCurves(bool unc) {

  RooMsgService::instance().setGlobalKillBelow(ERROR);

  ConfigAnalysis( 0,0, 0, unc, 1);
 
  vector<TGraphAsymmErrors*> graphs;
  TGraphAsymmErrors* uncGraph(0);

  TLegend*leg=new TLegend(0.5,0.5,0.7,0.7);
  leg->SetFillColor(0); leg->SetLineColor(0);
  leg->SetShadowColor(0);

  for(size_t i=0;i<metNames.size();i++) {

    //FIXME
    loaded=false;
    if(i%2==1) {dat=true;}
    else dat=false;
    vector<TGraph*> tmpgs = computeResponse(metNames[i]);
    graphs.push_back( (TGraphAsymmErrors*)tmpgs[0] );
  
    graphs.back()->SetLineColor(colors[i/2]);
    graphs.back()->SetMarkerColor(colors[i/2]);
    graphs.back()->SetMarkerStyle(mtypes[i]);
    
    if(i%2==1) graphs.back()->SetName("data");
    else  {
      graphs.back()->SetName("mc");
      cout<<" ============= uncertainty graph"<<endl;
      cout<<tmpgs[1]<<endl;
      if(computeUnc && !noPU)
	uncGraph = (TGraphAsymmErrors*)tmpgs[1];
      if(computeUnc && noPU) 
	uncGraph= ComputeChristianUncertainties("resp","",tmpgs[0]);
     
      if(sumUnc)
	quadSumUnc(graphs.back(), uncGraph );

    }
    // if(i%2==0)
    //   leg->AddEntry(graphs.back(), (metNames[i]+" data").c_str(), "pl");
    // else 
    //   leg->AddEntry(graphs.back(), (metNames[i]+" MC").c_str(), "pl");

    map<string,string>::const_iterator iter;
    iter = nameMap.find(metNames[i]);
    if(iter!=nameMap.end()) {
      if(i%2==1)
	leg->AddEntry(graphs.back(), (iter->second+" #slash{E}_{T} data").c_str(), "pl");
      else 
	leg->AddEntry(graphs.back(), (iter->second+" #slash{E}_{T} MC").c_str(), "pl");
    }
  //leg->AddEntry(graphs.back(), (iter->second+" #slash{E}_{T}").c_str(), "pl");

    // if(iter->second=="pfType1") {
    //   graphs.back()->SetLineColor(kOrange);
    //   graphs.back()->SetMarkerColor(kOrange);
    //   graphs.back()->SetMarkerStyle(25);
    // }

 }
  
  
  TCanvas* c=NULL;
  vector<TPad*>pads = PreparePads(c);
  
  pads[0]->cd();
  

  TH1F* histo= new TH1F("h","h",150,0,400);

  histo->GetXaxis()->SetTitle("Z q_{T} [GeV]");
  histo->GetYaxis()->SetTitle("-<u_{||}>/q_{T}");
  histo->GetYaxis()->SetRangeUser(0,1.1);

  histo->Draw();
  
  //uncertainty
  if(computeUnc) {
    uncGraph->SetMarkerStyle(1);
    uncGraph->SetMarkerColor(0);
    uncGraph->SetLineColor(0);
    uncGraph->SetFillColor(kGray+1);
    uncGraph->SetFillStyle(3001);
    uncGraph->Draw("P E3 same");
  }

  for(size_t i=0;i<graphs.size();i++) {
    graphs[i]->Draw("p");
  }

  leg->Draw("same");

  TLatex lat;
  lat.SetTextSize(0.04);
  lat.SetNDC();
  lat.DrawLatex(0.597,0.9653,"CMS Preliminary 2012");
  lat.DrawLatex(0.145,0.9598,"12.2 fb^{-1} at #sqrt{8} TeV");

  pads[1]->cd();
  TH1F* histo2=(TH1F*)histo->Clone();
  histo2->GetYaxis()->SetRangeUser(0.7,1.3);
  histo2->GetYaxis()->SetRangeUser(0.7,1.3);
  histo2->GetYaxis()->SetNdivisions(4,5,0);
  histo2->GetYaxis()->SetLabelSize(0.13);
  histo2->GetYaxis()->SetTitle("Data/MC");
  histo2->GetYaxis()->SetTitleSize(0.19);
  histo2->GetYaxis()->SetTitleOffset(0.32);
  histo2->GetXaxis()->SetLabelSize(0.20);
  histo2->GetXaxis()->SetTitleSize(0.21);
  histo2->Draw();

   // and uncertainties
  if(computeUnc) {
  TGraphAsymmErrors* ratioUnc=new TGraphAsymmErrors();
  ratioUnc->SetName("ratioUnc");
  double x,y,exl,eyl,exh,eyh;
  double xd,yd;
  int n=0;
  for(int i=0;i<uncGraph->GetN();i++) {

    uncGraph->GetPoint(i,x,y);
    if(i<graphs[1]->GetN()) {
      graphs[1]->GetPoint(i,xd,yd);
    }
    
    exh = uncGraph->GetErrorXhigh(i);
    exl = uncGraph->GetErrorXlow(i);
    eyh = uncGraph->GetErrorYhigh(i);
    eyl = uncGraph->GetErrorYlow(i);
    
    if( fabs(x-xd)<xd/10. || i==uncGraph->GetN()-1 ) {
      cout<<x<<"  "<<xd<<endl;
      ratioUnc->SetPoint(n,x,1.);
      ratioUnc->SetPointError(n,exl,exh,eyl/y*yd/y,eyh/y*yd/y);
      n++;
    }
  }

  ratioUnc->SetMarkerStyle(1);
  ratioUnc->SetMarkerColor(0);
  ratioUnc->SetLineColor(0);
  ratioUnc->SetFillColor(kGray+2);
  ratioUnc->SetFillStyle(3002);
  ratioUnc->Draw("p E3 same");
  }  

  TLine* line=new TLine(0,1,400,1);
  line->SetLineColor(kGray+2);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  

  vector<TGraphAsymmErrors*> ratios = getRatios(graphs);
  for(size_t ig=0;ig<ratios.size();ig++) {
    ratios[ig]->Draw("P");
  }
  
  //and save the plots
  string base = "pf";
  if(metNames[0].find("MVA")!=(size_t)-1)
    base = "mva";
  if(metNames[0].find("NoPileUp")!=(size_t)-1)
    base = "noPU";
  if(metNames[0].find("Calo")!=(size_t)-1)
    base = "calo";

  string name=base+"responseVsQt"+(fit?"_doFit":"");
  
  c->SaveAs( ("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/Plots/"+name+".png").c_str() );
  c->SaveAs( ("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/Plots/"+name+".pdf").c_str() );
  c->SaveAs( ("/afs/cern.ch/user/m/mmarionn/workspace/private/METStudies/Plots/"+name+".root").c_str() );
  //c->SaveAs( ("Plots/"+name+".png").c_str() );





}



void CompaResoCurves(string comp, bool vsVtx, bool Fit=false, bool doqt=false, bool unc=false, bool sumEt=false) {

  RooMsgService::instance().setGlobalKillBelow(ERROR);

  _sumEt = sumEt;

  if(vsVtx || sumEt) noPU=false;

  ConfigAnalysis( 0,Fit, doqt, unc, 1); //FIXME ttbar

  vector<TGraphAsymmErrors*> graphs;
  TGraphAsymmErrors* uncGraph(0);
  
  TLegend*leg=new TLegend(0.54,0.16,0.84,0.403);
  leg->SetFillColor(0); leg->SetLineColor(0);
  leg->SetShadowColor(0);

  for(size_t i=0;i<metNames.size();i++) {
    
    //FIXME
    loaded=false;
    if(i%2==1) {dat=true;}
    else dat=false;
    vector<TGraph*> tmpgs = computeResolution(comp,metNames[i], !vsVtx || sumEt, vsVtx||sumEt);
    graphs.push_back( (TGraphAsymmErrors*)tmpgs[0] );
    
    graphs.back()->SetLineColor(colors[i/2]);
    graphs.back()->SetMarkerColor(colors[i/2]);
    graphs.back()->SetMarkerStyle(mtypes[i]);
    
    if(i%2==1) graphs.back()->SetName("data");
    else {
      graphs.back()->SetName("mc");
      cout<<" ============= uncertainty graph"<<endl;
      cout<<tmpgs[1]<<endl;
      if(computeUnc && !noPU)
	uncGraph = (TGraphAsymmErrors*)tmpgs[1];
      if(computeUnc && noPU) 
	uncGraph= ComputeChristianUncertainties("resoQt",comp,tmpgs[0]);
	
      if(sumUnc)
	quadSumUnc(graphs.back(), uncGraph );

    }

    // if(i%2==0)
    //   leg->AddEntry(graphs.back(), (metNames[i]+" data").c_str(), "pl");
    // else
    //   leg->AddEntry(graphs.back(), (metNames[i]+" MC").c_str(), "pl");
    
    map<string,string>::const_iterator iter;
    iter = nameMap.find(metNames[i]);
    if(iter!=nameMap.end()) {
      if(i%2==1)
	leg->AddEntry(graphs.back(), (iter->second+" #slash{E}_{T} data").c_str(), "pl");
      else 
	leg->AddEntry(graphs.back(), (iter->second+" #slash{E}_{T} MC").c_str(), "pl");
    }


    //   leg->AddEntry(graphs.back(), (iter->second+" #slash{E}_{T}").c_str(), "pl");

    // if(iter->second=="pfType1") {
    //   graphs.back()->SetLineColor(kOrange);
    //   graphs.back()->SetMarkerColor(kOrange);
    //   graphs.back()->SetMarkerStyle(25);
    // }
    
  }
  
  TCanvas* c=NULL;
  vector<TPad*>pads = PreparePads(c);
  
  pads[0]->cd();
  

  TH1F* histo;
  if(!vsVtx && !sumEt)
    histo= new TH1F("hr","hr",150,0,400);
  else if(sumEt)
    histo= new TH1F("hr","hr",30,0,3000);
  else
    histo= new TH1F("hr","hr",35,0,35);

  if(!vsVtx)
    histo->GetXaxis()->SetTitle("Z q_{T} [GeV]");
  else
    histo->GetXaxis()->SetTitle("NPV  ");

  if(comp=="para")
    histo->GetYaxis()->SetTitle( Fit?"#sigma(u_{||}) [GeV] ":"RMS(u_{||}) [GeV] ");
  else
    histo->GetYaxis()->SetTitle( Fit?"#sigma(u_{#perp}  ) [GeV] ":"RMS(u_{#perp}  ) [GeV] ");  

  if(!vsVtx)
    histo->GetYaxis()->SetRangeUser(1.,40.);
  else {
    histo->GetYaxis()->SetRangeUser(1.,30.);
    histo->GetYaxis()->SetTitleOffset(0.9);
  }

  histo->Draw();
  
  //uncertainty
  if(computeUnc) {
    uncGraph->SetMarkerStyle(1);
    uncGraph->SetMarkerColor(0);
    uncGraph->SetLineColor(0);
    uncGraph->SetFillColor(kGray+2);
    uncGraph->SetFillStyle(3002);
    uncGraph->Draw("P E3 same");
  }

  for(size_t i=0;i<graphs.size();i++) {
    graphs[i]->Draw("p same");
  }
  
  leg->Draw("same");


  TLatex lat;
  lat.SetTextSize(0.04);
  lat.SetNDC();
  lat.DrawLatex(0.597,0.9653,"CMS Preliminary 2012");
  lat.DrawLatex(0.636,0.878,"12.2 fb^{-1} at #sqrt{8} TeV");
  

  pads[1]->cd();
  TH1F* histo2=(TH1F*)histo->Clone();
  histo2->GetYaxis()->SetRangeUser(0.7,1.3);
  histo2->GetYaxis()->SetNdivisions(4,5,0);
  histo2->GetYaxis()->SetLabelSize(0.13);
  histo2->GetYaxis()->SetTitle("Data/MC");
  histo2->GetYaxis()->SetTitleSize(0.19);
  histo2->GetYaxis()->SetTitleOffset(0.32);
  histo2->GetXaxis()->SetLabelSize(0.20);
  histo2->GetXaxis()->SetTitleSize(0.21);
  histo2->Draw();
  
  //  TPolyLine* l = ne

  // and uncertainties
 if(computeUnc) {
  TGraphAsymmErrors* ratioUnc=new TGraphAsymmErrors();
  double x,y,exl,eyl,exh,eyh;
  double xd,yd;
  int n=0;
  for(int i=0;i<uncGraph->GetN();i++) {

    uncGraph->GetPoint(i,x,y);
    if(i<graphs[1]->GetN()) {
      graphs[1]->GetPoint(i,xd,yd);
    }
    
    exh = uncGraph->GetErrorXhigh(i);
    exl = uncGraph->GetErrorXlow(i);
    eyh = uncGraph->GetErrorYhigh(i);
    eyl = uncGraph->GetErrorYlow(i);
    
    if( fabs(x-xd)< xd/5.|| i==uncGraph->GetN()-1) {
      cout<<x<<"  "<<xd<<endl;
      ratioUnc->SetPoint(n,x,1.);
      ratioUnc->SetPointError(n,exl,exh,eyl/y*yd/y,eyh/y*yd/y);
      n++;
    }
  }

  ratioUnc->SetMarkerStyle(1);
  ratioUnc->SetMarkerColor(0);
  ratioUnc->SetLineColor(0);
  ratioUnc->SetFillColor(kGray+2);
  ratioUnc->SetFillStyle(3002);
  ratioUnc->Draw("p E3 same");
 }

  // and the end now
  TLine* line=new TLine(0,1,vsVtx?35:400,1);
  line->SetLineColor(kGray+2);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();
  
  vector<TGraphAsymmErrors*> ratios = getRatios(graphs);
  for(size_t ig=0;ig<ratios.size();ig++) {
    ratios[ig]->Draw("P");
  }
  

  //and save the plots
  string base = "pf";
  if(metNames[0].find("MVA")!=(size_t)-1)
    base = "mva";
  if(metNames[0].find("NoPileUp")!=(size_t)-1)
    base = "noPU";
  if(metNames[0].find("Calo")!=(size_t)-1)
    base = "calo";

  string name=base+"rms"+((comp=="para")?"Upara":"Uperp")
    +"Vs"+(vsVtx?"NPV":(sumEt?"SumEt":"Qt"))+(doqt?"_W":"")+(fit?"_doFit":"");
  
  cout<<name<<"   "<<endl;
  cout<<c<<endl;

  c->SaveAs( ("Plots/"+name+".png").c_str() );
  c->SaveAs( ("Plots/"+name+".pdf").c_str() );
  c->SaveAs( ("Plots/"+name+".root").c_str() );
  //c->SaveAs( ("Plots/"+name+".png").c_str() );
  
}



void drawDistribution(string var, bool b, int NB=100, 
		      float mB=0,float MB=400) {

  ConfigAnalysis(0,0,0,0,0 );
  

  TTree* D;
  if(!b) {
    if(!loaded)
      D = LoadTree();
    else
      D = Data;
  }
  else
    D = LoadTree(b);

  vector<TH1F*> vh;
  TLegend*leg=new TLegend(0.5,0.5,0.7,0.7);
  leg->SetFillColor(0); leg->SetLineColor(0); leg->SetShadowColor(0);
  for(size_t i=0;i<metNames.size();i++) {

    TString sn="name_"; sn+=i;
    TH1F* htmp=new TH1F(sn,sn,NB,mB,MB);
    htmp->SetLineColor(colors[i]);
    htmp->SetLineWidth(2);

    
    if(var=="met") {
      D->Draw( (metNames[i]+"_pt>>"+(string)sn).c_str(),"","goff");
      htmp->GetXaxis()->SetTitle("#slash{E}_{T} [GeV] ");
      htmp->GetYaxis()->SetTitle("1/N#times#partialN/#partial#slash{E}_{T} ");
    }
    else
      D->Draw( (metNames[i]+"_"+var+">>"+(string)sn).c_str(),"","goff");

    htmp->Scale(1./htmp->Integral());

    map<string,string>::const_iterator iter;
    iter = nameMap.find(metNames[i]);
    if(iter!=nameMap.end())
      leg->AddEntry(htmp, (iter->second+" #slash{E}_{T}").c_str(), "l");
    else
      leg->AddEntry(htmp,(metNames[i]).c_str(),"l");
    vh.push_back(htmp);

    if(iter->second=="pfType1") {
      htmp->SetLineColor(kOrange);
      htmp->SetMarkerColor(kOrange);
      htmp->SetMarkerStyle(25);
    }

  }
  
  TCanvas* c = new TCanvas("cc","cc");
  c->cd();
  for(size_t i=0;i<vh.size();i++) {
    if( i==0 ) vh[i]->Draw();
    else vh[i]->Draw("same");
  }

  leg->Draw("same");


}




TVector2 
phiCorrection(TVector2 met, int Nvtx, bool isD, string type) {

  float corX,corY;
 
 if(isD) {
   
    //CHristian corrections 2012ABC

    if(type.find("mva")!=(size_t)-1 || type.find("MVA")!=(size_t)-1) {
      corX = 0.2792 + 0.0530*Nvtx; 
      corY = -0.1661 -0.0223*Nvtx;
    }
    else if(type.find("noPU")!=(size_t)-1){
      corX = 0.1554 + 0.1361*Nvtx; 
      corY = -0.2539 -0.0943 *Nvtx;
    }
    else if(type.find("calo")!=(size_t)-1 || type.find("Calo")!=(size_t)-1){
      corX =  -0.7622 + 0.2308*Nvtx; 
      corY =  0.3364 -0.1221 *Nvtx;
    }
    else {
      corX = 0.2661 + 0.3217*Nvtx; 
      corY = -0.2251 -0.1747*Nvtx;
    }
  }
  else {
    
   if(type.find("mva")!=(size_t)-1) {
      corX = -0.3009 -0.0251*Nvtx; 
      corY = -0.2129 -0.0261 *Nvtx;
    }
    else if(type.find("noPU")!=(size_t)-1){
      corX = 0.1900  - 0.0030*Nvtx; 
      corY = -0.2539 - 0.0943*Nvtx;
    }
    else if(type.find("calo")!=(size_t)-1 || type.find("Calo")!=(size_t)-1){
      corX = 0.3120 -0.2755 *Nvtx; 
      corY = 0.6399 -0.2384 *Nvtx;
    }
    else {
      corX = 0.1166 + 0.0200*Nvtx; 
      corY = 0.2764 - 0.1280*Nvtx;
    }

  }


  TVector2 corMET(0,0);
  corMET.Set( met.X() - corX, met.Y() -corY );

  // cout<<met.X()<<"  "<<corX<<"  "<<corMET.X()<<" :::: "
  //     <<met.Y()<<"  "<<corY<<"  "<<corMET.Y()<<endl;

  return corMET;

}




//
// 2Dimentional fits
// 


void ResoParametrization(string comp, bool data=false, bool Fit=false) {
  RooMsgService::instance().setGlobalKillBelow(ERROR);

  corResp=true;
  loadDB=true;

  ConfigAnalysis( data, Fit,0,0,1);

  gStyle->SetErrorX(0.);

  //Default, first met used
  string metType=metNames[0];
  recMap oMap;
  if(!loadDB) {
    oMap = GetSampleStat(comp,metType, false, false, false );

    string name="DBAscii/DB_"+metType+"_"+(data?"data":"MC")+"_"+comp;
    if(Fit) name += "_Fitted";
    name +=".dat";

    //if db exists, skip writing
    ifstream f( name.c_str(), ios::in );
    if(!f) { 
      //Create the DB
      int nBin = 23;
      int NVMax=35;
      int offset=1;
    
      ofstream file( name.c_str() , ios::out | ios::trunc );
      if(file) {
	for(int iNv=offset;iNv<NVMax;iNv++) {
	  for(int iqt=0;iqt<nBin;iqt++) {
	    if(oMap[ "centerQt" ][iNv-offset][iqt] != -1000) {
	      file << iNv<<"   "<<iqt<<"    "
		   <<oMap[ "centerQt" ][iNv-offset][iqt]<<"    "
		   <<oMap[ "u"+comp ][iNv-offset][iqt]<<"    "
		   <<oMap[ "centerQtErr" ][iNv-offset][iqt]<<"    "
		   <<oMap[ "u"+comp+"Err" ][iNv-offset][iqt]<<"   + "
		   <<0/*[iNv-offset][iqt][2]*/<<"   - "
		   <<0/*[iNv-offset][iqt][3]*/<<endl;
	    }
	  }
	}
	file.close();
      } else {
	cout<<" Error, no DB written !!!! "<<endl;
      }
      f.close();
    }
  }
  else {
    oMap = LoadDB(metType, comp, data?"data":"MC" , Fit );
    // oMap = BuildFromDB( comp, Fit );
  }

  //Do the graphs from the values ====
  //
  
  map< int , TGraphErrors* > graphs;
  int NptBin = (oMap[ "centerQt" ])[0].size();
  int NvtxBin = (oMap[ "centerQt" ]).size();
  
  for(int in=0;in<NvtxBin;in++) {

    TGraphErrors* tmp=new TGraphErrors();
    TString s="g2D_"; s+=in;
    tmp->SetName(s);
    tmp->SetTitle(s);
    int n=0;
    
    for(int ip=0;ip<NptBin;ip++) {
      float x  = (oMap[ "centerQt" ])[in][ip];
      float y = (oMap[ "u"+comp ])[in][ip];
      
      float ex = (oMap[ "centerQtErr" ])[in][ip];
      float ey = (oMap[ "u"+comp+"Err" ])[in][ip];

      if(x!=0 && x<150) { //FIXME 150
	tmp->SetPoint(n, x, y );
	tmp->SetPointError(n, ex, ey );
	n++;
      }
    }
    graphs[ in ]=tmp;
  }
  //
  //==================================
  
  GetGeneralPar(comp, graphs );

  for(int i=0;i<4;i++) {
    if(comp=="para") {
      dataFixPar[i] = FitPar[0][0][i];
      dataFixParError[i] = FitParError[0][0][i];
    }
    else {
      dataFixPar[i] = FitPar[0][1][i];
      dataFixParError[i] = FitParError[0][1][i];
    }
  }
  
  //simultaneous fit
  FitSigmaP0(metType, comp, graphs);

}


recMap LoadDB(string metname, string comp, string ds, bool FIT) {

  loadDB =true;

  int offset=1;
  int NVMax=35;
  int nBin = 23;

  recMap oMap;
  vector<vector<float> > cqt(NVMax,vector<float>(nBin,0));
  vector<vector<float> > cqte(NVMax,vector<float>(nBin,0));
  vector<vector<float> > upx(NVMax,vector<float>(nBin,0));
  vector<vector<float> > upxe(NVMax,vector<float>(nBin,0));
  
  //FIXME
  string name="DBAscii/DB_"+metname+"_"+ds+"_"+comp;
  if(FIT) name += "_Fitted";
  name +=".dat";
  cout<<name<<"   loading"<<endl;
  ifstream file(name.c_str(), ios::in);
  if(file) {

    int nV, iqT;
    string nVs;   

    while(!file.eof() ) {
     
      file >> nVs ;
      if(file.eof() ) break;
      file>> iqT;
      
      if(nVs!="#") {
	nV = atoi(nVs.c_str());
	file >> cqt[nV-offset][iqT]
	     >> upx[nV-offset][iqT]
	     >> cqte[nV-offset][iqT]
	     >> upxe[nV-offset][iqT];
	
	file >> nVs >> nVs >> nVs >> nVs;
	
      }
      else {
	file >> nVs >> nVs >> nVs >> nVs >> nVs >> nVs >> nVs >> nVs >> nVs;
      }
    }

    file.close();
  } else {
    cout<<" Error, DB not loaded !!!! "<<endl;
  }
  cout<<" End loading "<<endl;

  //filling output ==============
  oMap [ "centerQt" ]     = cqt;
  oMap [ "centerQtErr" ]  = cqte;
  oMap [ "u"+comp ]       = upx;
  oMap [ "u"+comp+"Err" ] = upxe;

  return oMap;
}



void GetGeneralPar(string comp, map<int, TGraphErrors*> graphs) {
  
  //true data first to get B parameter
  TF1* f1fit = new TF1("f1fit",fitfunc,0,200,7);
  f1fit->FixParameter(3,0); //First fit do not fit sigmaPU
  // f1fit->FixParameter(6,0); //liNFact
  // f1fit->FixParameter(5,0); //SigOOT
  f1fit->SetParNames("sqrt factor", "sqrt biais", "noise","sig_PU","NVertex","sig_PU_OOT","linear factor");

  cout<<endl<<" ========= Fit on Data ==========="<<endl<<endl;
  
  f1fit->SetParameter(0,1.56);
  f1fit->SetParameter(1,0);
  f1fit->SetParameter(2,0.5);
  f1fit->FixParameter(3,0);
  
  f1fit->FixParameter(4,0); //N vertex PU
  
  int n=10;
  
  graphs[n]->Fit("f1fit","N0","",0,200);
    
  if(comp=="para") {
 
    FitPar[0][0][0] = f1fit->GetParameter(0);
    FitPar[0][0][1] = f1fit->GetParameter(1);
    FitPar[0][0][2] = f1fit->GetParameter(2);
    FitPar[0][0][3] = f1fit->GetParameter(3);
    FitParError[0][0][0] = f1fit->GetParError(0);
    FitParError[0][0][1] = f1fit->GetParError(1);
    FitParError[0][0][2] = f1fit->GetParError(2);
    FitParError[0][0][3] = f1fit->GetParError(3);
  }
  else {
    FitPar[0][1][0] = f1fit->GetParameter(0);
    FitPar[0][1][1] = f1fit->GetParameter(1);
    FitPar[0][1][2] = f1fit->GetParameter(2);
    FitPar[0][1][3] = f1fit->GetParameter(3);
    FitParError[0][1][0] = f1fit->GetParError(0);
    FitParError[0][1][1] = f1fit->GetParError(1);
    FitParError[0][1][2] = f1fit->GetParError(2);
    FitParError[0][1][3] = f1fit->GetParError(3);
  }
  
}


vector<TObject*> FitSigmaP0(string metname, string comp, map<int, TGraphErrors*> graphs ){

  bool isMC = !dat;
  
  int NC=4;
  int NVwanted[4]={5,10,15,20};
  int mType[4] = {20,24,25,21};
  int cols[4] = {kOrange+9,38,896,814};

  //int colors2[15] = {1000,kOrange-3,kOrange+7,kOrange+10,33,kRed+1,1,kMagenta+3,3,4,5,6,7,8,11};
  UseSimultaneousFit* simFit= new UseSimultaneousFit();
  
  simFit->Parametrize(graphs,dataFixPar, dataFixParError,dat,corResp);
  c35= new TCanvas("RooFitCanvasData","RooFitCanvasData",600,600);
  c35->cd();
  
  simFit->LoadLaw();
  simFit->SetFitParameters();

  TFitterMinuit myMinimizer(4); //7
  myMinimizer.SetMinuitFCN(simFit);
  myMinimizer.CreateMinimizer();
  myMinimizer.SetPrintLevel(100);
  myMinimizer.SetErrorDef( simFit->Up() );
 
  cout<<myMinimizer.GetNumberTotalParameters()<<endl;
  cout<<" Datafixpar "<<dataFixPar[0]<<"   "<<dataFixPar[1]<<"   "<<dataFixPar[2]<<"   "<<dataFixPar[3]<<endl;
  myMinimizer.SetParameter(0,"factor",dataFixPar[0],0.0001,0,10);
  myMinimizer.SetParameter(1,"bias",dataFixPar[1],0.0001,0,20);
  myMinimizer.SetParameter(2,"sig_noise",dataFixPar[2],0.0001,0,20);
  myMinimizer.SetParameter(3,"sig_PU",3.5,0.0001,0,10);
  //cout<<myMinimizer.GetNumberTotalParameters()<<endl;
  // myMinimizer.SetParameter(4,"NPU",0,1,0,100);
  cout<<" End loading parameters "<<endl;

  // myMinimizer.FixParameter(0);
  // myMinimizer.FixParameter(1);
  // myMinimizer.FixParameter(3);


  //Minimization
  myMinimizer.Minimize(100000,0.0001);
 
  //int nk=0, nj=0;

  for(int k=0;k<4;k++) {
    //nj=0;

    for(int j=0;j<4;j++) {
     
      if(isMC) {
	matrixMC[k][j] = myMinimizer.GetCovarianceMatrixElement(k,j);
      }
      else {
	matrixData[k][j] = myMinimizer.GetCovarianceMatrixElement(k,j);
      }
    }
  }

  //==================== simul fit

  double sigPU = myMinimizer.GetParameter(3);
  double sigPUErr = myMinimizer.GetParError(3);

  //  double sigPUOOT =0;// = myMinimizer.GetParameter(5);
  // double sigPUOOTErr =0;// = myMinimizer.GetParError(5);

  for(int i=0;i<4;i++)
    {
      dataFixPar[i] = myMinimizer.GetParameter(i);
      dataFixParError[i] = myMinimizer.GetParError(i);

      if(isMC) {
	valpMC[i] =  myMinimizer.GetParameter(i);
	errpMC[i] = myMinimizer.GetParError(i);
      }
      else {
	valpData[i] =  myMinimizer.GetParameter(i);
	errpData[i] = myMinimizer.GetParError(i);
      }
    }

  //plot simul fit =======================
  cout<<endl;
  cout<<" factor --> "<<dataFixPar[0]<<" +/- "<<dataFixParError[0]<<endl;
  cout<<" bias   --> "<<dataFixPar[1]<<" +/- "<<dataFixParError[1]<<endl;
  cout<<" noise  --> "<<dataFixPar[2]<<" +/- "<<dataFixParError[2]<<endl;

  cout<<" Sigma PU "<<sigPU<<" +- "<<sigPUErr<<endl;
  // cout<<" Sigma PU OOT "<<sigPUOOT<<" +- "<<sigPUOOTErr<<endl;
  //cout<<"  size !! "<<data.size()<<endl;
  TLegend* leg=new TLegend(0.186,0.71,0.574,0.9);
  leg->SetFillColor(0);

  vector<TF1*> f;
  vector<TF1*> f2; int cnt=0;

  //out objects
  vector<TObject*> oObj;

  for(int ik=0;ik<NC;ik++) {
     
    int il = NVwanted[ik]-1;
    int ih=il;
    // if(compare && !dat)
    //   ih+=NVMax-1;

    ostringstream os;
    os << ih;
    cout<<ik<<"   "<<il<<"   "<<ih<<endl;
    TF1* f1fit2 = new TF1((os.str()).c_str(),fitfunc,0,200,7);
    //  TF1* tmpf1fit2 = new TF1(("tmp"+os.str()).c_str(),fitfunc,0,200,7);
    for(int i=0;i<4;i++)
      {
	f1fit2->FixParameter(i,dataFixPar[i]);
	//	tmpf1fit2->FixParameter(i,SavePar[0][i]);
      }
    //tmpf1fit2->FixParameter(4,SavePar[0][4]);
    //tmpf1fit2->FixParameter(5,SavePar[0][5]);
    f.push_back(f1fit2);
    //  f2.push_back(tmpf1fit2);
 
    string nf = "fMC";
    if(dat) nf = "fdata";

    f[ik]->SetName( (nf+os.str()).c_str() ); 
    f[ik]->FixParameter(4,il);
    f[ik]->FixParameter(3,sigPU);
    f[ik]->SetParError(3,sigPUErr);
    f[ik]->FixParameter(5,0);
    f[ik]->SetLineColor(cols[ik]);
    f[ik]->SetLineWidth(2);

    if(il==0) cout<<" evaluation at 0 : "<<f[ik]->Eval(0.)<<endl;

    graphs[ih]->SetLineWidth(2);
    graphs[ih]->SetMarkerColor(cols[ik]); //ih
    graphs[ih]->SetMarkerStyle(mType[ik]);
    graphs[ih]->GetYaxis()->SetRangeUser(0,30);
    graphs[ih]->GetXaxis()->SetRangeUser(0,300);
    // if(il==3)
    //   graphs[ih]->SetMarkerSize(1.2);
    
    if(ik==0) {
      graphs[ih]->Draw("AP");
      string s="#sigma(u_{#perp}  )  [GeV]";
      if(comp=="para") s="#sigma(u_{||})  [GeV]";

      graphs[ih]->GetYaxis()->SetTitle( s.c_str() );
      graphs[ih]->GetXaxis()->SetTitle("Z q_{T}  [GeV]");
    }
    else
      graphs[ih]->Draw("Psame");

    //if( !(il==5 || il==10 || il==15 || il==20) ) continue;
    //cout<<ih<<endl;
    
    f[ik]->DrawCopy("Lsame+");

    oObj.push_back( graphs[ih] ); 
    oObj.push_back( f[ik] );
    
    string legt = "1 vertex (no PU)";
    if(il!=0)
      {
	ostringstream os2, os3;
	os2 << il;
	os3 << il+1;
	legt = os3.str() + " vertices"; 
      }

    leg->AddEntry(graphs[ih],legt.c_str(),"pl");
    cnt++;
  }
  leg->Draw("same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.042); //cout<<sigPUErr<<" bordel de merde "<<endl;

  if(comp=="perp") {
    latex.DrawLatex(0.50,0.213,Form("#sigma_{PU} = %.2f #pm %.2f GeV",sigPU,sigPUErr) );
    // latex.DrawLatex(0.57,0.702,Form("#sigma_{PU OOT} = %.2f #pm %.2f GeV",sigPUOOT,sigPUOOTErr) );
  }
  else {
    latex.DrawLatex(0.50,0.213,Form("#sigma_{PU} = %.2f #pm %.2f GeV",sigPU,sigPUErr) ); 
    // latex.DrawLatex(0.57,0.302,Form("#sigma_{PU OOT} = %.2f #pm %.2f GeV",sigPUOOT,sigPUOOTErr) );
  }
  string name = metname+"_"+comp;

  //  cmsPrel(1500);
  string type="MC";
  if(dat)
    type="Data";    
  
  if(fit)
    type+="_Fitted";
  else
    type+="_RMS";

  if(corResp) type+="_corResp";
  
  c35->SaveAs( ("Plots/"+name+"_"+type+".png").c_str() );
  c35->SaveAs( ("Plots/"+name+"_"+type+".pdf").c_str() );
  c35->SaveAs( ("Plots/"+name+"_"+type+".eps").c_str() );
  c35->SaveAs( ("Plots/"+name+"_"+type+".C").c_str() );
  c35->SaveAs( ("Plots/"+name+"_"+type+".root").c_str() );
  
  return oObj;

}
 
double fitfunc(double *x, double *par) {

  double response =  1./respfunc(x[0]); 
  if(!corResp) response=1.;
  double nVtx = par[4]+1;
  int nIT = (nVtx-1)/0.7;
  
  double f = sqrt( pow( sqrt(x[0])*par[0] + par[1] ,2 )   +par[2]*par[2]*pow(response, 2 ) + (par[3]*par[3]*nIT) *pow(response, 2 ) );
  //double f = sqrt( (x[0])*par[0] + par[1]  +par[2]*par[2]*pow(response, 2 ) + (par[3]*par[3]*nIT) *pow(response, 2 ));
  return f;

}

double respfunc(double pt) {
  
  float r =1;
  if(dat) {
    r = 0.993 + 4.603/(126.2+pt) - 36.47/(111.70 +pt*pt);
    //r = 1.00678+ 1.17817/(15.7148 +pt) - 38.5088/(108.963+pt*pt);
 }
  else {
    r = 1.00 + 3.81/(59.42+pt) - 26.80/(84.71+pt*pt);
    // r = 1.00578+ 1.92346/(17.2487+pt)  -33.5747/(91.5004+pt*pt);
  }
  return r;
}




void MakeAllPlots() {

  string comp[2]={"para","perp"};
  for(int ic=0;ic<2;ic++) {
    for(int id=0;id<2;id++) {
      for(int ik=0;ik<2;ik++) {
      
      ResoParametrization(comp[ic], id, ik );
	
	}
  }
}


}



void LoadVtxReweight() {

  TH1F *ratioNVertex__1 = new TH1F("ratioNVertex__1","92",50,0,50);
  ratioNVertex__1->SetBinContent(1,4.352382);
   ratioNVertex__1->SetBinContent(2,0.912782);
   ratioNVertex__1->SetBinContent(3,0.9450765);
   ratioNVertex__1->SetBinContent(4,0.9245538);
   ratioNVertex__1->SetBinContent(5,0.9409372);
   ratioNVertex__1->SetBinContent(6,0.9554647);
   ratioNVertex__1->SetBinContent(7,0.9564364);
   ratioNVertex__1->SetBinContent(8,0.962751);
   ratioNVertex__1->SetBinContent(9,0.9730481);
   ratioNVertex__1->SetBinContent(10,0.9775884);
   ratioNVertex__1->SetBinContent(11,0.987715);
   ratioNVertex__1->SetBinContent(12,0.9920352);
   ratioNVertex__1->SetBinContent(13,1.002308);
   ratioNVertex__1->SetBinContent(14,1.001181);
   ratioNVertex__1->SetBinContent(15,1.006853);
   ratioNVertex__1->SetBinContent(16,1.009103);
   ratioNVertex__1->SetBinContent(17,1.014458);
   ratioNVertex__1->SetBinContent(18,1.01718);
   ratioNVertex__1->SetBinContent(19,1.020113);
   ratioNVertex__1->SetBinContent(20,1.02148);
   ratioNVertex__1->SetBinContent(21,1.029908);
   ratioNVertex__1->SetBinContent(22,1.027123);
   ratioNVertex__1->SetBinContent(23,1.028488);
   ratioNVertex__1->SetBinContent(24,1.038353);
   ratioNVertex__1->SetBinContent(25,1.028244);
   ratioNVertex__1->SetBinContent(26,1.032906);
   ratioNVertex__1->SetBinContent(27,1.041499);
   ratioNVertex__1->SetBinContent(28,1.061431);
   ratioNVertex__1->SetBinContent(29,1.049645);
   ratioNVertex__1->SetBinContent(30,1.027846);
   ratioNVertex__1->SetBinContent(31,1.042164);
   ratioNVertex__1->SetBinContent(32,1.079537);
   ratioNVertex__1->SetBinContent(33,1.151077);
   ratioNVertex__1->SetBinContent(34,1.081402);
   ratioNVertex__1->SetBinContent(35,1.116077);
   ratioNVertex__1->SetBinContent(36,1.025585);
   ratioNVertex__1->SetBinContent(37,1.196205);
   ratioNVertex__1->SetBinContent(38,1.240716);
   ratioNVertex__1->SetBinContent(39,0.8237543);
   ratioNVertex__1->SetBinContent(40,0.9180249);
   ratioNVertex__1->SetBinContent(41,0.7386603);
   ratioNVertex__1->SetBinContent(42,0.4086069);
   ratioNVertex__1->SetBinContent(43,0.3347449);
   ratioNVertex__1->SetBinContent(44,0.3824325);
   ratioNVertex__1->SetBinContent(45,0.2143044);
   ratioNVertex__1->SetBinContent(48,13.12211);


   vtxW = (TH1F*)ratioNVertex__1->Clone();
   vtxW->Reset("ICEM");
   for(int i=0;i<ratioNVertex__1->GetNbinsX()+2;i++) {
     vtxW->SetBinContent(i, ratioNVertex__1->GetBinContent(i));
    //  if(ratioNVertex__1->GetBinContent(i)!=0)
   //     vtxW->SetBinContent(i, 1./ratioNVertex__1->GetBinContent(i));
   //   else
   //     vtxW->SetBinContent(i, 0.);
   }
}

void LoadQtReweight() {

  qTW = new TH1F("qTW","qTW",50,0,400);
   qTW->SetBinContent(13,1.631756);
   qTW->SetBinContent(14,1.080202);
   qTW->SetBinContent(15,1.007457);
   qTW->SetBinContent(16,0.9652824);
   qTW->SetBinContent(17,0.7878233);
   qTW->SetBinContent(18,0.712092);
   qTW->SetBinContent(19,0.7016643);
   qTW->SetBinContent(20,0.7437111);
   qTW->SetBinContent(21,0.6738006);
   qTW->SetBinContent(22,0.7243134);
   qTW->SetBinContent(23,0.6183339);
   qTW->SetBinContent(24,0.5528829);
   qTW->SetBinContent(25,0.5928317);
   qTW->SetBinContent(26,0.6035337);
   qTW->SetBinContent(27,0.5171213);
   qTW->SetBinContent(28,0.5957144);
   qTW->SetBinContent(29,0.5552399);
   qTW->SetBinContent(30,0.6545408);
   qTW->SetBinContent(31,0.4946534);
   qTW->SetBinContent(32,0.5904603);
   qTW->SetBinContent(33,0.7835543);
   qTW->SetBinContent(34,0.4619016);
   qTW->SetBinContent(35,0.4381354);
   qTW->SetBinContent(36,0.5330759);
   qTW->SetBinContent(37,0.4976588);
   qTW->SetBinContent(38,0.4986609);
   qTW->SetBinContent(39,0.4439213);
   qTW->SetBinContent(40,0.5714276);
   qTW->SetBinContent(41,0.5304265);
   qTW->SetBinContent(42,0.4869582);
   qTW->SetBinContent(43,0.427651);
   qTW->SetBinContent(44,0.4483841);
   qTW->SetBinContent(45,0.4146701);
   qTW->SetBinContent(46,0.4631154);
   qTW->SetBinContent(47,0.6068047);
   qTW->SetBinContent(48,0.5008188);
   qTW->SetBinContent(49,0.5521687);
   qTW->SetBinContent(50,0.3970008);
   qTW->SetEntries(38);
   qTW->SetLineStyle(0);
   qTW->SetMarkerStyle(20);
   qTW->GetXaxis()->SetLabelFont(42);
   qTW->GetXaxis()->SetLabelOffset(0.007);
   qTW->GetXaxis()->SetLabelSize(0.05);
   qTW->GetXaxis()->SetTitleSize(0.06);
   qTW->GetXaxis()->SetTitleOffset(0.9);
   qTW->GetXaxis()->SetTitleFont(42);
   qTW->GetYaxis()->SetLabelFont(42);
   qTW->GetYaxis()->SetLabelOffset(0.007);
   qTW->GetYaxis()->SetLabelSize(0.05);
   qTW->GetYaxis()->SetTitleSize(0.06);
   qTW->GetYaxis()->SetTitleOffset(1.1);
   qTW->GetYaxis()->SetTitleFont(42);
   qTW->GetZaxis()->SetLabelFont(42);
   qTW->GetZaxis()->SetLabelOffset(0.007);
   qTW->GetZaxis()->SetLabelSize(0.05);
   qTW->GetZaxis()->SetTitleSize(0.06);
   qTW->GetZaxis()->SetTitleFont(42);


}

float GetQtW(float qT) {
  int bin = qTW->GetXaxis()->FindBin( qT );
  
  if(bin>50) return 1.;
  else return qTW->GetBinContent(bin);
}


float GetVtwW(int nvtx) {
  //if(nvtx==44)cout<<vtxW->GetBinContent( vtxW->GetXaxis()->FindBin( (float)nvtx)+1 )<<endl;
  //cout<<vtxW->GetBinContent( vtxW->GetXaxis()->FindBin( (float)nvtx) )<<endl;
  //  cout<<vtxW<<endl;
  if(nvtx>= vtxW->GetNbinsX()) return 0;//nvtx=vtxW->GetNbinsX();
  return vtxW->GetBinContent( vtxW->GetXaxis()->FindBin( (float)nvtx)+1 );
}


vector<TPad*> PreparePads(TCanvas *& c2) {
  
  //Prepare pads
  double t = 35;
  double h = 535;
  double he = 5;
  double b = 70;
  double hl = 75;
  double H = t + h + b + hl + he;

  double w = 535;
  double g = 90;
  double d = 25;
  double W = w + g + d;

  c2=new TCanvas("c2","c2",W,H);
  TPad* pHigh = new TPad( "phigh", "phigh", 
			  0, (hl+b+he)/H , W/W, 1. );
  pHigh->SetLeftMargin(  g/W );
  pHigh->SetRightMargin( d/W );
  pHigh->SetTopMargin(  t/H );
  pHigh->SetBottomMargin( he/H );
  c2->cd();
  TPad* pLow = new TPad( "plow", "plow", 
			 0, 0 , W/W, (hl+b+he)/H );
  pLow->SetLeftMargin(  g/W );
  pLow->SetRightMargin( d/W );
  pLow->SetTopMargin(  he/H );
  pLow->SetBottomMargin( b/(hl+b+he) );
  
  c2->cd();
  pHigh->Draw();
  c2->cd();
  pLow->Draw();
  c2->cd();

  vector<TPad*> v;
  v.push_back(pHigh);
  v.push_back(pLow);
  
  return v;
}

vector<TGraphAsymmErrors*> getRatios(vector<TGraphAsymmErrors*> graphs) {

  vector<TGraphAsymmErrors*> ratios;
  
  for(size_t i=0;i<graphs.size();i+=2) {
    TGraphAsymmErrors* g=new TGraphAsymmErrors();
    double xd,yd,eyd;
    double xm,ym,eymh, eyml;
    int n=0;
    for(int j=0;j<graphs[i]->GetN();j++) {
      graphs[i+1]->GetPoint(j,xd,yd);
      eyd = graphs[i+1]->GetErrorYhigh(j);

      graphs[i]->GetPoint(j,xm,ym);
      eymh = graphs[i]->GetErrorYhigh(j);
      eyml = graphs[i]->GetErrorYlow(j);

      if(ym!=0) {
	g->SetPoint(n, xd,yd/ym );

	cout<<yd<<" +/- "<<eyd<<" ; "<<ym<<" +/- "<<eyml<<" --> "<<yd/ym<<" ; "<<sqrt( pow(eyd/yd ,2) + pow(eyml/ym ,2) )<<"  "<<(yd/ym)*sqrt( pow(eyd/yd ,2) + pow(eyml/ym ,2) )<<endl;

	g->SetPointError(n, 0,0, (yd/ym)*sqrt( pow(eyd/yd ,2) + pow(eyml/ym ,2) ),(yd/ym)*sqrt( pow(eyd/yd ,2) + pow(eyml/ym ,2) ));
	n++;
      }

    } //points

    g->SetLineColor( graphs[i+1]->GetLineColor() );
    g->SetMarkerColor( graphs[i+1]->GetMarkerColor() );
    g->SetMarkerStyle( graphs[i+1]->GetMarkerStyle() );
    ratios.push_back(g);
  }//graphs

  return ratios;
}




TPolyLine* GetSystBand(string comp, float xmax, int NV, int nqt, bool loaddb) {

  // if(loaddb)
  //   LoadDBUnc(comp,"",1);
  
  if(comp=="ff" || loaddb) cout<<"gloubi"<<endl;

  int nbin = 23;// centerData[0].size()-1;
  
  double BVar[24] = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,230,260,300};
  //double BVar[15] = {0,5,10,20,30,40,60,80,100,120,140,170,200,250,300};

  vector<double> BinVar;

  int vtx=0, qt=0;
  float factor=1.;
  if(nqt == -1) {
    nbin = 23;
    vtx = NV;
    factor=1.;
    for(int i=0;i<24;i++)
      BinVar.push_back(BVar[i] );
  }
  if(NV==-1) {
    nbin = 34;
    qt = nqt;
    factor = sqrt(2);
    for(int i=0;i<35;i++)
      BinVar.push_back( (double)i );
  }
  
  TPolyLine* band = new TPolyLine(4*nbin+1); 

  vector< float > vecx_(4*nbin+2);
  vector< float > vecy_(4*nbin+2); float lsyst;



  for( int ibin=0; ibin<nbin; ibin++ )
    {
      if(nqt == -1) {
	qt = ibin;
      }
      if( NV==-1) {
	vtx =ibin;
      }   

      //  cout<<vtx<<"   "<<qt<<"   "<<(BinVar[ibin]+BinVar[ibin+1])/2.<<"   "<<Systematics[vtx][qt]<<endl; 
      if(xmax > 1.0001*(BinVar[ibin]+BinVar[ibin+1])/2. && ibin<18) {
	vecx_[2*ibin]   = 1.0001*(BinVar[ibin]+BinVar[ibin+1])/2.;
	vecy_[2*ibin]   = Systematics[vtx][qt]/factor;
	lsyst = Systematics[vtx][qt]/factor;
      }
	
      else {
	if(NV==-1) {
	  vecx_[2*ibin]   = xmax;
	  vecy_[2*ibin]   = lsyst;
	}
	else {
	  vecx_[2*ibin]   = 1.0001*(BinVar[ibin]+BinVar[ibin+1])/2.;
	  vecy_[2*ibin]   = lsyst;
	  //  cout<<" bili "<<vtx<<"   "<<qt<<"   "<<(BinVar[ibin]+BinVar[ibin+1])/2.<<"   "<<Systematics[vtx][qt]<<endl; 
	}
      }
	
      if(xmax > 0.9999*(BinVar[ibin+1]+BinVar[ibin+2])/2. && ibin<18) {


	if(nqt == -1) {
	  qt = ibin+1;
	}
	if( NV==-1) {
	  vtx =ibin+1;
	}   

	vecx_[2*ibin+1] = 0.9999*(BinVar[ibin+1]+BinVar[ibin+2])/2.;
	vecy_[2*ibin+1] = Systematics[vtx][qt]/factor;
	lsyst = Systematics[vtx][qt]/factor;
      } else {
	if(NV==-1) {
	  vecx_[2*ibin+1] = xmax;
	  vecy_[2*ibin+1] = lsyst;
	}
	else {
	  vecx_[2*ibin+1]   = 0.9999*(BinVar[ibin+1]+BinVar[ibin+2])/2.;
	  vecy_[2*ibin+1]   = lsyst;
	}
      }

	
      vecx_[4*nbin-2-2*ibin]   = vecx_[2*ibin+1];
      vecx_[4*nbin-2-2*ibin+1] = vecx_[2*ibin];
	
      vecy_[4*nbin-2-2*ibin]   = -vecy_[2*ibin+1];
      vecy_[4*nbin-2-2*ibin+1] = -vecy_[2*ibin];

    }
      //vecx_[0] = 0;
      //vecy_[0] = vecy_[4*nbin+1];
      //Border
  // vecx_[0] = 0;vecx_[1] = 0;
  vecx_[4*nbin] = 0; vecy_[4*nbin] = vecy_[4*nbin-1];
  vecx_[4*nbin+1] = 0; 
  vecy_[4*nbin+1] = -vecy_[4*nbin-1];
  // vecy_[4*nbin+2] = vecy_[4*nbin-1];

  for( size_t ipt=0; ipt<vecx_.size(); ipt++ )
    {
      // cout<<ipt<<"   "<<vecx_[ipt]<<"   "<<vecy_[ipt]<<"   "<<xmax<<endl;
      band->SetPoint(ipt, vecx_[ipt], 1+vecy_[ipt] );

    }
  band->SetLineColor(kGray+1);
  band->SetFillColor(kGray+1);
  band->SetFillStyle(3002);

  return band;

}



float ComputeUncertainties( map<string, TH1*> uncs ) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  map<string, float> vals;

  map<string, TH1*>::iterator itVar; 
  for(itVar=uncs.begin();itVar!=uncs.end();itVar++) {
    
    if(fit) {
      RooRealVar x("x","x",0,-150,150);
      RooDataHist Hist("Hist","Hist",x, (TH1*)itVar->second->Clone() );
      RooRealVar g_w("g_w","width Gaus", 10., 0., 40., "GeV");
      RooRealVar gamma_Z0( "gamma_Z0_U", "Z0 width",2.3, 0, 20, "GeV" );
      RooRealVar v_m("v_m","v_m",0, -40 , 40 , "GeV" );
      RooVoigtian* voigt=new RooVoigtian("voigt","Voightian",x, v_m, gamma_Z0, g_w);
      
      RooFitResult* result = voigt->fitTo( (Hist) ,RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );
      
      double sigma = g_w.getVal();
      double gamma = gamma_Z0.getVal();
      double fwhm = FWHM( sigma, gamma );
      double val=fwhm/2.3548;
      vals[ itVar->first ] = val;

      cout<<val<<"   "<<itVar->first<<"   "<<v_m.getVal()<<endl;
    }
    else {
      cout<<itVar->second->GetName()<<"   "<<itVar->first<<" --> "<<itVar->second->GetMean()<<endl;
      vals[ itVar->first ] = _resp?itVar->second->GetMean():itVar->second->GetRMS();
    }
  } //all uncs

  float central = vals["central"];
  float totUnc=0;
  //cout<<" central : "<<central<<endl;
  
  for(map<string, float>::iterator it=vals.begin();
      it!=vals.end();it++) {

    if(it->first=="central") continue;

    if( (it->first.find("Up")==(size_t)-1 &&  
	 it->first.find("Do")==(size_t)-1 ) 
	|| it->first.find("JetResDown") //fixme, this is ugly
	) {
      totUnc += pow(it->second - central, 2);
    }
    else { //symztrisation
      totUnc += pow(it->second - central, 2)/4.;
    } 
  }
  
  return sqrt(totUnc);
  
}



int main() {

  CompaRespCurves(0);

}




void DoSingleFit(bool isdata, string comp, string metType, int qtbin) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  
  apf=true;
  bool AllVtx=true;
  bool AllPts=false;
  bool resp=true;
  // bool fit=true;
  // bool dat=false;
  
  ConfigAnalysis(isdata, 1, 0, 1, 1);

  if(resp) fit=false;

  metNames.clear();
  metNames.push_back(metType);

  _resp=resp;
  bool debug_=false;
  _sumEt = AllVtx && AllPts;
  if(!loaded)
    Data = LoadTree();

  cout<<" Getting sample :: component = "<<comp<<" for "<<metType
      <<"  allvtx="<<AllVtx<<" allpts="<<AllPts<<"  "<<dat<<endl;

  cout<<" number of entries "<<Data->GetEntries()<<endl;

  TTree* ZData = (TTree*)Data->Clone();
  cout<<" Tree cloned "<<endl;
  //int nBin = 13;
  //double BinVar[17] = {0,5,10,20,30,40,50,60,70,80,90,100,110,120,140,160,200};
  //double BinVar[14] = {0.,10.,20.,30.,40.,60.,80.,100.,120.,140.,170.,200.,250.,300.};
  //double BinVar[35] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,230,260,300};
  // double BinVar[24] = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,230,260,300};
  
  int nBin=18;
  double bVar[19] = {0,5,10,15,20,25,30,40,60,85,100,120,145,165,200,250,300,400,600};
  double *BinVar(0);
  //  int nBin=19;
  //double BinVar[20] = {0,5,10,15,20,25,30,40,60,85,100,120,145,165,200,250,300,350,400,600};

  //int mType[16] = {100,20,24,25,26,27,21,2,2,2,2,2,2,2,2,2};
  //int colors[16] = {1000,kOrange+9,1,38,896,814,kMagenta+3,2,3,4,5,6,7,8,9,11};
  
  // nBin=1;
  // double BinVar[2] = {10,20};

  int nVtxMax=35;

  int offset=1;
  int NVMax = nVtxMax;
  
  //bool sumEt=false;

  if(AllVtx) {
    offset=0;
    NVMax=1;
  }

  if(AllPts) {
    nBin=1;
  }

  if(AllVtx && AllPts) {
    offset=0;
    NVMax=1;
    nBin=23;
    //sumEt=true;
    BinVar=new double[24];
    for(int i=0;i<24;i++)
      BinVar[i]= 100*i;
  }
  else {
    BinVar=new double[19];
    for(int i=0;i<19;i++) BinVar[i]=bVar[i];
  }

  //prepare output ====================================
  map<string, vector<vector<float> > > oMap;
  vector<vector<float> > cqt(NVMax,vector<float>(nBin,0));
  vector<vector<float> > cqte(NVMax,vector<float>(nBin,0));
  vector<vector<float> > upx(NVMax,vector<float>(nBin,0));
  vector<vector<float> > upxe(NVMax,vector<float>(nBin,0));
  vector<vector<float> > upxu(NVMax,vector<float>(nBin,0));

  // ====================================================
  map<string, TH1* > uncHistosU;
 
  int Nvtx=qtbin;
  //for(int Nvtx=offset;Nvtx<NVMax;Nvtx++) {
  RooRealVar weight("weight","weight",-10,100);
  RooRealVar Vtxcat("nVtx","nVtx", AllVtx?-0.1:(Nvtx - 0.1) , AllVtx?100.:(Nvtx+0.1) ); //NVertex have to be changed 
  //with the name of the
  //corresponding tree branch

  for(int iqt=qtbin;iqt<qtbin+1;iqt++) {
    cout<<" Start bin "<<iqt<<"   "<<BinVar[iqt]<<"    "<<BinVar[iqt+1]<<"   "<<Nvtx<<"  "<<AllVtx<<"   with sumET?"<<_sumEt<<"   "<<comp<<endl;

    string n=_sumEt?((metType+"_sumEt").c_str()):"Zpt"; //cout<<n<<"  "<<AllVtx<<"  "<<AllPts<<"  "<<sumEt<<endl;
    RooRealVar qTcat(n.c_str(),n.c_str(),(AllPts && !_sumEt)?0.:BinVar[iqt],(AllPts && !_sumEt)?1000.:BinVar[iqt+1]);
    //RooRealVar qTcat("Zpt","Zpt",(AllPts)?0.:BinVar[iqt],(AllPts)?1000.:BinVar[iqt+1]);
    double m,M,c;
    if(comp=="para")
      {
	c = (BinVar[iqt]+BinVar[iqt+1])/2.;
	m = AllPts? -400 : -400+(-10*iqt)+5; // -c-100;
	M = AllPts? 400  : (-10*iqt)+400+5; //-c+100
      }
    else {
      //	c = 0;
      m = -100; //FIXME 100
      M = 100;
    }
    //cout<<iqt<<"   "<<BinVar[iqt]<<"    "<<BinVar[iqt+1]<<"   "<<m<<"    "<<M<<endl;
    RooRealVar upara((metType+"_u"+comp).c_str(),"upara",m,M);
    RooRealVar resoP( (metType+"_redupara").c_str(),"resoPar", -150,150 ); //100

    //by pass for debugging
    if(debug_) {
      cqt[Nvtx-offset][iqt] = (BinVar[iqt]+BinVar[iqt+1])/2.;
      cqte[Nvtx-offset][iqt] = 0;

      upx[Nvtx-offset][iqt] = Nvtx;
      upxe[Nvtx-offset][iqt] = 0;
      continue;
    }

    RooDataSet* Dpara;
    if(dat)
      { 
	if(!_sumEt)
	  Dpara = new RooDataSet("data","Data",ZData, RooArgSet(Vtxcat,qTcat,upara,resoP) );
	else
	  Dpara = new RooDataSet("data","Data",ZData, RooArgSet(Vtxcat,qTcat,upara,weight,resoP),"","weight" );
      }
    else
      Dpara = new RooDataSet("data","Data",ZData, RooArgSet(Vtxcat,qTcat,upara,weight,resoP),"","weight" );
	
    if(computeUnc && !dat && !noPU) {
      for(size_t iu=0;iu<metUncNames.size();iu++) {
	delete uncHistosU[ metUncNames[iu] ];
	delete uncHistosU[ "central" ];
	RooRealVar uncU((metUncNames[iu]+"_"+metType+"_u"+comp).c_str(),"upara",m,M);
	RooRealVar uncRU( (metUncNames[iu]+"_"+metType+"_redupara").c_str(),"resoPar", -150,150 );
	RooDataSet* uncDS = new RooDataSet("data","Data",ZData, RooArgSet(Vtxcat,qTcat,uncU,weight,uncRU),"","weight" );
	cout<<uncDS->sumEntries()<<"   "<<uncDS->numEntries()<< "---------- "<<uncDS->meanVar(uncU)->getVal()<<endl;
	if(comp=="para" && !resp) {
	  uncHistosU[ metUncNames[iu] ] = uncDS->createHistogram( (metUncNames[iu]+"_"+metType+"_redupara").c_str(), 150 );
	  uncHistosU[ "central" ] = Dpara->createHistogram( (metType+"_redupara").c_str(), 150 );
	}
	else if(resp) {
	  uncHistosU[ metUncNames[iu] ] = uncDS->createHistogram( (metUncNames[iu]+"_"+metType+"_upara").c_str(), 150 );
	  uncHistosU[ "central" ] = Dpara->createHistogram( (metType+"_upara").c_str(), 150 );
	}
	else {
	  uncHistosU[ metUncNames[iu] ] = uncDS->createHistogram( (metUncNames[iu]+"_"+metType+"_u"+comp).c_str(), 150 );
	  uncHistosU[ "central" ] = Dpara->createHistogram( (metType+"_u"+comp).c_str(), 150 );
	}
	delete uncDS;
      }//unc loop
    }

    //cout<<" gloubi "<<endl;
    RooRealVar* meanqT = (*Dpara).meanVar(qTcat);
    RooRealVar* rms = (*Dpara).rmsVar(upara);
    RooRealVar* mean = (*Dpara).meanVar(upara);
    if(comp=="para" && !resp)
      rms = (*Dpara).rmsVar(resoP);

	
    double fwhm=-1, efwhm=-1;
    double peak=-1, epeak=-1;
    cout<<" rms over, now fit or store "<<endl;
    if(fit) {
      double qt= 0;//-1* ( BinVar[iqt+1] + BinVar[iqt] )/2.;
      double qtm = -30;//qt -100;
      double qtM = 30;//qt +100;
      // double qt= -1* ( BinVar[iqt+1] + BinVar[iqt] )/2.;
      // double qtm = qt -100;
      // double qtM = qt +100;

      if(comp=="perp") {
	qt =0;
	qtm = -20; //-100
	qtM = 20; //100
      }

      // if(Dpara->sumEntries() < 0.001 || Dpara->numEntries() ==0 ) 
      //   cout<< " No Entry for bin ["<<m<<":"<<M<<"] / ["<<(AllVtx?0:Nvtx)<<":"<<(AllVtx?100.:Nvtx)<<"]"<<endl;
      
      RooRealVar g_w("g_w","width Gaus", 12., 0., 40., "GeV");
      RooRealVar gamma_Z0( "gamma_Z0", "Z0 width",3.2, 0, 10, "GeV" );
      RooRealVar v_m("v_m","v_m",7, qtm , qtM , "GeV" );
      RooVoigtian* voigt= new RooVoigtian("voigt","Voightian",upara, v_m, gamma_Z0, g_w);
      if(comp=="para")
	voigt= new RooVoigtian("voigt","Voightian",resoP, v_m, gamma_Z0, g_w);

      RooFitResult* result = voigt->fitTo( (*Dpara) ,RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE));//, RooFit::PrintLevel(-1) 
      //  cout<<qt<<"    "<<Nvtx<<" ==================> The status "<<result->status()<<endl;
	  
      ostringstream os;
      os<<iqt;
      TCanvas* c1 =new TCanvas("c1","ZPeak",600,600);
      RooPlot* frame = upara.frame() ;
      if(comp=="para")
	frame = resoP.frame() ;
      Dpara->plotOn(frame) ;
      voigt->plotOn(frame) ;
      voigt->paramOn(frame);
      frame->Draw();
      c1->SaveAs( ("projs/singlefit.root") );

      //Get the FWHM
      double sigma = g_w.getVal();
      double gamma = gamma_Z0.getVal();
      double esigma = g_w.getError();
      double egamma = gamma_Z0.getError();

      double Vsg = result->correlation(g_w,gamma_Z0) ;
      double Vgs = result->correlation(gamma_Z0,g_w) ;
      double Vss = result->correlation(g_w,g_w) ;
      double Vgg = result->correlation(gamma_Z0,gamma_Z0) ;

      fwhm = FWHM( sigma, gamma );
      efwhm = FWHMError(/*sigma, gamma,*/ esigma, egamma,
			Vss, Vsg, Vgs, Vgg);

      peak = v_m.getVal();
      epeak = v_m.getError();
      //delete voigt;
    }
	
    if( (*Dpara).sumEntries()>20) {
	
      // cqt[Nvtx-offset][iqt] = meanqT->getVal();
      // cqte[Nvtx-offset][iqt] = meanqT->getError();

      if(!fit) {
	// if(dat)
	//   upx[Nvtx-offset][iqt] = resp?mean->getVal():rms->getVal();// /(0.985 + 0.000387* meanqT->getVal() ); //FIXME
	// else
	//   upx[Nvtx-offset][iqt] = resp?mean->getVal():rms->getVal();
	// upxe[Nvtx-offset][iqt] = resp?mean->getError():rms->getError ();

	cout<<" store over "<<mean->getVal()<<"  "<<mean->getError()<<endl;
	if(isnan(rms->getVal())) {
	  // TH1* ht(0);
	  // ht = Dpara->createHistogram( (metType+"_u"+comp).c_str(),150);
	  // if(comp=="para" && !resp)
	  //   ht = Dpara->createHistogram( (metType+"_redupara").c_str(),150);
	  // TH1* ht2 = Dpara->createHistogram("Zpt",100);
	  // TFile* tr=new TFile("fr.root","RECREATE");
	  // ht->Write();
	  // ht2->Write();
	  // tr->Close();
	  //cout<<rms->getVal()<<"   "<<rms->getError()<<"   "<<mean->getVal()<<"  "<<ht->GetRMS()<<"   "<<ht2->GetMean()<<endl;

	  // cqt[Nvtx-offset][iqt] = ht2->GetMean();
	  // cqte[Nvtx-offset][iqt] =ht2->GetMeanError();

	  // upx[Nvtx-offset][iqt] =resp?ht->GetMean():ht->GetRMS();
	  // upxe[Nvtx-offset][iqt] =resp?ht->GetMeanError():ht->GetRMSError();

	}

	cout<<" ---> "<<mean->getVal()<<"  "<<mean->getError()<<"  "<<rms->getVal()<<"  "<<Dpara->numEntries()<<"  "<<(*Dpara).sumEntries()<<"   "<<meanqT->getVal()<<"   "<<fwhm/2.3548<<"   "<<efwhm/2.3548<<endl;
      }
      else {
	// upx[Nvtx-offset][iqt] = resp?peak:fwhm/2.3548;
	// upxe[Nvtx-offset][iqt] = resp?epeak:efwhm/2.3548;
	cout<<fwhm/2.3548<<"   "<<efwhm/2.3548<<endl;
	// TH1* ht(0);
	// ht = Dpara->createHistogram( (metType+"_u"+comp).c_str(),100);
	// if(comp=="para" && !resp)
	//   ht = Dpara->createHistogram( (metType+"_redupara").c_str(),100);
	// TH1* ht2 = Dpara->createHistogram("Zpt",100);
	// TFile* tr=new TFile("tmp.root","RECREATE");
	// ht->Write();
	// ht2->Write();
	// tr->Close();
      }
	  
      if(computeUnc && !dat && !noPU) {
	float unc = ComputeUncertainties(uncHistosU);
	//	upxu[Nvtx-offset][iqt] = resp?unc:unc;
	cout<<" Uncertainty : "<<unc<<endl;
	//upx[Nvtx-offset][iqt] = upx[Nvtx-offset][iqt]-unc;
      }	  
	  

    }// number of event > 20
    //cout<<Dpara<<endl;
    //delete Dpara;
  }//qt
  
}



void ManualFit() {

  TFile* file=new TFile("tmp.root","READ");
  TH1F* h =(TH1F*)file->Get("data__pat_patPFMetMVAPhi_redupara");
  cout<<h<<endl;
  RooRealVar x("x","x",0,-100,100);
  RooDataHist Hist("Hist","Hist",x, h );
  RooRealVar g_w("g_w","width Gaus", 10., 0., 40., "GeV");
  RooRealVar gamma_Z0( "gamma_Z0_U", "Z0 width",2.3, 0, 20, "GeV" );
  RooRealVar v_m("v_m","v_m",0, -10 , 10 , "GeV" );
  RooVoigtian* voigt=new RooVoigtian("voigt","Voightian",x, v_m, gamma_Z0, g_w);
  
  RooFitResult* result = voigt->fitTo( (Hist) ,RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE) );//, RooFit::PrintLevel(-1)
      
  TCanvas* c1 =new TCanvas("c1","ZPeak",600,600);
  RooPlot* frame = x.frame() ;
  
  Hist.plotOn(frame) ;
  voigt->plotOn(frame) ;
  frame->Draw();
  
  double sigma = g_w.getVal();
  double gamma = gamma_Z0.getVal();
  double esigma = g_w.getError();
  double egamma = gamma_Z0.getError();

  double Vsg = result->correlation(g_w,gamma_Z0) ;
  double Vgs = result->correlation(gamma_Z0,g_w) ;
  double Vss = result->correlation(g_w,g_w) ;
  double Vgg = result->correlation(gamma_Z0,gamma_Z0) ;
  double fwhm = FWHM( sigma, gamma );
  double val=fwhm/2.3548;

  double efwhm = FWHMError(/*sigma, gamma,*/ esigma, egamma,
			   Vss, Vsg, Vgs, Vgg);
      
  cout<<val<<"   "<<efwhm/2.3548<<endl;
}


TGraphAsymmErrors* ComputeChristianUncertainties(string type, string comp,TGraph* mc) {

  TFile* scvfile=new TFile("/home/mmarionn/Documents/CMS/METStudies/Combination/ChristianSystematics/plotZllRecoilCorrection_pfMEtNoPileUpSmeared_beforeAddPUreweight_graphs.root","READ");

  string base;
  if(type=="resp")
    base = "uParl_div_qT_mean_mc_signal";
  if(type=="resoQt")
    base = ((comp=="para")?"uParl":"uPerp")+(string)"_rms_mc_signal";

  string bname="graph_"+base;

  TGraphAsymmErrors* mean=(TGraphAsymmErrors*)scvfile->Get(bname.c_str());

  TGraphAsymmErrors* uJES=(TGraphAsymmErrors*)scvfile->Get((bname+"_jetEnUp").c_str() );
  TGraphAsymmErrors* dJES=(TGraphAsymmErrors*)scvfile->Get((bname+"_jetEnDown").c_str() );
  TGraphAsymmErrors* uJER=(TGraphAsymmErrors*)scvfile->Get((bname+"_jetResUp").c_str() );
  TGraphAsymmErrors* dJER=(TGraphAsymmErrors*)scvfile->Get((bname+"_jetResDown").c_str() );
  // TGraphAsymmErrors* uMuo=;
  // TGraphAsymmErrors* dMuo=;
  TGraphAsymmErrors* uUnc=(TGraphAsymmErrors*)scvfile->Get((bname+"_unclEnUp").c_str() );
  TGraphAsymmErrors* dUnc=(TGraphAsymmErrors*)scvfile->Get((bname+"_unclEnDown").c_str() );

  TGraphAsymmErrors* errors=new TGraphAsymmErrors( mc->GetN() );
  double x,y;
  double eyh,eyl,jeh,jel,jrh,jrl,ueh,uel,tmp;
  cout<<bname<<endl;
  cout<<mean<<"  "<<uJES<<"  "<<dJES<<"   "<<uJER<<"   "<<dJER<<"   "<<uUnc<<"   "<<dUnc<<endl;
  double x2,y2,xm,ym,d=10000; int it=-1;

  for(int ii=0;ii<mc->GetN();ii++) {
    mc->GetPoint(ii,x2,y2);
    errors->SetPoint(ii,x2,y2);

    d=10000;
    for(int ip=0;ip<mean->GetN()-1;ip++) {

      mean->GetPoint(ip,x,y);

      if( fabs(x2-x) < d) {
	d=fabs(x2-x);
	it=ip;
	ym=y;
	cout<< " min "<<ip<<"   "<<x2<<"    "<<ym<<"   "<<d<<"   "<<it<<endl;
      }

    }

   

    uJES->GetPoint(it,tmp,jeh);
    dJES->GetPoint(it,tmp,jel);
    uJER->GetPoint(it,tmp,jrh);
    dJER->GetPoint(it,tmp,jrl);
    uUnc->GetPoint(it,tmp,ueh);
    dUnc->GetPoint(it,tmp,uel);

    cout<<" test "<<it<<"  "<<x2<<"  "<<tmp<<"  --> "<<jeh<<"   "<<jel<<"   "<<jrh<<"   "<<jrl<<"  " <<ueh<<"   "<<uel<<endl;
   
      
    // if( fabs(x2-x) < d) {
    // 	d=fabs(x2-x);
    // 	xm = x2;
    // 	ym = y2;
    // 	it=ii;
    // }
      
    
    eyh =  pow( (ym-jeh)/2. ,2);
    eyh += pow( (ym-jel)/2. ,2);
    eyh += pow( (ym-jrh)/2. ,2);
    eyh += pow( (ym-jrl)/2. ,2);
    eyh += pow( (ym-ueh)/2. ,2);
    eyh += pow( (ym-uel)/2. ,2);
    
    cout<<ii<<"   "<<x<<"   "<<x2<<"  " <<it<<"   "<<eyh<<"  "<<ym<<"   "<<y2<<endl;

    //if(it!=-1) {
    //errors->SetPoint(ii,xm,ym);
    errors->SetPointError(ii,0,0,sqrt(eyh)/ym*y2,sqrt(eyh)/ym*y2);
    //}
  }

  errors->SetName("unc");

  return errors;
}
