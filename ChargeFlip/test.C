#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <map>
#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TGraphAsymmErrors.h>
//#include <TLorentzVector.h>
#include <TMinuit.h>
#include <TKey.h>

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooArgSet.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumConvPdf.h>
#include <RooMsgService.h>

#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooPolynomial.h>

#include <Math/Functor.h>
#include <Fit/Fitter.h>



using namespace std;
using namespace RooFit;




void doFits(string file, string mcfile, bool isData, 
	    bool appendDb, string dbName, string singleCateg) {

  bool appDb=appendDb;

  TFile* f=new TFile(file.c_str(), "read");
  TFile* mcf=new TFile(mcfile.c_str(), "read");
  
  //scan the file content
  
  string name;
  vector<float> vs;
  map<string, vector<float> > vals;

  TIter nextkey(f->GetListOfKeys());
  TKey *key;
  while ((key = ((TKey*)nextkey()))) {
   
    TObject* obj = key->ReadObj(); 
   
    if( obj==nullptr ) continue;
    
    if(((string)obj->IsA()->GetName()).substr(0,5)=="TTree") {
      string name=(string)(obj->GetName());
      
      if(name.find("Pt_")==string::npos) continue;
      if(name.find("Fake")!=string::npos) continue;
      if(name.find("Prompt")!=string::npos) continue;
      if(singleCateg!="" && name.find(singleCateg)==string::npos) continue;

      //cout<<name<<" --->>> "<<f<<"  "<<mcf<<" -> "<<obj<<"  "<<f->Get( (name+"_Fake").c_str() )<<"  "<<mcf->Get( (name+"_Prompt").c_str() )<<endl;
      // TTree* promptTree = (TTree*)mcf->Get( (name+"_Prompt").c_str() );
      TTree* mainTree = (TTree*)obj;
      //TTree* fakeTree = (TTree*)f->Get( (name+"_Fake").c_str() );
      // cout<<" ====>>> "<<mainTree<<"  "<<fakeTree<<"   "<<promptTree<<"  "<<f<<endl;
      RooRealVar mZ("mZ","m_{ll}",50,120,"GeV");
   
      //cout<<" -> youp000 "<<mainTree->GetEntries()<<"  "<<fakeTree->GetEntries()<<"  "<<promptTree->GetEntries()<<endl;
      vs = doSingleFit("prout");//, nullptr, nullptr, nullptr, nullptr , true); // mainTree->GetName(), mainTree, fakeTree, promptTree, &mZ
      cout<<" -> youpi"<<endl;
      vals[ name ] = vs;
      cout<<" -> youpi2"<<endl;
      if(appendDb) appendDataBase(name, vs, appDb, dbName);
      if(appDb==false) appDb=true; //otherwise we overwrite the file
      cout<<" -> youpi3"<<endl;
      //delete promptTree, mainTree, fakeTree;
    }
    
  }
  cout<<"end of function "<<endl;
  // f->Close();
  // mcf->Close();


  cout<<"end /// "<<endl;
 
}



