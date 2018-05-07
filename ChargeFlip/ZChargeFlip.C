#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <map>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include "TMinuit.h"

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumConvPdf.h>

#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooPolynomial.h>

//#include "ProbaFit.h"

using namespace std;
using namespace RooFit;

vector<std::pair<vector<float>, vector<int> > > _vals;


void setPoint(float val, float eval, float p1, float p2) {

  vector<float> vals(2,0);
  vector<int> bins(2,0);
  
  vals[0]=val;
  vals[1]=eval;
  bins[0]=p1;
  bins[1]=p2;

  std::pair<vector<float>, vector<int> > p(vals,bins);
  _vals.push_back(p);

}


void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *param, Int_t iflag)
{
 
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

  //return chi2;
  f=chi2;
}




RooNumConvPdf* shapeZ(string tag, RooRealVar* x) {

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

  RooRealVar* mZ0=new RooRealVar( ("m_Z0_"+tag).c_str(),"Z0 mass", 91.188, "GeV/c^{2}" );
  RooRealVar* gammaZ0=new RooRealVar( ("gamma_Z0_"+tag).c_str(), "Z0 width",2.4952, "GeV/c^{2}" );
  RooBreitWigner* bw0=new RooBreitWigner( ("bw0"+tag).c_str(),"true BW",*x, *mZ0, *gammaZ0);
   
  RooRealVar* cb_bias=new RooRealVar( ("cbb_"+tag).c_str(), "bias",0.07, -3.0, 3.0 );
  RooRealVar* cb_width=new RooRealVar( ("cbw_"+tag).c_str(),"width", 1.,0.0,5 ); 
  RooRealVar* cb_alpha=new RooRealVar( ("cba_"+tag).c_str(),"alpha", 1.2,0.03,2.0 ); 
  RooRealVar* cb_power=new RooRealVar( ("cbn_"+tag).c_str(),"power", 5 ); 

  RooCBShape* cb_pdf=new RooCBShape( ("cb_pdf_"+tag).c_str(), "CB shape", 
				     *x,*cb_bias, *cb_width, *cb_alpha, *cb_power );
  
  RooNumConvPdf* bw=new RooNumConvPdf( ("bw_"+tag).c_str(),"Convolution", *x, *cb_pdf, *bw0 );

  return bw;
}


RooAbsPdf* shapeSB(string tag, RooRealVar* x, RooRealVar* nSig, RooDataSet* fakeDS) {


  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

  RooNumConvPdf* bw = shapeZ(tag, x);

  RooRealVar* p0=new RooRealVar( ("pol_p0_"+tag).c_str(), "p0", -10000, 10000., 0.);
  RooRealVar* p1=new RooRealVar( ("pol_p1_"+tag).c_str(), "p1", -10000, 10000., 0.);
  RooRealVar* p2=new RooRealVar( ("pol_p2_"+tag).c_str(), "p2", -10000, 10000., 0.);
  RooRealVar* p3=new RooRealVar( ("pol_p3_"+tag).c_str(), "p3", -10000, 10000., 0.);
  RooRealVar* p4=new RooRealVar( ("pol_p4_"+tag).c_str(), "p4", -10000, 10000., 0.);
  
  if(fakeDS!=NULL) {
    RooPolynomial* pol_pdf=new RooPolynomial( ("pol_pdf_"+tag).c_str(), "bkg shape", *x, RooArgList(*p0,*p1,*p2,*p3,*p4),0 );
  
    RooFitResult* result = pol_pdf->fitTo(*fakeDS,RooFit::SumW2Error(kFALSE) );
  
    p0->setConstant(kTRUE);
    p1->setConstant(kTRUE);
    p2->setConstant(kTRUE);
    p3->setConstant(kTRUE);
    p4->setConstant(kTRUE);
  }
  RooRealVar* n_bkg=new RooRealVar( ("N_{bkg} "+tag).c_str(),"n bkg events", fakeDS?(fakeDS->sumEntries()):0 );
  n_bkg->setConstant(kTRUE);
  

  RooArgList* listPdf=new RooArgList( *pol_pdf, *bw );
  RooArgList* listPdfVal=new RooArgList( *n_bkg, *nSig );
  RooAddPdf* bw_tot=new RooAddPdf( "bw_EBEB_MC", "PDF ee", *listPdf, *listPdfVal );

  return bw_tot;

}

vector<float> doSingleFit(TH1* histo, bool isData) {

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval); // 1 for INFO
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);

  RooRealVar mass("mass","m_{ll}",60,120,"GeV");
  RooRealVar ch("channel","channel",channel-0.1,channel+0.1);
  RooRealVar cat("categ","categ",categ-0.1,categ+0.1);
  RooRealVar w("w","w",std::numeric_limits<double>::min(),std::numeric_limits<double>::max());
  
  RooDataSet* data(0);
  // if(isData)
  //   data = new RooDataSet("data","data",tree,RooArgSet(mass,ch,cat));
  // else
  data = new RooDataSet("data","data",tree,RooArgSet(mass,ch,cat,w),"w==1","w");

  RooRealVar wf("w","w",std::numeric_limits<double>::min(),std::numeric_limits<double>::max());
  RooDataSet* dataFake = new RooDataSet("dataFake","dataFake",tree,RooArgSet(mass,ch,cat,w),"w!=1","w");

  cout<<"category "<<categ<<" -> "<<data->numEntries()<<"  "<<data->sumEntries()<<endl;
  vector<float> v(4,0);
  //return v;
  if(data->sumEntries() < 0.001 || data->numEntries() <10 ) {
    //cout<<"result:\t"<<categ<<"\t"<<data->numEntries()<<"\t"<<data->sumEntries()<<endl;
    
    
    RooDataSet* dataRed(0);
    RooDataSet* dataRedOut(0);
    if(isData) {
      dataRed= new RooDataSet("datared","datared",data,RooArgSet(mass,ch,cat),"mass<120 && mass >50");
      dataRedOut= new RooDataSet("datared","datared",data,RooArgSet(mass,ch,cat),"mass>120 || mass <50");
    }
    else {
      dataRed= new RooDataSet("datared","datared",data,RooArgSet(mass,ch,cat,w),"mass<120 && mass >50","weight");
      dataRedOut= new RooDataSet("datared","datared",data,RooArgSet(mass,ch,cat,w),"mass>120 || mass <50","weight");
    }

    //if(isData) {
    v[0]=dataRed->sumEntries();
    v[1]=sqrt( dataRed->numEntries()*dataRed->weight() );
    v[2]=dataRedOut->sumEntries();
    v[3]=sqrt( dataRed->numEntries()*dataRed->weight() );
      //}
    // else {

    // }
    return v;
  }

  RooRealVar nSig("nSig","nSig",data->sumEntries(),0,10000000);
  //RooRealVar nBkg("nBkg","nBkg",1.,0,10000000);

  ostringstream os;
  os<<categ;
  RooAbsPdf* shape=shapeSB( ("s"+os.str()),&mass, &nSig, dataFake);

 

  RooFitResult* result;
  if(data->numEntries() < 1000 ) { //reasonnable number of entries for a unbinned fit
    result = shape->fitTo( *data ,RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );
  }
  else {
    RooDataHist hist("hist","hist",mass,(TH1*)(data->createHistogram("mass",80))->Clone());
    result = shape->fitTo( hist ,RooFit::SumW2Error(kFALSE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );
  }

  double N=nSig.getVal();
  double eN=nSig.getError();

  double NB=nBkg.getVal();
  double eNB=nBkg.getError();

  cout<<"result:\t"<<categ<<"\t"<<N<<"\t"<<eN<<endl;

  TCanvas* c=new TCanvas( ("c"+os.str()).c_str(),("c"+os.str()).c_str());
  RooPlot* frame=mass.frame();
  data->plotOn(frame);
  shape->plotOn(frame);
  frame->Draw();

  string name="plots/fitData_";
  if(!isData) name="plots/fitMC_";
  c->SaveAs( (name+os.str()+".png").c_str() );

  delete frame;
  delete c;
  
  delete shape;
  delete data;


  v[0]=N;
  v[1]=eN;
  v[2]=NB;
  v[3]=eNB;
  
  return v;
}



void doFits(bool isData=false, int min=0, int max=288) {

  
  TChain* tree= new TChain("tree");
  
  if(isData) {
    tree->Add("/home/mmarionn/Documents/CMS/MPAF/workdir/skims/DoubleEG_Run2015D_v3_runs_256630_257599.root");
  }
  else {
    tree->Add("/home/mmarionn/Documents/CMS/MPAF/workdir/skims/DYJetsToLL_M50.root");
    tree->Add("/home/mmarionn/Documents/CMS/MPAF/workdir/skims/TT_pow.root");
  }
  
  //TFile* f=new TFile("/home/mmarionn/Documents/CMS/MPAF/workdir/skims/DYJetsToLL_M50.root","read");
  //TTree* tree=(TTree*)f->Get("tree");
  
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("mass",1);
  tree->SetBranchStatus("channel",1);
  tree->SetBranchStatus("categ",1);
  tree->SetBranchStatus("weight",1);

  cout<<"Number of entries "<<tree->GetEntries()<<endl;

  gROOT->SetBatch(kTRUE);

  ostringstream osm,osM;
  osm<<min;
  osM<<max;
  string name="fitResults.txt";
  if(!isData) name="fitResultsMC_"+osm.str()+"_"+osM.str()+".txt";
  ofstream ofile(name.c_str(),ios::out | ios::trunc);
  vector<float> v;
  for(int i=min;i<max;i++) {
    v = doSingleFit(i, 22, tree, isData);
    ofile<<i<<"\t"<<v[0]<<"\t"<<v[1]<<"\t"<<v[2]<<"\t"<<v[3]<<endl;
  }
  
  gROOT->SetBatch(kFALSE);
}




void getDataES() {

  TFile* dF=new TFile("/home/mmarionn/Documents/CMS/MPAF/workdir/skims/dataFull.root","read");
  // TFile* mcF=new TFile("DYJetsToLL_ZMass.root","read");

  TTree* dT=(TTree*)dF->Get("tree");
  // TTree* mcT=(TTree*)mcF->Get("Zee");
  
  RooRealVar cBB("mass","m_{ll} ",60,120, "GeV");
  RooRealVar cEE("mass","m_{ll} ",60,120, "GeV");

  RooRealVar chan("channel","channel",21,23);
  RooRealVar charge("charge","charge",-0.1,0.1);
  RooRealVar bb("fidu","fidu",1.9,2.1);
  RooRealVar ee("fidu","fidu",-0.1,0.1);

  RooDataSet dataBB("dataBB","dataBB",dT, RooArgSet(cBB,chan,bb,charge) ); // ,chan,bb,charge
  RooDataSet dataEE("dataEE","dataEE",dT, RooArgSet(cEE,chan,ee,charge) );

  // RooDataSet mcBB("mcBB","mcBB",mcT, RooArgSet(cBB) );
  // RooDataSet mcEE("mcEE","mcEE",mcT, RooArgSet(cEE) );
  
  RooDataHist datahBB("datahBB","datahBB", cBB, dataBB.createHistogram("mass" , 60) );
  RooDataHist datahEE("datahEE","datahEE", cEE, dataEE.createHistogram("mass" , 60) );

  // RooDataHist mchBB("mchBB","mchBB", cBB, mcBB.createHistogram("Z_massB" , 240) );
  // RooDataHist mchEE("mchEE","mchEE", cEE, mcEE.createHistogram("Z_massE" , 240) );

  RooRealVar nSigBB("nSigBB","nSig",500,0,10000000);
  RooRealVar nSigEE("nSigEE","nSig",500,0,10000000);
  RooRealVar nBkgBB("nBkgBB","nBkg",500,0,10000000);
  RooRealVar nBkgEE("nBkgEE","nBkg",500,0,10000000);

  RooAbsPdf* zDBB=shapeSB("dBB",&cBB, &nSigBB, NULL);
  RooAbsPdf* zDEE=shapeSB("dEE",&cEE, &nSigEE, NULL);
  
  // RooAddPdf* zMcBB=shapeSB("mcBB",&cBB);
  // RooAddPdf* zMcEE=shapeSB("mcEE",&cEE);

  //Data BB ====================================
  cout<<" data BB "<<zDBB<<"   "<<dataBB.numEntries()<<" /  "<<dataEE.numEntries()<<endl;

  RooFitResult* rDBB = zDBB->fitTo(datahBB,RooFit::SumW2Error(kFALSE) );
  TCanvas* cDBB =new TCanvas("cDBB","DBB",600,600);
  // cout<<" gloubi "<<endl;
  RooPlot* fDBB = cBB.frame() ;
  datahBB.plotOn(fDBB) ;
  zDBB->plotOn(fDBB) ;
  zDBB->plotOn(fDBB,Components("exp_pdf_dBB"),LineStyle(kDashed),LineColor(kRed)) ;
  zDBB->paramOn(fDBB,Layout(0.181,0.52,0.92) );
  fDBB->getAttText()->SetTextSize(0.026);
  fDBB->Draw();

  //Data EE ====================================
  cout<<" data EE "<<endl;
  RooFitResult* rDEE = zDEE->fitTo(datahEE,RooFit::SumW2Error(kFALSE) );
  TCanvas* cDEE =new TCanvas("cDEE","DEE",600,600);
   
  RooPlot* fDEE = cEE.frame() ;
  datahEE.plotOn(fDEE) ;
  zDEE->plotOn(fDEE) ;
  zDEE->plotOn(fDEE,Components("exp_pdf_dEE"),LineStyle(kDashed),LineColor(kRed)) ;
  zDEE->paramOn(fDEE,Layout(0.582,0.921,0.916) );
  fDEE->getAttText()->SetTextSize(0.026);
  fDEE->Draw();



  // //MC BB ====================================
  // cout<<" mc BB "<<endl;
  // RooFitResult* rMcBB = zMcBB->fitTo(mchBB,RooFit::SumW2Error(kFALSE) );
  // TCanvas* cMcBB =new TCanvas("cMcBB","McBB",600,600);
   
  // RooPlot* fMcBB = cBB.frame() ;
  // mcBB.plotOn(fMcBB) ;
  // zMcBB->plotOn(fMcBB) ;
  // zMcBB->paramOn(fMcBB);
  // fMcBB->Draw();

  // //MC EE ====================================
  // cout<<" mc EE "<<endl;
  // RooFitResult* rMcEE = zMcEE->fitTo(mchEE,RooFit::SumW2Error(kFALSE) );
  // TCanvas* cMcEE =new TCanvas("cMcEE","McEE",600,600);
   
  // RooPlot* fMcEE = cEE.frame() ;
  // mcEE.plotOn(fMcEE) ;
  // zMcEE->plotOn(fMcEE) ;
  // zMcEE->paramOn(fMcEE);
  // fMcEE->Draw();


  // //get the parameters ==================================
  // float biasDBB = zDBB->getParameters(datahBB)->getRealValue("cb_w_dBB");
  // float biasDEE = zDEE->getParameters(datahEE)->getRealValue("cb_w_dEE");
  // float biasMcBB = zMcBB->getParameters(mchBB)->getRealValue("cb_w_mcBB");
  // float biasMcEE = zMcEE->getParameters(mchEE)->getRealValue("cb_w_mcEE");

  // cout<<biasDBB<<"   "<<biasDEE<<"   "<<biasMcBB<<"   "<<biasMcEE<<endl;

  // float esBB = 91.18/(biasMcBB+91.18- biasDBB);
  // float esEE = 91.18/(biasMcEE+91.18- biasDEE);

  // cout<<" Barrel : "<<esBB<<endl;
  // cout<<" Endcap : "<<esEE<<endl;


}


float getErr(float m1, float m2, float e1, float e2) {

  float e=pow(e1/m1,2)+pow(e2/m2,2);
  if(m1==0) return 0;
  //cout<<m1<<"  "<<m2<<"  "<<e1<<"   "<<e2<<" ===>  "<<m1/m2<<"   "<<pow(e1/m1,2)<<"   "<<pow(e2/m2,2)<<" --->  "<<(m1/m2)*sqrt(e)<<endl;
  return (m1/m2)*sqrt(e);
}


void makeDB() {

  ifstream categs("categs.txt", ios::in);
  ifstream vals("fitResults.txt", ios::in);
  string line;
  
  map<int,vector<float> > catDef;
  map<int,vector<float> > cont;

  map<int,bool> skim;
  map<int,bool>::const_iterator sit;
  skim[144]=1;
  skim[157]=1;
  skim[170]=1;
  skim[183]=1;
  skim[196]=1;
  skim[209]=1;
  skim[222]=1;
  skim[235]=1;
  skim[248]=1;
  skim[261]=1;
  skim[274]=1;
  skim[287]=1;


  while(getline(categs, line)) 
    {
      istringstream iss(line);
      vector<string> tks;
      copy(istream_iterator<string>(iss),
	   istream_iterator<string>(),
	   back_inserter<vector<string> >(tks));
      
      if(tks[1]=="0") continue;
      
      vector<float> v;
      for(unsigned int i=2;i<tks.size();i++)
	v.push_back(atof(tks[i].c_str()) );
      catDef[ atoi(tks[0].c_str()) ] = v;
      //cout<<atoi(tks[0].c_str())<<endl;
    }

  vector<vector<float> > ssVals;
  int n=0,k=0,pN=0;

  vector<vector<float> > Xprobs(4,vector<float>(3,0));

  while(getline(vals, line)) 
    {
      //cout<<"gloubdwi"<<endl;
      
      istringstream iss(line);
      vector<string> tks;
      copy(istream_iterator<string>(iss),
	   istream_iterator<string>(),
	   back_inserter<vector<string> >(tks));
  
      //cout<<n<<"  "<<ssVals.size()<<endl;
      //cout<<n<<"  "<<k<<"   "<<tks[0]<<"   "<<catDef.size()<<endl;
      if(n>(catDef.size()-1)) {
	  
	float val=atof(tks[1].c_str());
	float eval=atof(tks[2].c_str());
	  
	float mean=0;
	float err=0;
	//cout<<"gloubi"<<endl;
	if(val!=0) {
	  mean=ssVals[k][0]/val;
	  err=getErr(ssVals[k][0], val, ssVals[k][1], eval);
	}	  
	//cout<<n<<"  "<<k<<"   "<<tks[0]<<"   "<<mean<<"  "<<val<<"   "<<ssVals[k][0]<<"   "<<err<<endl;
	vector<float> v(2,0);
	v[0]=mean;
	v[1]=err;


	// cout<<"gloubi1"<<endl;
	// cout<<n<<"  "<<k<<endl;
	cont[n] = v;
	sit=skim.find(n);
	if(sit!=skim.end()) { 
	  cout<<n<<" ===>  "<<mean/2.<<"   "<<err/2.<<"    ("<<val<<","<<ssVals[k][0]<<") "<<endl;

	  int eb=pN%3;
	  int pb=pN/3;//pN-eb*3;
	  //cout<<pN<<"  "<<pb<<"  "<<eb<<endl;
	  Xprobs[ pb ][ eb ]=mean/2.;
	  pN++;
	}
	k++;
	// cout<<"gloubi2"<<endl;
	// if(n==287) break;//boh...
      }
      else {
	//cout<<" coin? "<<endl;
	vector<float> v;
	for(unsigned int i=1;i<tks.size();i++)
	  v.push_back(atof(tks[i].c_str()) );

	//cout<<n<<"   "<<v[0]<<"   "<<v[1]<<endl;
	ssVals.push_back(v);
      }
   
      n++;
    }
  // return;
  //cout<<"pouet"<<endl;
  categs.close();
  vals.close();
  //cout<<"pouet"<<endl;
  // ofstream odb("db.txt", ios::out | ios::trunc);
  
  // map<int,vector<float> >::const_iterator it;
  // map<int,vector<float> >::const_iterator it2;
  // for(it = cont.begin();it!=cont.end();++it) {

  //   it2=catDef.find(it->first);
  //   for(size_t i=0;i<it2->second.size();i++)
  //     odb<<it2->second[i]<<"\t";
    
  
  //   odb<<it->second[0]<<"\t"<<it->second[1]<<endl;
  // }
 
  // //odb.write()
  // odb.close();


  float p00=0;
  float p01=1;
  float p02=2;
  float p10=3;
  float p11=4;
  float p12=5;
  float p20=6;
  float p21=7;
  float p22=8;
  float p30=9;
  float p31=10;
  float p32=11;

  map<int, pair<int,int> > assoc;
  assoc[144]=std::make_pair( p00, p00 );
  assoc[145]=std::make_pair( p00, p01 );
  assoc[146]=std::make_pair( p00, p02 );
  assoc[157]=std::make_pair( p01, p01 );
  assoc[158]=std::make_pair( p01, p02 );
  assoc[170]=std::make_pair( p02, p02 );
  assoc[180]=std::make_pair( p10, p00 );
  assoc[181]=std::make_pair( p10, p01 );
  assoc[182]=std::make_pair( p10, p02 );
  assoc[183]=std::make_pair( p10, p10 );
  assoc[184]=std::make_pair( p10, p11 );
  assoc[185]=std::make_pair( p10, p12 );
  assoc[192]=std::make_pair( p11, p00 );
  assoc[193]=std::make_pair( p11, p01 );
  assoc[194]=std::make_pair( p11, p02 );
  assoc[196]=std::make_pair( p11, p11 );
  assoc[197]=std::make_pair( p11, p12 );
  assoc[204]=std::make_pair( p12, p00 );
  assoc[205]=std::make_pair( p12, p01 );
  assoc[206]=std::make_pair( p12, p02 );
  assoc[209]=std::make_pair( p12, p12 );
  assoc[216]=std::make_pair( p20, p00 );
  assoc[217]=std::make_pair( p20, p01 );
  assoc[218]=std::make_pair( p20, p02 );
  assoc[219]=std::make_pair( p20, p10 );
  assoc[220]=std::make_pair( p20, p11 );
  assoc[221]=std::make_pair( p20, p12 );
  assoc[222]=std::make_pair( p20, p20 );
  assoc[223]=std::make_pair( p20, p21 );
  assoc[224]=std::make_pair( p20, p22 );
  assoc[228]=std::make_pair( p21, p00 );
  assoc[229]=std::make_pair( p21, p01 );
  assoc[230]=std::make_pair( p21, p02 );
  assoc[231]=std::make_pair( p21, p10 );
  assoc[232]=std::make_pair( p21, p11 );
  assoc[233]=std::make_pair( p21, p12 );
  assoc[235]=std::make_pair( p21, p21 );
  assoc[236]=std::make_pair( p21, p22 );
  assoc[240]=std::make_pair( p22, p00 );
  assoc[241]=std::make_pair( p22, p01 );
  assoc[242]=std::make_pair( p22, p02 );
  assoc[243]=std::make_pair( p22, p10 );
  assoc[244]=std::make_pair( p22, p11 );
  assoc[245]=std::make_pair( p22, p12 );
  assoc[248]=std::make_pair( p22, p22 );
  assoc[252]=std::make_pair( p30, p00 );
  assoc[253]=std::make_pair( p30, p01 );
  assoc[254]=std::make_pair( p30, p02 );
  assoc[255]=std::make_pair( p30, p10 );
  assoc[256]=std::make_pair( p30, p11 );
  assoc[257]=std::make_pair( p30, p12 );
  assoc[258]=std::make_pair( p30, p20 );
  assoc[259]=std::make_pair( p30, p21 );
  assoc[260]=std::make_pair( p30, p22 );
  assoc[261]=std::make_pair( p30, p30 );
  assoc[262]=std::make_pair( p30, p31 );
  assoc[263]=std::make_pair( p30, p32 );
  assoc[264]=std::make_pair( p31, p00 );
  assoc[265]=std::make_pair( p31, p01 );
  assoc[266]=std::make_pair( p31, p02 );
  assoc[267]=std::make_pair( p31, p10 );
  assoc[268]=std::make_pair( p31, p11 );
  assoc[269]=std::make_pair( p31, p12 );
  assoc[270]=std::make_pair( p31, p20 );
  assoc[271]=std::make_pair( p31, p21 );
  assoc[272]=std::make_pair( p31, p22 );
  assoc[274]=std::make_pair( p31, p31 );
  assoc[275]=std::make_pair( p31, p32 );
  assoc[276]=std::make_pair( p32, p00 );
  assoc[277]=std::make_pair( p32, p01 );
  assoc[278]=std::make_pair( p32, p02 );
  assoc[279]=std::make_pair( p32, p10 );
  assoc[280]=std::make_pair( p32, p11 );
  assoc[281]=std::make_pair( p32, p12 );
  assoc[282]=std::make_pair( p32, p20 );
  assoc[283]=std::make_pair( p32, p21 );
  assoc[284]=std::make_pair( p32, p22 );
  assoc[287]=std::make_pair( p32, p32 );

  // ProbaFit* fit=new ProbaFit();
  map<int, pair<int,int> >::const_iterator asIt;
  map<int,vector<float> >::const_iterator it;
  for(asIt=assoc.begin();asIt!=assoc.end();asIt++) {
    int key=asIt->first;
    it = cont.find(key);

    //fit->setPoint(it->second[0], it->second[1], asIt->second.first, asIt->second.second);
    setPoint(it->second[0], it->second[1], asIt->second.first, asIt->second.second);
  }

  


  TMinuit *gMinuit = new TMinuit(12); //initialize TMinuit with a maximum of 5 params
  gMinuit->SetFCN(chi2);

  Double_t arglist[300];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist , 1, ierflg);



// Set starting values and step sizes for parameters
//static
 double step[12] = {0.00001,0.000001,0.00001,0.000001,
		     0.00001,0.000001,0.00001,0.000001,
		     0.00001,0.000001,0.00001,0.000001};
  double variable[12] = {0.0001,0.0001,0.0001,0.0001,
			 0.0001,0.0001,0.0001,0.0001,
			 0.0001,0.0001,0.0001,0.0001};
  
  gMinuit->mnparm(0, "p0", variable[0], step[0], 0, 0, ierflg);
  gMinuit->mnparm(1, "p1", variable[1], step[1], 0, 0, ierflg);
  gMinuit->mnparm(2, "p2", variable[2], step[2], 0, 0, ierflg);
  gMinuit->mnparm(3, "p3", variable[3], step[3], 0, 0, ierflg);
  gMinuit->mnparm(4, "p4", variable[4], step[4], 0, 0, ierflg);
  gMinuit->mnparm(5, "p5", variable[5], step[5], 0, 0, ierflg);
  gMinuit->mnparm(6, "p6", variable[6], step[6], 0, 0, ierflg); return;
  gMinuit->mnparm(7, "p7", variable[7], step[7], 0, 0, ierflg);
  gMinuit->mnparm(8, "p8", variable[8], step[8], 0, 0, ierflg);
  gMinuit->mnparm(9, "p9", variable[9], step[9], 0, 0, ierflg);
  gMinuit->mnparm(10, "p10", variable[10], step[10], 0, 0, ierflg);
  gMinuit->mnparm(11, "p11", variable[11], step[11], 0, 0, ierflg);
 


  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1;//1
  gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);

// Print results
   Double_t amin, edm, errdef;
   Int_t nvpar, nparx, icstat;
   gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
   gMinuit->mnprin(12, amin);









}
