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
#include <TError.h>

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooArgSet.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooSimultaneous.h>
#include <RooNumConvPdf.h>
#include <RooFFTConvPdf.h>
#include <RooMsgService.h>

#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooPolynomial.h>


#include <Math/Functor.h>
#include <Fit/Fitter.h>

bool shutup=false;

using namespace std;
using namespace RooFit;



struct binSt{

  binSt() {};
  binSt(vector<float> pt, vector<float> eta) {
    binsPt=pt;
    binsEta=eta;

    nBPt=binsPt.size();
    nBEta=binsEta.size();
  }

  vector<float> binsPt;
  vector<float> binsEta;

  int nBPt;
  int nBEta;

  int getNProb() {return nBPt*nBEta;}
  int getNBPt() {return nBPt;}
  int getNBEta() {return nBEta;}

};

//special function to determine which fake pdf to use
bool
usePolPdf(const string& tag, float np, float nm) {

  //asymmetry
  if( (np-nm)/(np+nm) > 0.95 ) return true;
  if( (np-nm)/(np+nm) < -0.95 ) return false;

  if(tag=="SC") return false;

  std::vector<std::string> elems;
  std::stringstream ss(tag);
  std::string item;
  while (std::getline(ss, item, '_')) {
    elems.push_back(item);
  }
  
  float pt1=atoi(elems[1].c_str());
  float eta1=atoi(elems[3].c_str());
  float pt2=atoi(elems[5].c_str());

  if(pt1>10) return true;
  if(std::abs(eta1)>1.479) return true;
  
  if(pt2>20) return true;

  return false;
}


//template signal
RooAbsPdf*
getTemplateSignal(string tag, RooRealVar* x, TTree* tempTree) {
  TH1F* histo=new TH1F( ("hTemp_"+tag).c_str(), ("hTemp_"+tag).c_str(), 40, 50, 120); //80
  tempTree->Draw( ("mZ>>hTemp_"+tag).c_str(),"w" );
  RooDataHist* sigDH = new RooDataHist( ("TSig_"+tag).c_str(), "data", *x, histo);
  RooHistPdf* sigPdf = new RooHistPdf( ("sigPdf_"+tag).c_str(), 
				       ("sigPdf_"+tag).c_str(), *x, *sigDH, 0 );
  return sigPdf;
}


//=============================================================================
void makeFit(RooAbsPdf* pdf, RooAbsData* dataset, bool fixPar=false, vector<RooRealVar*> vars=vector<RooRealVar*>(0,NULL), int offset=0 ) {
  
  RooFitResult* res = pdf->fitTo( *dataset,RooFit::SumW2Error(kFALSE),RooFit::PrintEvalErrors(-1) );
  //mares->minNll() 
  if(fixPar) {
    for(size_t iv=offset;iv<vars.size();iv++) { //skipping the four firsts
      vars[iv]->setConstant(kTRUE);
    }
  }

}


void add2Pdf(string tag, RooAbsPdf* pdf1, RooAbsPdf* pdf2, RooRealVar* y1, RooRealVar* y2, RooAddPdf*& outPdf) {

  RooArgList* listPdf=new RooArgList( *pdf1, *pdf2 );
  RooArgList* listPdfVal=new RooArgList( *y1, *y2 );
  outPdf=new RooAddPdf( tag.c_str(), "PDF ee", *listPdf, *listPdfVal );

}


void add3Pdf(string tag, RooAbsPdf* pdf1, RooAbsPdf* pdf2, RooAbsPdf* pdf3, 
	     RooRealVar* y1, RooRealVar* y2, RooRealVar* y3,
	     RooAddPdf*& outPdf) {

  RooArgList* listPdf=new RooArgList( *pdf1, *pdf2, *pdf3 );
  RooArgList* listPdfVal=new RooArgList( *y1, *y2, *y3 );
  outPdf=new RooAddPdf( tag.c_str(), "PDF ee", *listPdf, *listPdfVal );

}


RooAbsPdf* shapeZ(string tag, RooRealVar* x,
		  vector<RooRealVar*>& vars, RooAbsPdf* extPdf=nullptr) {

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);

  RooAbsPdf* bw0=nullptr;
  
  float width=3;
  if(extPdf==nullptr) {
    RooRealVar* mZ0=new RooRealVar( ("m_Z0_"+tag).c_str(),"Z0 mass", 91.188, "GeV/c^{2}" );
    RooRealVar* gammaZ0=new RooRealVar( ("gamma_Z0_"+tag).c_str(), "Z0 width",2.4952, "GeV/c^{2}" );
    RooBreitWigner* bw=new RooBreitWigner( ("bw0"+tag).c_str(),"true BW",*x, *mZ0, *gammaZ0);
    bw0=bw;
  } else {
    bw0 = extPdf;
    width=1;
  }   

  RooRealVar* cb_bias=new RooRealVar( ("cbb_"+tag).c_str(), "bias",3, -10.0, 10.0 );
  RooRealVar* cb_width=new RooRealVar( ("cbw_"+tag).c_str(),"width", width,0.00,30 ); //1
  RooRealVar* cb_alpha=new RooRealVar( ("cba_"+tag).c_str(),"alpha", 0.07,0.02,3.0 ); 
  RooRealVar* cb_power=new RooRealVar( ("cbn_"+tag).c_str(),"power", 5 ); 

  RooCBShape* cb_pdf=new RooCBShape( ("cb_pdf_"+tag).c_str(), "CB shape", 
				     *x,*cb_bias, *cb_width, *cb_alpha, *cb_power );
  //RooNumConvPdf* bw=new RooNumConvPdf( ("bw_"+tag).c_str(),"Convolution", *x, *cb_pdf, *bw0 );
  
  x->setBins(50000,("cache"+tag).c_str() ) ;
  RooFFTConvPdf*  bw=new RooFFTConvPdf( ("bw_"+tag).c_str(),"Convolution", *x, *bw0, *cb_pdf );

  vars.push_back(cb_bias);
  vars.push_back(cb_width);
  vars.push_back(cb_alpha);
  vars.push_back(cb_power);
  
  return (RooAbsPdf*)bw;
}


//=============================================================================

//Prompt constribution
RooAbsPdf* 
getPromptInFake(string tag, TTree* tree, RooRealVar* x, double& yield,  vector<RooRealVar*>& vars, bool temppdf=false) {
  
  RooAbsPdf* promptShape=nullptr;
  RooRealVar* n_Z=new RooRealVar( ("N_{sig}_p "+tag).c_str(),"n Z events",200000, 0., 10000000.);
  RooRealVar* n_bkg=new RooRealVar( ("N_{bkg}_p "+tag).c_str(),"n bkg events", 10., 0., 6000000.);
  TH1F* histo=new TH1F( ("hP_"+tag).c_str(), ("hP_"+tag).c_str(), 40, 50, 120); //80
  tree->Draw( ("mZ>>hP_"+tag).c_str(),"w" );
  RooDataHist* histPrompt=new RooDataHist( ("promptDs"+tag).c_str(), ("promptDs"+tag).c_str(), *x, histo );

  if(!temppdf) {
    
    RooAbsPdf* bw = shapeZ(tag, x, vars);
  
    RooRealVar* exp_tau=new RooRealVar( ("exp_tau_p_"+tag).c_str(), "tau", -0.05, -40., 0.);
    RooExponential* exp_pdf=new RooExponential( ("exp_pdf_p_"+tag).c_str(), "bkg shape", *x, *exp_tau );

    vars.push_back(exp_tau);
    vars.push_back(n_Z);

  
    RooAddPdf* tmpprtpdf;
    add2Pdf("prompt_"+tag, exp_pdf, bw, n_bkg, n_Z, tmpprtpdf);
    promptShape=tmpprtpdf;
  
 
    makeFit(promptShape, histPrompt , true, vars,3);
    if(n_bkg->getVal()/(n_Z->getVal()+n_bkg->getVal()) >0.9 )
      makeFit(promptShape, histPrompt , true, vars,0);
  
    cout<<" prompt bin : "<<histPrompt->sumEntries()<<" <bkg> "<<n_bkg->getVal()<<" <Z> "<<n_Z->getVal()<<endl;
  } else {

    RooAbsPdf* tmpPrtShape = getTemplateSignal("prt"+tag,x, tree);
    promptShape = shapeZ("Zprt"+tag, x, vars, tmpPrtShape);
    makeFit(promptShape, histPrompt , true, vars,3);
    n_Z->setVal(0);
    n_bkg->setVal(histo->Integral(0,100000) );
  }

  n_bkg->setConstant(kTRUE);
  vars.push_back(n_bkg);
  
  TCanvas* c=new TCanvas( "c","c");
  RooPlot* frame=x->frame();
  histPrompt->plotOn(frame);
 
  promptShape->plotOn(frame);
  //promptShape->plotOn(frame,RooFit::Components(*exp_pdf),RooFit::LineStyle(kDashed), RooFit::LineColor(kRed+1));
  if(tag=="SC")
    promptShape->paramOn(frame);
  frame->Draw();
  
  if(tree->GetEntries()>10) {
    c->SaveAs( ("plots/promptInFake_"+tag+".png").c_str() );
    c->SaveAs( ("plots/promptInFake_"+tag+".pdf").c_str() );
    c->SaveAs( ("plots/promptInFake_"+tag+".root").c_str() );
  }

  //delete frame; delete c;
  yield=vars.back()->getVal() + n_Z->getVal();
 
  if(tree->GetEntries()<10) yield=0;

  // delete histPrompt;
  return promptShape;
}

//Fake contribution
RooAbsPdf* 
getFakes(string tag, TTree* tree, RooRealVar* x, RooAbsPdf* promptPdf, double promptYield, double& fakeYield) {

  //get sum of weights
  float sumWgts=0;
  float iw,mZ;
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("w",1);
  tree->SetBranchStatus("mZ",1);
  tree->SetBranchAddress("w",&iw);
  tree->SetBranchAddress("w",&iw);
  float np=0,nm=0;
  for(int ie=0;ie<tree->GetEntries();ie++) {
    tree->GetEntry(ie);
    sumWgts+=iw;
    if(mZ>85)
      np+=iw;
    else
      nm+=iw;
  }
  tree->SetBranchStatus("*",1);

  vector<RooRealVar*> vars;
  cout<<" ========================= prompt yield : "<<promptYield<<"  entries="<<tree->GetEntries()<<"/ sumW="<<sumWgts<<" ==> fakes="<<(sumWgts - promptYield)<<"  /// "<<tag<<endl; //*0.1393807*0.62
  RooRealVar* p0=new RooRealVar( ("pol_p0_"+tag).c_str(), "p0", -1, -100, 100.); //-1
  RooRealVar* p1=new RooRealVar( ("pol_p1_"+tag).c_str(), "p1", 0.05, -10, 10.); //0.05
  RooRealVar* p2=new RooRealVar( ("pol_p2_"+tag).c_str(), "p2", 0.01, -10, 10.); //0.01
  RooRealVar* p3=new RooRealVar( ("pol_p3_"+tag).c_str(), "p3", 0.007, -1, 1.); //0.007
  RooRealVar* p4=new RooRealVar( ("pol_p4_"+tag).c_str(), "p4", -0.0002, -1, 1.); //-0.0002
  RooChebychev* pol_pdf=new RooChebychev( ("pol_pdf_"+tag).c_str(), "bkg shape", *x, RooArgList(*p0,*p1,*p2,*p3,*p4));

  RooRealVar* exp_tau=new RooRealVar( ("exp_tau_f_"+tag).c_str(), "tau", -0.05, -40., 0.);
  RooExponential* exp_pdf=new RooExponential( ("exp_pdf_f_"+tag).c_str(), "bkg shape", *x, *exp_tau );

  bool usePol=usePolPdf(tag, np,nm); //(tag.find("l1Pt_10")==string::npos && tag.find("l2Pt_10")==string::npos) || tag=="SC";
  RooAbsPdf* bk_pdf=usePol?((RooAbsPdf*)pol_pdf):((RooAbsPdf*)exp_pdf);

  RooRealVar* nBkg=new RooRealVar( ("N_{bkg} "+tag).c_str(),"n bkg events", std::max(0.,sumWgts - promptYield), 0, sumWgts );// ,0, sumWgts );// *0.1393807 *0.1393807*0.62
  RooRealVar* nPrompt=new RooRealVar( ("N_{prompt} "+tag).c_str(),"n prompt events", promptYield*0.65, promptYield*0.3, promptYield);//, promptYield*0.10, promptYield*1.1 ); //*0.1393807 //0.9636 *0.1393807*0.62
  
  vars.push_back(p0);
  vars.push_back(p1);
  vars.push_back(p2);
  vars.push_back(p3);
  vars.push_back(p4);
  vars.push_back(exp_tau);
  vars.push_back(nBkg);

  RooAddPdf* bkgShapeTmp=nullptr;
  RooAbsPdf* bkgShape=nullptr;
  if(promptYield!=0) {
    cout<<promptYield<<" <>  "<<sumWgts<<endl;
    if(promptYield<sumWgts) {
      cout<<" coin coin "<<endl;
      add2Pdf("bkg_"+tag, bk_pdf, promptPdf, nBkg, nPrompt, bkgShapeTmp);
      bkgShape=bkgShapeTmp;
    } else {
      bkgShape = promptPdf;
    }
    
  }
  else
    bkgShape=bk_pdf;

  // cout<<" input "<<vars[0]->getVal()<<" "<<vars[1]->getVal()<<" "<<vars[2]->getVal()<<" "<<vars[3]->getVal()<<" "<<vars[4]->getVal()<<" "<<vars[5]->getVal()<<" "<<vars[6]->getVal()<<" "<<usePol<<endl;

  RooRealVar w("w","w",0,100000,"w");
  RooDataSet* dataFake=nullptr;
  RooDataHist* dataFakeH=nullptr;
  if(tree->GetEntries()<2000) {
    dataFake= new RooDataSet( ("fake_"+tag).c_str(), "data",
			      tree, RooArgSet(*x, w),"","w");
    //cout<<" fake unbin : "<<dataFake->sumEntries()<<"  "<<dataFake->numEntries()<<" // "<<promptYield<<endl;
    makeFit(bkgShape, dataFake, true, vars);
  } else {
    TH1F* histo=new TH1F( ("hF_"+tag).c_str(), ("hF_"+tag).c_str(), 40, 50, 120); //80
    tree->Draw( ("mZ>>hF_"+tag).c_str(),"w" );
    dataFakeH = new RooDataHist( ("fake_"+tag).c_str(), "data", *x, histo);
    //cout<<" fake binned : "<<dataFakeH->sumEntries()<<"  "<<tree->GetEntries()<<" // "<<promptYield<<endl;
    makeFit(bkgShape, dataFakeH, true, vars);
  }
  // cout<<"   --> "<<vars[0]->getVal()<<" "<<vars[1]->getVal()<<" "<<vars[2]->getVal()<<" "<<vars[3]->getVal()<<" "<<vars[4]->getVal()<<" "<<vars[5]->getVal()<<" "<<vars[6]->getVal()<<" "<<endl;
 

  TCanvas* c=new TCanvas( "c2","c2");
  RooPlot* frame=x->frame();
  if(tree->GetEntries()<2000)
    dataFake->plotOn(frame);
  else
    dataFakeH->plotOn(frame);
  bkgShape->plotOn(frame);
  bkgShape->plotOn(frame, RooFit::Components(*promptPdf),RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen+1) );
  bkgShape->plotOn(frame, RooFit::Components(*bk_pdf),RooFit::LineStyle(kDashed), RooFit::LineColor(kRed+1) );
  if(tag=="SC")
    bkgShape->paramOn(frame);
  frame->Draw();
  if(tree->GetEntries()>20 || tag=="SC") {
    c->SaveAs( ("plots/fakes_"+tag+".png").c_str() );
    c->SaveAs( ("plots/fakes_"+tag+".pdf").c_str() );
    c->SaveAs( ("plots/fakes_"+tag+".root").c_str() );
  }
  delete frame, c;

  //delete dataFake, dataFakeH;
  fakeYield=nBkg->getVal();
  cout<<" ======================= fake yield "<<fakeYield<<" ("<<nPrompt->getVal()<<" / "<<promptYield<<")"<<endl;

  if(tree->GetEntries()<20) fakeYield=0;
  if(promptYield==0) fakeYield=(dataFake?dataFake->sumEntries():dataFakeH->sumEntries());
  if(fakeYield<=0) fakeYield=0;
  if(fakeYield<2) {
    RooRealVar* p0f=new RooRealVar( ("pol_p0f_"+tag).c_str(), "p0f", fakeYield);
    RooPolynomial* ppdf=new RooPolynomial( ("polpdf"+tag).c_str(),("polpdf"+tag).c_str(), *x, RooArgSet(*p0f),0 );
  }

  return bk_pdf;

}


//final contribution 



vector<RooRealVar*>
getYield(string tag, RooAbsData* data, RooRealVar* x, RooAbsPdf* fakePdf, 
	 float bkgYield, RooAbsPdf*& tt_pdf, RooAbsPdf*& Z_pdf, RooAbsPdf*& mainPdf, bool validTT) {
  cout<<" ======================= data yields --> "<<data->sumEntries()<<"  bkg="<<bkgYield<<" -> expSig="<<data->sumEntries()-bkgYield<<" / "<<validTT<<endl;
  vector<RooRealVar*> vars;
  RooAbsPdf* sigZ = shapeZ("Z"+tag, x, vars, Z_pdf);
  //RooNumConvPdf* tt_pdf = shapeZ(tag, x, vars, tt_pdf);
  tt_pdf = shapeZ("tt"+tag, x, vars, tt_pdf);


  RooRealVar* nTT=new RooRealVar( ("N_{tt} "+tag).c_str(),"n bkg events", std::max(0.,(data->sumEntries()-bkgYield)*0.1),0, data->sumEntries() );
  RooRealVar* nZ=new RooRealVar( ("N_{Z} "+tag).c_str(),"n prompt events", std::max(0.,(data->sumEntries()-bkgYield)), 0, data->sumEntries() );

  // RooAddPdf* sigshape=nullptr;
  // add2Pdf("tot_"+tag, tt_pdf, sigZ, nTT, nZ, sigshape );
 
  RooRealVar* nBkg=new RooRealVar( ("N_{bkg}"+tag).c_str(),"n bkg events", bkgYield );
  //RooRealVar* nSig=new RooRealVar( ("N_{sig}"+tag).c_str(),"n sig events", data->sumEntries()-bkgYield, 0., 1000000. );
  RooAddPdf* tmpshape=nullptr;
  RooAbsPdf* shape=nullptr;
  if(bkgYield>1)
    add3Pdf("tot_"+tag, fakePdf, tt_pdf, sigZ, nBkg, nTT, nZ, tmpshape );
  else if(validTT)
    add2Pdf("tot_"+tag, tt_pdf, sigZ, nTT, nZ, tmpshape );
  else
    shape=sigZ;

  if(tmpshape) shape=tmpshape;

  makeFit(shape, data, false);

  mainPdf=(RooAbsPdf*)shape;
  vector<RooRealVar*> out({nZ,nTT,nBkg});
  return out;
}



//=============================================================================

vector<float> doSingleFit(string tag, TTree* mainTree, TTree* fakeTree, 
			  TTree* promptTree, TTree* tempTreeZ, TTree* tempTreeTT, 
			  RooRealVar* x, bool isData, bool cutAndCount) {
  cout<<" ========================= "<<tag<<" =============================== "<<endl;
  vector<float> v(4,0);

  //temporary
  float factor=1;

  if(!shutup) {
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval); // 1 for INFO
    RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::Caching);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
    RooMsgService::instance().getStream(0).removeTopic(RooFit::InputArguments);
    RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);
    RooMsgService::instance().setSilentMode(true); 
    shutup=true;
  }  
 
  double N, eN, NTT, eNTT, NB, eNB;
  if(cutAndCount) {
    N=mainTree->GetEntries();
    eN=sqrt(mainTree->GetEntries());
  
    NTT=promptTree->GetEntries();
    eNTT=sqrt(promptTree->GetEntries());

    NB=fakeTree->GetEntries();
    eNB=sqrt(fakeTree->GetEntries());

    //test with no fakes
    TH1F* histo=new TH1F( ("hP_"+tag).c_str(), ("hP_"+tag).c_str(), 40, 50, 120); //80
    promptTree->Draw( ("mZ>>hP_"+tag).c_str(),"w" );
    TH1F* histoF=new TH1F( ("hF_"+tag).c_str(), ("hF_"+tag).c_str(), 40, 50, 120); //80
    fakeTree->Draw( ("mZ>>hF_"+tag).c_str(),"w" );
    N+=histo->Integral()-histoF->Integral();
  }
  if(!cutAndCount) {

    //first, prompt contribution
    double promptYield;vector<RooRealVar*> vars;
    RooAbsPdf* promptPdf=getPromptInFake(tag, promptTree, x, promptYield, vars, promptTree->GetEntries()>1000);
    //return v;
    // //second, fake contribution
    double fakeYield;
    RooAbsPdf* fakePdf=getFakes(tag, fakeTree, x, promptPdf, promptYield, fakeYield);

    for(int i=0;i<vars.size();i++) {
      vars[i]->setConstant(kTRUE);
    }

    cout<<" checking MC stats : "<<tempTreeZ->GetEntries()<<" / "<<tempTreeTT->GetEntries()<<endl;
    //third, get the template from simulation
    RooAbsPdf* sigPdf_tt = (tempTreeTT->GetEntries()>5)?getTemplateSignal("TT"+tag, x, tempTreeTT):nullptr;
    RooAbsPdf* sigPdf_Z  = (tempTreeZ->GetEntries()>100)?getTemplateSignal("Z"+tag, x, tempTreeZ):nullptr;
  
    //and final fit
    RooRealVar w("w","w",0,100000,"w");
    RooDataSet* data=nullptr;
    RooDataHist* dataH=nullptr;
    if(mainTree->GetEntries()<2000) {
      data= new RooDataSet( ("data_"+tag).c_str(), "data",
			    mainTree, RooArgSet(*x, w),"","w");
      //cout<<" main unbinned "<<data->sumEntries()<<endl;
    } else {
      TH1F* histo=new TH1F( ("h_"+tag).c_str(), ("h_"+tag).c_str(), 40, 50, 120); //80
      mainTree->Draw( ("mZ>>h_"+tag).c_str(),"w" );
      dataH = new RooDataHist( ("data_"+tag).c_str(), "data", *x, histo);
      //cout<<" main binned "<<dataH->sumEntries()<<endl;
    }
    cout<<" ============================>>>>>>>>>>>>>>>> "<<promptYield<<" <> "<<fakeYield<<" (factor "<<factor<<")"<<endl;
    RooAbsPdf* mainPdf; RooAbsPdf* ttPdf;
    vector<RooRealVar*> yields=getYield(tag, (dataH?(RooAbsData*)dataH:(RooAbsData*)data),
					x, fakePdf, fakeYield, 
					sigPdf_tt, sigPdf_Z, mainPdf, sigPdf_tt!=nullptr);
  
    N=yields[0]->getVal()/factor;
    eN=yields[0]->getError()/factor; if(eN==0) eN=sqrt(N);
  
    NTT=yields[1]->getVal()/factor;
    eNTT=yields[1]->getError()/factor;

    NB=yields[2]->getVal()/factor;
    eNB=yields[2]->getError()/factor;
  
    cout<<"result:\t"<<tag<<"\t nZ= "<<N<<" +- "<<eN<<"   (nTT= "<<NTT<<" +- "<<eNTT<<" )"<<"   (nbkg= "<<NB<<" +- "<<eNB<<" )"<<endl;
  
    string os=tag;
    TCanvas* c=new TCanvas( ("c"+os).c_str(),("c"+os).c_str());
    RooPlot* frame=x->frame();

    if(data)
      data->plotOn(frame,Binning(20));
    else
      dataH->plotOn(frame);

    mainPdf->plotOn(frame);
    mainPdf->plotOn(frame,RooFit::Components( ("pol_pdf_"+tag).c_str() ),RooFit::LineColor(kGray+3),RooFit::LineStyle(kDashed)); //fakePdf
    mainPdf->plotOn(frame,RooFit::Components( ("exp_pdf_f_"+tag).c_str() ),RooFit::LineColor(kGray+3),RooFit::LineStyle(kDashed)); //fakePdf
    mainPdf->plotOn(frame,RooFit::Components( RooArgSet(*fakePdf,*sigPdf_tt) ), RooFit::LineColor(kRed+1),RooFit::LineStyle(kDashed));
    mainPdf->plotOn(frame,RooFit::Components( ("bw_Z"+tag).c_str() ), RooFit::LineColor(kGreen+2),RooFit::LineStyle(kDashed));
    if(tag=="SC")
      mainPdf->paramOn(frame);
    frame->Draw();

    FILE *test=fopen( "plots", "r" );
    if( test==0 ) system( "mkdir plots");
    else fclose( test );
    string name="plots/fitData_";
    if(!isData) name="plots/fitMC_";
    if(mainTree->GetEntries()>0) {
      c->SaveAs( (name+os+".png").c_str() );
      c->SaveAs( (name+os+".pdf").c_str() );
      c->SaveAs( (name+os+".root").c_str() );
    }
    delete c;
  }
  cout<<" ======================================================================"<<N<<endl;
  v[0]=N;
  v[1]=eN;
  v[2]=NB;
  v[3]=eNB;
  
  if(v[0]<0) { 
    v[0]=0;
    v[1]=0;
    v[2]=0;
    v[3]=0;
  }

 
  return v;
}


void appendDataBase(string name, vector<float> vs, bool appendDb, string dbname) {
  ofstream ofile(dbname.c_str(), ios::out | (appendDb?(ios::app):(ios::trunc)) );
  ofile<<name<<"\t"<<vs[0]<<"\t"<<vs[1]<<"\t"<<vs[2]<<"\t"<<vs[3]<<endl;
}


map<string, vector<float> > doFits(string file, string mcfile, bool isData, 
				   bool appendDb, string dbName, string singleCateg, bool cutAndCount) {

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
      if(name.find("TempPrompt")!=string::npos) continue;
      if(singleCateg!="" && name.find(singleCateg)==string::npos) continue;

      TTree* promptTree = (TTree*)mcf->Get( (name+"_Prompt").c_str() );
      TTree* tempTreeZ  = (TTree*)mcf->Get( (name+"_PromptZ").c_str() );
      TTree* tempTreeTT = (TTree*)mcf->Get( (name+"_PromptTT").c_str() );
      TTree* mainTree = (TTree*)obj;
      TTree* fakeTree = (TTree*)f->Get( (name+"_Fake").c_str() );
      cout<<mainTree->GetEntries()<<"  "<<fakeTree->GetEntries()<<"  "<<promptTree->GetEntries()<<endl;
      RooRealVar mZ("mZ","m_{ll}",50,120,"GeV");
      string tag=(singleCateg!="" && !appendDb)?"SC":mainTree->GetName();
      vs = doSingleFit(tag, mainTree, fakeTree, 
		       promptTree, tempTreeZ, tempTreeTT, &mZ, true, cutAndCount);
      vals[ name ] = vs;
    
      if(appendDb) appendDataBase(name, vs, appDb, dbName);
      if(appDb==false) appDb=true; //otherwise we overwrite the file
    
      //delete promptTree, mainTree, fakeTree;
    }
    
  }

  f->Close();
  mcf->Close();

  return vals;

}


vector<float> parseCateg(string cat) {
  
  std::replace( cat.begin(), cat.end(), '_', ' ');
  istringstream iss(cat);
  vector<string> tks;
  copy(istream_iterator<string>(iss),
       istream_iterator<string>(),
       back_inserter<vector<string> >(tks));

  //cout<<cat<<endl;
  int offset=1-(tks.size()==10);

  float pt1=atof( tks[2-offset].c_str() );
  float pt2=atof( tks[6-offset].c_str() );
  float eta1=atof( tks[4-offset].c_str() );
  float eta2=atof( tks[8-offset].c_str() );
  int os=atoi( tks[9-offset].c_str() );

  float tmp;
  // if(pt2>pt1) { //re-ordering if needed for pt
  //   tmp=pt2;
  //   pt2=pt1;
  //   pt1=tmp;
    
  // }
  // if(pt2==pt1 && eta2<eta1) { //re-ordering if needed for eta
  //   tmp=eta2;
  //   eta2=eta1;
  //   eta1=tmp;
  // }

  vector<float> bin(5,0);
  bin[0]=pt1;
  bin[1]=eta1;
  bin[2]=pt2;
  bin[3]=eta2;
  bin[4]=os;

  return bin;
}



float getErr(float m1, float m2, float e1, float e2) {

  float e=pow(e1/m1,2)+pow(e2/m2,2);
  if(m1==0) return 1;
  return (m1/m2)*sqrt(e);
}


// function Object to be minimized
struct Chi2 {
  
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

  // implementation of the function to be minimized
  double operator() (const double * param) {
    double chi2=0;

    float val,eval;
    int p1,p2;
    for(unsigned int ip=0;ip<_vals.size();ip++) {
      val=_vals[ip].first[0];
      eval=_vals[ip].first[1];
      p1=_vals[ip].second[0];
      p2=_vals[ip].second[1];
      
      float prob=param[p1]+param[p2]-param[p1]*param[p2];
      if(eval==0) eval=val;
      //chi2 += pow( val-(param[p1]+param[p2]), 2)/pow(eval,2);
      chi2 += pow( val-(prob), 2)/pow(eval,2);
    }
    //cout<<" chi2: "<<chi2<<endl;
    return chi2;
  }

};

void fillPoints(vector<float>& binsPt, vector<float>& binsEta,
		const map<string, float> yields, const map<string, float> eyields,
		const map<string, vector<float> >& bins, Chi2& chi2, 
		map<int,float[2]>& fixProbs ) {

  std::sort(binsPt.begin(), binsPt.end() );
  std::sort(binsEta.begin(), binsEta.end() );
  

  //build the number of probabilities needed
  //  int nProb=binsPt.size()*binsEta.size(); //lower boundary exists but not higher boundary
  int nBEta=binsEta.size();

  //then, fill
  int bPt1, bEta1, bPt2, bEta2;
  int p1, p2;
  float prob, eprob, yd, yn, eyn, eyd;
  map<string, vector<float> >::const_iterator itv;
  map<string, float>::const_iterator it,it2;
  int np=0;
  for(itv=bins.begin();itv!=bins.end();++itv) {
    
    if(itv->first.find("SS")!=string::npos) continue; //only looping over the OS categs

    //super ugly...
    bPt1=-1;bPt2=-1;
    for(size_t i=0;i<binsPt.size()-1;++i) {
      if(itv->second[0]>=binsPt[i] && itv->second[0]<binsPt[i+1])
	{ bPt1=i; }
      if(itv->second[2]>=binsPt[i] && itv->second[2]<binsPt[i+1])
	{ bPt2=i; }
    }
    if(bPt1==-1) bPt1=binsPt.size()-1;
    if(bPt2==-1) bPt2=binsPt.size()-1;

    bEta1=-1;bEta2=-1;
    for(size_t i=0;i<binsEta.size()-1;++i) {
      if(itv->second[1]>=binsEta[i] && itv->second[1]<binsEta[i+1])
	{ bEta1=i; }
      if(itv->second[3]>=binsEta[i] && itv->second[3]<binsEta[i+1])
	{ bEta2=i; }
    }
  
    if(bEta1==-1) bEta1=binsEta.size()-1;
    if(bEta2==-1) bEta2=binsEta.size()-1;

    p1=bPt1*nBEta + bEta1;
    p2=bPt2*nBEta + bEta2;
    
    it=yields.find( itv->first );
    it2=yields.find( (itv->first.substr(0, itv->first.size()-2)+"SS") );
    yd=it->second;
    yn=it2->second;

  
    
    it=eyields.find( itv->first );
    it2=eyields.find( (itv->first.substr(0, itv->first.size()-2)+"SS") );
    eyd=it->second;
    eyn=it2->second;

    prob=yn/(yn+yd);
    eprob=getErr(yn, yd, eyn, eyd);

    // cout<<itv->first<<" ===>>> ("<<itv->second[0]<<"/"<<itv->second[1]<<" -> "<<bPt1<<"/"<<bEta1<<")    ("
    // 	<<itv->second[2]<<"/"<<itv->second[3]<<" -> "<<bPt2<<"/"<<bEta2<<")"<<" --->> "<<yd<<"  "<<p1<<"/"<<p2<<endl;

    //FIXME
    if(yd==0) {prob=0; eprob=1;}
    
    
    // if(p1==p2) {
    //   cout<<" CHECK probability ("<<itv->second[0]<<","<<itv->second[1]
    // 	  <<") : "<<yn<<"+-"<<eyn<<" / "<<yd<<"+-"<<eyd<<" ==> "<<prob<<" +- "<<eprob<<endl;
    // }

    cout<<" CHECK probability ("<<setw(3)<<itv->second[2]<<","<<setw(5)<<itv->second[3]
     	<<" / "<<setw(3)<<itv->second[0]<<","<<setw(5)<<itv->second[1]
     	<<") : "<<setw(8)<<yn<<"+-"<<setw(8)<<eyn<<" / "<<setw(8)<<yd<<"+-"<<setw(8)<<eyd<<" ==> "<<setw(10)<<prob<<" +- "<<setw(10)<<eprob;
    

    //if(yd<10000) {cout<<endl; continue;} //limiting categories to be used
    cout<<" ======>>>> SELECTED "<<endl;


    chi2.setPoint( prob, eprob, p1, p2);
    np++;
    
    if(p1==p2 && yd>10000) {
      float d=4-4*prob;
      float x1=(2-sqrt(d))/2;
      float x2=(2+sqrt(d))/2;

      float de=4-4*(prob+eprob);
      float x1e=(2-sqrt(de))/2 - x1;
      float x2e=(2+sqrt(de))/2 - x2;

      cout<<"("<<itv->second[0]<<","<<itv->second[1]<<") fixing -> "<<(x1<1?x1:x2)<<" +-  "<<(x1e<1?x1e:x2e)<<endl; //<<eprob/2<<"  "
      fixProbs[p1][0]=(x1<1?x1:x2);
      fixProbs[p1][1]=(x1e<1?x1e:x2e);//eprob/2; //aproximate
    }

    // cout<<itv->first<<" ==> "<<bPt1<<"/"<<bEta1<<":"<<p1<<"   "<<bPt2<<"/"<<bEta2<<":"<<p2<<endl;
    //" :: "<<yn<<"  "<<yd<<"  "<<yn/yd<<"      "
    // <<itv->second[0]<<"  "<<itv->second[1]<<"  "
    // <<itv->second[2]<<"  "<<itv->second[3]<<endl;

  }
  cout<<" number of points in the fit: "<<np<<endl;
}

binSt setPointsFromDB(string file, Chi2& chi2, map<int,float[2]>& fixProbs) {

  ifstream categs(file.c_str(), ios::in);
  string line;

  vector<float> binsPt;
  vector<float> binsEta;
  
  map<string, float> yields;
  map<string, float> eyields;
  map<string, vector<float> > bins;

  //first read
  vector<float> bin;
  float y, ey;
  while(getline(categs, line)) 
    {
      istringstream iss(line);
      vector<string> tks;
      copy(istream_iterator<string>(iss),
	   istream_iterator<string>(),
	   back_inserter<vector<string> >(tks));

      bin = parseCateg(tks[0]);
      y= atof(tks[1].c_str() );
      ey= atof(tks[2].c_str() );

      bool flag=true;
      for(size_t i=0;i<binsPt.size();i++) {
	if(binsPt[i]==bin[0]) { flag=false; break;}
      }
      if(flag) binsPt.push_back(bin[0]);
	 
      flag=true;
      for(size_t i=0;i<binsEta.size();i++) {
	if(binsEta[i]==bin[1]) { flag=false; break;}
      }
      if(flag) binsEta.push_back(bin[1]);

      bins[tks[0]] = bin;
      yields[tks[0]] = y;
      eyields[tks[0]] = ey;
    }

  fillPoints( binsPt, binsEta, yields, eyields, bins, chi2, fixProbs);
  
  binSt binstruct(binsPt, binsEta);
  return binstruct;
 
}


binSt setPoints(map<string, vector<float> > vals, Chi2& chi2, map<int,float[2]>& fixProbs) {

  //first read the categories
  vector<float> binsPt;
  vector<float> binsEta;
  
  map<string, float> yields;
  map<string, float> eyields;
  map<string, vector<float> > bins;

  vector<float> bin;
  float y, ey;

  map<string, vector<float> >::const_iterator it;
  for(it=vals.begin();it!=vals.end();++it) {

    bin = parseCateg(it->first);
    y= it->second[0];
    ey= it->second[1];

    bool flag=true;
    for(size_t i=0;i<binsPt.size();i++) {
      if(binsPt[i]==bin[0]) { flag=false; break;}
    }
    if(flag) binsPt.push_back(bin[0]);
    
    flag=true;
    for(size_t i=0;i<binsEta.size();i++) {
      if(binsEta[i]==bin[1]) { flag=false; break;}
    }
    if(flag) binsEta.push_back(bin[1]);
    
    bins[it->first] = bin;
    yields[it->first] = y;
    eyields[it->first] = ey;
  }

  
  fillPoints( binsPt, binsEta, yields, eyields, bins, chi2, fixProbs);

  binSt binsstruct(binsPt, binsEta);
  return binsstruct;
}


int main(int argc, char* argv[]) {

  string file;
  string mcfile;
  bool isData=false;
  bool isRootFile=true;
  string singleCateg="";
  bool appendDb=false;
  string dbName="database.db";
  bool cutAndCount=false;
  
  char c;

  while ((c = getopt(argc, argv, "f:m:d:s:D:a:n:c:h")) != -1 ) {
    switch (c) {
      //case 'd': { file=optarg; break;}
    case 'f': { file=string(optarg); break;}
    case 'm': { mcfile=string(optarg); break;}
    case 'd': { isRootFile=bool(1-atoi(optarg)); break;}
    case 's': { singleCateg=string(optarg); break;}
    case 'D': { isData=bool(atoi(optarg)); break;}
    case 'a': { appendDb=bool(atoi(optarg)); break;}
    case 'n': { dbName=string(optarg); break;}
    case 'c': { cutAndCount=true; break;}
    case 'h': { 
      cout<<"configuration options:\n -f : file to read (root or ASCII) \n -m : MC file to read (root) \n -d proceed with a database reading instead of making fits (0 per default). \n -s <categ> perform a fit over a single Z category. \n -D run on data (0 per default). \n -a store the numbers into an existing database (1 per default). \n -n set the database file name for reading (database.db per default). \n -h help \n"<<endl;
      return 0; }
    default : { 
      cout<<"configuration options:\n -f : file to read (root or ASCII) \n -d proceed with a database reading instead of making fits (0 per default). \n -s <categ> perform a fit over a single Z category. \n -D run on data (0 per default). \n -a store the numbers into an existing database (1 per default). \n -n set the database file name for reading (database.db per default). \n -h help \n"<<endl;
      return 0; }
    }
  }

  if(singleCateg=="") appendDb=true;

  cout<<" ======= Configuration options ======== "<<endl;
  cout<<"\tfile     : "<<file<<endl;
  cout<<"\tfile MC  : "<<mcfile<<endl;
  cout<<"\trootfile : "<<isRootFile<<endl;
  cout<<"\tsingle   : "<<singleCateg<<endl;
  cout<<"\tis data  : "<<isData<<endl;
  cout<<"\t append  : "<<appendDb<<endl;
  cout<<"\t db name : "<<dbName<<endl;
  cout<<" ===================================== "<<endl<<endl;

  //==============================================
  gErrorIgnoreLevel = kWarning;


  Chi2 chi2;binSt bins;map<int,float[2]> fixProbs;
  if(!isRootFile) { //read the DB file
    bins=setPointsFromDB(file, chi2, fixProbs);
  } 
  else { //read root file 
    map<string, vector<float> > vals=doFits(file, mcfile, isData, appendDb,
					    dbName, singleCateg, cutAndCount);
    return 0;
    bins= setPoints(vals,chi2, fixProbs);
  }

  if(singleCateg!="") return 0;
  
  //perform the final fit ====================
  int nvars=bins.getNProb();

  ROOT::Fit::Fitter  fitter;
  ROOT::Math::Functor fcn(chi2,nvars);

  //bloody ROOT and lack of vector handling
  double* vars= new double[nvars];
  fitter.SetFCN(fcn,vars);

  // set step sizes and limits
  for (int i=0; i<nvars; ++i) {
    fitter.Config().ParSettings(i).SetStepSize(0.000000001);
    fitter.Config().ParSettings(i).SetLimits(0,2);
    fitter.Config().ParSettings(i).SetValue(0.00001);
  }

  for(map<int,float[2]>::const_iterator it=fixProbs.begin();it!=fixProbs.end();it++) {
    int ip=it->first;
    float val=it->second[0];
    float eval=it->second[1];
    fitter.Config().ParSettings(ip).SetValue(val);
    fitter.Config().ParSettings(ip).Fix();
  }

  bool ok = fitter.FitFCN();
  if (!ok) {
    cout<<" The final fit did not converged properly, please check your data and the read/written database "<<endl; 
    return 1;
  }

  fitter.CalculateMinosErrors();
  fitter.CalculateHessErrors();
  const ROOT::Fit::FitResult & result = fitter.Result();
  result.Print(std::cout);

  const double * parFit = result.GetParams();
  const double * parErrsTmp = result.GetErrors();
  vector<double> parErrs;
  for(int i=0;i<nvars;i++) {
    map<int,float[2]>::const_iterator it=fixProbs.find(i);
    if(it!=fixProbs.end()) parErrs.push_back(it->second[1]);
    else parErrs.push_back( parErrsTmp[i] );
  }

  cout<<"probabilities (i/ptbin/etabin) and fit result"<<endl;
  for(int i=0;i<nvars;i++) {
    cout<<i<<"    "<<i/bins.nBEta<<"   "<<i%bins.nBEta<<" ==> "<<parFit[i]<<" +- "<<parErrs[i]<<endl;
  }

  size_t p=file.find(".root");
  if(!isRootFile) p=file.find(".db");
  string tag=file.substr(0, p);

  string fname=tag+(isData?"_data":"_MC")+".root";
  TFile* f=new TFile(fname.c_str(),"recreate");
  f->cd();

  //adding the last bin boundaries, bloody root...
  double* binsPt=new double[ bins.nBPt+1 ];
  double* binsEta=new double[ bins.nBEta+1 ];

  for(size_t ib=0;ib<max( bins.nBPt, bins.nBEta );ib++) {
    if(ib<bins.nBPt) binsPt[ib]= bins.binsPt[ib];
    if(ib<bins.nBEta) binsEta[ib]= bins.binsEta[ib];
  }
  binsPt[bins.nBPt] = 300;
  binsEta[bins.nBEta] = 2.5;
  
  TH2F* h=new TH2F("chargeMisId","chargeMisId;p_{T}(e) [GeV];#eta(e) ",bins.nBPt,binsPt,bins.nBEta,binsEta);
  for(int i=0;i<nvars;i++) {
    h->SetBinContent((i/bins.nBEta)+1,(i%bins.nBEta)+1, parFit[i]);
    h->SetBinError((i/bins.nBEta)+1,(i%bins.nBEta)+1, parErrs[i]);
    map<int,float[2]>::const_iterator it=fixProbs.find(i);
    // if(it!=fixProbs.end()) 
    //   h->SetBinError((i/bins.nBEta)+1,(i%bins.nBEta)+1, it->second[1]);
  }
 
  f->Write();
  f->Close();

}
