#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TChain.h"

using namespace std;


// void getBoundaries(TTree* tree) {

//   vector<float> etaBins({0,0.8,1.479,2.5});
//   vector<float> ptBins({10});

//   float ptUp=11;
//   while(ptUp<1000) {

//     int nEvts=tree->GetEntries("");


//   }



// }




void makeFlipPlots() {

  //TFile* file=TFile::Open("DYJetsToLL_M50_LO.root","read");
  TFile* filett=TFile::Open("TT_pow_ext4.root","read");
  TTree* tree=(TTree*)filett->Get("tree");

  tree->SetWeight(12500*832./169720600. ); 
  tree->AutoSave(); 

  TFile* filedy=TFile::Open("DYJetsToLL_M50_LO.root","read");
  TTree* treedy=(TTree*)filedy->Get("tree");

  treedy->SetWeight( 12500*6021./49619140); 
  treedy->AutoSave(); 

  // TChain* tree = new TChain("tree");
  // tree->Add("TT_pow_ext4.root");
  // tree->Add("DYJetsToLL_M50_LO.root"); 


  vector<string> tags({
      "OS",
	"OSflips",
	"OSfakes",
	"SS",
	"SSflips",
	"SSfakes",
	"SSfakesPred",
	"SSprompt"
	});

  map<string, TH1F*> histos;
  for(size_t ii=0;ii<tags.size();ii++) {
    histos[ tags[ii] ] = new TH1F( ("h"+tags[ii]).c_str(), ("h"+tags[ii]).c_str(), 80, 50, 130 );
  }

  map<string, string> sels;
  sels[ "OS" ] = "Sum$(LepGood_charge)==0 && aux==1 && nLepGood==2";
  sels[ "OSflips" ] = "Sum$(LepGood_charge)==0 && ( (LepGood_mcUCSXMatchId[0]==0&&LepGood_mcUCSXMatchId[1]==0) || (LepGood_mcUCSXMatchId[0]==0&&LepGood_mcUCSXMatchId[1]==0) ) && aux!=1 && nLepGood==2";
  sels[ "OSfakes" ] = "Sum$(LepGood_charge)==0 && ( (LepGood_mcUCSXMatchId[0]>1&&LepGood_mcUCSXMatchId[1]==0) || (LepGood_mcUCSXMatchId[0]==0&&LepGood_mcUCSXMatchId[1]>1) ) && aux!=1 && nLepGood==2";
  sels[ "SS" ] = "Sum$(LepGood_charge)!=0 && aux==1&& nLepGood==2";
  sels[ "SSflips" ] = "Sum$(LepGood_charge)!=0 && ( (LepGood_mcUCSXMatchId[0]==1&&LepGood_mcUCSXMatchId[1]==0) || (LepGood_mcUCSXMatchId[0]==0&&LepGood_mcUCSXMatchId[1]==1) ) && aux==1 && nLepGood==2";
  sels[ "SSfakes" ] = "Sum$(LepGood_charge)!=0 && ( (LepGood_mcUCSXMatchId[0]>1&&LepGood_mcUCSXMatchId[1]==0) || (LepGood_mcUCSXMatchId[0]==0&&LepGood_mcUCSXMatchId[1]>1) ) && aux==1 && nLepGood==2";
  sels[ "SSfakesPred" ] = "(Sum$(LepGood_charge)!=0 && ( (LepGood_mcUCSXMatchId[0]>1&&LepGood_mcUCSXMatchId[1]==0) || (LepGood_mcUCSXMatchId[0]==0&&LepGood_mcUCSXMatchId[1]>1) ) && w!=1 && nLepGood==2)*w"; //
  sels[ "SSprompt" ] = "Sum$(LepGood_charge)!=0 && aux==1&& (LepGood_mcUCSXMatchId[0]==0 && LepGood_mcUCSXMatchId[1]==0) && nLepGood==2";

  map<string, float> integs;

  for(size_t ii=0;ii<tags.size();ii++) {
    tree->Draw( ("mZ>>h"+tags[ii]).c_str(), (sels[ tags[ii] ]).c_str(), "goff" );
    cout<<tags[ii]<<" -> "<<histos[ tags[ii] ]->Integral(0,100000)<<endl;
    treedy->Draw( ("mZ>>+h"+tags[ii]).c_str(), (sels[ tags[ii] ]).c_str(), "goff" );

    cout<<tags[ii]<<" ==-> "<<histos[ tags[ii] ]->Integral(0,100000)<<endl;
    integs[ tags[ii] ]=histos[ tags[ii] ]->Integral(0,100000);

    histos[ tags[ii] ]->Rebin(2);
  }

  cout<<" All : "<<integs["SS"]<<" / "<<integs["OS"]<<" = "<<integs["SS"]/(float)integs["OS"]<<endl;
  cout<<" fakes only : "<<integs["SSfakes"]<<" ( "<<integs["SSfakesPred"]<<" ) / "<<integs["OSfakes"]<<" = "<<integs["SSfakes"]/(float)integs["OSfakes"]<<endl;
  cout<<" flips only : "<<integs["SSflips"]<<" / "<<integs["OSflips"]<<" = "<<integs["SSflips"]/(float)integs["OSflips"]<<endl;


  TCanvas* c=new TCanvas("c","c",1600,800);
  c->Divide(4,2);
  c->cd(1);
  histos[ "SS" ]->Draw();
  c->cd(2);
  histos[ "SSflips" ]->Draw();
  c->cd(3);
  histos[ "SSfakes" ]->Draw();
  c->cd(4);
  histos[ "SSfakesPred" ]->Draw();
  c->cd(5);
  histos[ "OS" ]->Draw();
  c->cd(6);
  histos[ "OSflips" ]->Draw();
  c->cd(7);
  histos[ "OSfakes" ]->Draw();
  // c->cd(8);
  // histos[ "SSprompt" ]->Draw();
}


void splitTree(string MVA="", TString sample="DoubleEG.root") {

  TFile* filedy=TFile::Open(sample,"read"); //MCSamples.root //DoubleEG.root
  TTree* treedy=(TTree*)filedy->Get("tree");

  string auxS="aux";
  string wS="w";
  string mz="mZ";
  if(MVA!="") {
    auxS="aux"+MVA;
    wS="w"+MVA;
    mz="mZ"+MVA;
  }

  treedy->SetBranchStatus("*",0);
  treedy->SetBranchStatus(mz.c_str(),1);
  treedy->SetBranchStatus(wS.c_str(),1);
  treedy->SetBranchStatus(auxS.c_str(),1);
  treedy->SetBranchStatus("sample",1);
  treedy->SetBranchStatus("isData",1);
  treedy->SetBranchStatus("nLepGood",1);
  treedy->SetBranchStatus("LepGood_charge",1);
  treedy->SetBranchStatus("LepGood_pt",1);
  treedy->SetBranchStatus("LepGood_eta",1);
  if(((string)sample).find("DoubleEG")==string::npos)
    treedy->SetBranchStatus("LepGood_mcUCSXMatchId",1);

  TFile* ofile=new TFile("FlipTrees_"+sample+MVA.c_str(),"recreate");
  

  vector<float> ptBins({10,20,50,100,200,300});
  vector<float> etaBins({0,0.8,1.479,2.5});

  // for(size_t ip1=0;ip1<ptBins.size()-1;ip1++) {
  //   for(size_t ie1=0;ie1<etaBins.size()-1;ie1++) {

  //     for(size_t ip2=ip1;ip2<ptBins.size()-1;ip2++) {
  // 	for(size_t ie2=0;ie2<etaBins.size()-1;ie2++) {

  for(size_t ie1=0;ie1<etaBins.size()-1;ie1++) {
    for(size_t ip1=0;ip1<ptBins.size()-1;ip1++) {
  
      for(size_t ie2=ie1;ie2<etaBins.size()-1;ie2++) {
	for(size_t ip2=((ie2>ie1)?0:ip1);ip2<ptBins.size()-1;ip2++) {


	  // if(ie2<=ie1 && ip2<ip1) continue;
	  // if(ip2<=ip1 && ie2<ie1) continue;


	  ostringstream pd, pu, ed, eu;
	  pd << ptBins[ip1];
	  pu << ptBins[ip1+1];
	  ed << etaBins[ie1];
	  eu << etaBins[ie1+1];

	  ostringstream pd2, pu2, ed2, eu2;
	  pd2 << ptBins[ip2];
	  pu2 << ptBins[ip2+1];
	  ed2 << etaBins[ie2];
	  eu2 << etaBins[ie2+1];

	  string isDataOK=((((string)sample).find("DoubleEG")==string::npos)?"(!isData && Sum$(LepGood_mcUCSXMatchId)==1) &&":"");
	  string isDataBAD=((((string)sample).find("DoubleEG")==string::npos)?"(!isData && Sum$(LepGood_mcUCSXMatchId)==0) &&":"");

	  string cutSS="(Sum$(LepGood_charge)!=0 && "+auxS+"==1 && nLepGood==2 && "+mz+"!=-1";
	  string cutOS="(Sum$(LepGood_charge)==0 && "+auxS+"==1 && nLepGood==2 && "+mz+"!=-1";

	  string cutSSFake=""+wS+"*(Sum$(LepGood_charge)!=0 && "+auxS+"!=1 && nLepGood==2 && "+mz+"!=-1";
	  string cutOSFake=""+wS+"*(Sum$(LepGood_charge)==0 && "+auxS+"!=1 && nLepGood==2 && "+mz+"!=-1";

	  string cutSSPrompt=""+wS+"*(Sum$(LepGood_charge)!=0 && "+isDataOK+" "+auxS+"!=1 && nLepGood==2 && "+mz+"!=-1";
	  string cutOSPrompt=""+wS+"*(Sum$(LepGood_charge)==0 && "+isDataBAD+" "+auxS+"!=1 && nLepGood==2 && "+mz+"!=-1";

	  string cutSSPromptZ=""+wS+"*(Sum$(LepGood_charge)!=0 && "+isDataOK+" "+auxS+"==1 && nLepGood==2 && sample==1 && "+mz+"!=-1";
	  string cutOSPromptZ=""+wS+"*(Sum$(LepGood_charge)==0 && "+isDataBAD+" "+auxS+"==1 && nLepGood==2 && sample==1 && "+mz+"!=-1";

	  string cutSSPromptTT=""+wS+"*(Sum$(LepGood_charge)!=0 && "+isDataOK+" "+auxS+"==1 && nLepGood==2 && sample==2 && "+mz+"!=-1";
	  string cutOSPromptTT=""+wS+"*(Sum$(LepGood_charge)==0 && "+isDataBAD+" "+auxS+"==1 && nLepGood==2 && sample==2 && "+mz+"!=-1";

	  string cstring=" && ( ( ( LepGood_pt[0]>"+pd.str()+" && LepGood_pt[0]<"+pu.str()+" ) ";
	  cstring += " && ( abs(LepGood_eta[0])>"+ed.str()+" && abs(LepGood_eta[0])<"+eu.str()+" ) ";
	  cstring += " && ( LepGood_pt[1]>"+pd2.str()+" && LepGood_pt[1]<"+pu2.str()+" ) ";
	  cstring += " && ( abs(LepGood_eta[1])>"+ed2.str()+" && abs(LepGood_eta[1])<"+eu2.str()+" ) ) ";

	  //counterpart
	  string cstringcc="|| ( ( LepGood_pt[1]>"+pd.str()+" && LepGood_pt[1]<"+pu.str()+" ) ";
	  cstringcc += " && ( abs(LepGood_eta[1])>"+ed.str()+" && abs(LepGood_eta[1])<"+eu.str()+" ) ";
	  cstringcc += " && ( LepGood_pt[0]>"+pd2.str()+" && LepGood_pt[0]<"+pu2.str()+" ) ";
	  cstringcc += " && ( abs(LepGood_eta[0])>"+ed2.str()+" && abs(LepGood_eta[0])<"+eu2.str()+" ) ) ) ";

	  //final strings
	  cutSS += cstring + cstringcc + ")";
	  cutOS += cstring + cstringcc + ")";
	  cutSSFake += cstring + cstringcc + ")";
	  cutOSFake += cstring + cstringcc + ")";
	  cutSSPrompt += cstring + cstringcc + ")";
	  cutOSPrompt += cstring + cstringcc + ")";

	  cutSSPromptZ += cstring + cstringcc + ")";
	  cutOSPromptZ += cstring + cstringcc + ")";

	  cutSSPromptTT += cstring + cstringcc + ")";
	  cutOSPromptTT += cstring + cstringcc + ")";
	  
	  string name="l1Pt_"+pd.str()+"_l1Eta_"+ed.str()+"_l2Pt_"+pd2.str()+"_l2Eta_"+ed2.str();

	  TTree* treeSS=treedy->CopyTree( cutSS.c_str() );   treeSS->SetName( (name+"_SS").c_str() );
	  TTree* treeOS=treedy->CopyTree( cutOS.c_str() );   treeOS->SetName( (name+"_OS").c_str() );
	  TTree* treeSSF=treedy->CopyTree( cutSSFake.c_str() );   treeSSF->SetName( (name+"_SS_Fake").c_str() );
	  TTree* treeOSF=treedy->CopyTree( cutOSFake.c_str() );   treeSSF->SetTitle( (name+"_SS_Fake").c_str() ); 
	  cout<<cutOSFake<<"   --->>>> "<<treedy->GetEntries(cutOSFake.c_str())<<"   "<<treeOSF->GetEntries()<<endl;
	  TTree* treeSSP=treedy->CopyTree( cutSSPrompt.c_str() );   treeSSP->SetName( (name+"_SS_Prompt").c_str() );
	  TTree* treeOSP=treedy->CopyTree( cutOSPrompt.c_str() );    treeOSP->SetName( (name+"_OS_Prompt").c_str() );
	  TTree* treeSSPZ=treedy->CopyTree( cutSSPromptZ.c_str() );   treeSSPZ->SetName( (name+"_SS_PromptZ").c_str() );
	  TTree* treeOSPZ=treedy->CopyTree( cutOSPromptZ.c_str() );   treeOSPZ->SetName( (name+"_OS_PromptZ").c_str() );
	  TTree* treeSSPTT=treedy->CopyTree( cutSSPromptTT.c_str() );   treeSSPTT->SetName( (name+"_SS_PromptTT").c_str() );
	  TTree* treeOSPTT=treedy->CopyTree( cutOSPromptTT.c_str() );   treeOSPTT->SetName( (name+"_OS_PromptTT").c_str() );
	  
	 

	
	
	  treeSS->SetTitle( (name+"_SS").c_str() );
	  treeOS->SetTitle( (name+"_OS").c_str() );
	  treeOSF->SetName( (name+"_OS_Fake").c_str() );
	  treeOSF->SetTitle( (name+"_OS_Fake").c_str() );

	  treeSSP->SetTitle( (name+"_SS_Prompt").c_str() );
	  treeOSP->SetTitle( (name+"_OS_prompt").c_str() );

	  treeSSPZ->SetTitle( (name+"_SS_PromptZ").c_str() );
	  treeOSPZ->SetTitle( (name+"_OS_promptZ").c_str() );

	  treeSSPTT->SetTitle( (name+"_SS_PromptTT").c_str() );
	  treeOSPTT->SetTitle( (name+"_OS_promptTT").c_str() );
	  
	  cout<<ptBins[ip1]<<"  "<<etaBins[ie1]
	      <<" // "<<ptBins[ip2]<<"  "<<etaBins[ie2]<<" :: "<<treeSS->GetEntries()<<" <> "<<treeOS->GetEntries()
	      <<" :: "<<treeSSF->GetEntries()<<" <> "<<treeOSF->GetEntries()
	      <<" :: "<<treeSSP->GetEntries()<<" <> "<<treeOSP->GetEntries()<<endl;


	}//eta2
      }//pt2
    }//eta1
  }//pt1

  
  ofile->Write();
  ofile->Close();
  
}


long factoriel(long n) {
   return n > 1?(n * factoriel(n-1)):1;
}


long sum(long n) {
   return n > 1?(n + sum(n-1)):1;
}


void test() {


  vector<float> ptBins({10,20,50,100,200,300});
  vector<float> etaBins({0,0.8,1.479,2.5});
  int N=0;
  for(size_t ie1=0;ie1<etaBins.size()-1;ie1++) {
    for(size_t ip1=0;ip1<ptBins.size()-1;ip1++) {
  
      for(size_t ie2=ie1;ie2<etaBins.size()-1;ie2++) {
	for(size_t ip2=((ie2>ie1)?0:ip1);ip2<ptBins.size()-1;ip2++) {

	  
	  //if(ie2<=ie1 && ip2<ip1) continue;
	  // if(ip2<=ip1 && ie2<ie1) continue;
	  //if(ip2<ip1 && ie2<ie1) continue;

	  cout<<setw(5)<<ptBins[ip1]<<","<<setw(6)<<etaBins[ie1]<<"   "<<setw(5)<<ptBins[ip2]<<","<<setw(6)<<etaBins[ie2]<<endl;
	  N++;
	}
      }
    }
  }
  cout<<" combinations : "<<N<<endl;
  cout<<" factorielle 15 --> "<<sum(15)<<endl;

}
