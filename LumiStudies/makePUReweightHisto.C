#include <iostream>
#include <string>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>

using namespace std;

void makePUReweightHisto(string dfname="pileup800pb.root", string mcfname="MCPileup.root", string tag="800pb"){

  TFile* dfile=TFile::Open(dfname.c_str(),"read");
  TFile* mcfile=TFile::Open(mcfname.c_str(),"read");
  
  if(tag!="") tag="_"+tag;
  
  TH1F* dh=(TH1F*)dfile->Get("pileup");
  TH1F* mch=(TH1F*)mcfile->Get("pileup");

  // dh->Scale(1./ dh->Integral(0,100000) );
  // mch->Scale(1./ mch->Integral(0,100000) );

  TH1F* ratio=(TH1F*)dh->Clone("pileup");
  ratio->Reset("ICEM");
  ratio->SetName("puw");

  TH1F* dhNorm=(TH1F*)dh->Clone("pileup");
  TH1F* mchNorm=(TH1F*)mch->Clone("pileup");
  dhNorm->Reset("ICEM");
  mchNorm->Reset("ICEM");

  float intD=dh->Integral(0,100000);
  float intMC=mch->Integral(0,100000);

  for(int i=0;i<dh->GetNbinsX()+2;i++) {
    float xu=dh->GetBinContent(i)/intD;
    float exu=dh->GetBinError(i)/intD;
    

    float xd=mch->GetBinContent(i)/intMC;
    float exd=mch->GetBinError(i)/intMC;
    //cout<<" --> "<<i<<endl;
    dhNorm->SetBinContent(i,xu);
    dhNorm->SetBinError(i,exu);
    mchNorm->SetBinContent(i,xd);
    mchNorm->SetBinError(i,exd);

    //cout<<dhNorm->GetBinContent(i)<<"  "<<mchNorm->GetBinContent(i)<<endl;

    if(xd==0 || xu==0) continue;
    //cout<<xu/xd<<"  "<<exu/xu<<"  "<<exd/xd<<" ==> "<<sqrt( pow(exu/xu,2) + pow(exd/xd,2) )*xu/xd<<" --> "<<exd<<"  "<<xd<<endl;

    ratio->SetBinContent(i, xu/xd);
    ratio->SetBinError(i, sqrt( pow(exu/xu,2) + pow(exd/xd,2) )*xu/xd );
  }
  
  TFile* ofile=TFile::Open(("puWeights"+tag+".root").c_str(),"recreate");
  ofile->cd();
  ratio->Write();
  ofile->Write();
  ofile->Close();
  

  TCanvas* c=new TCanvas("c","c");
  //c->cd();
  dhNorm->SetLineWidth(2);
  dhNorm->SetLineColor(kBlack);
  
  mchNorm->SetLineWidth(2);
  mchNorm->SetLineColor(kRed+1);

  dhNorm->Draw();
  mchNorm->Draw("same");




}
