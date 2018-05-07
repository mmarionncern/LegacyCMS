
void 
makeFlipPlot(TString fname="Moriond_MC.root") {

  TGaxis::SetMaxDigits(3);
  gStyle->SetPaintTextFormat("4.1e");
  gStyle->SetOptStat(0);
  //  gStyle->SetPalette(51,0);
  
  TFile* file=new TFile(fname,"read");
  TH2F* chMIdProb_el=(TH2F*)file->Get("chargeMisId");


  TCanvas* cEl=new TCanvas("cEl","Electron Charge MisId",1000,700);
  cEl->cd();
  cEl->SetLogz(1);
  chMIdProb_el->SetTitle("charge misId (fit)");
  chMIdProb_el->GetXaxis()->SetTitle("p_{T}(e) [GeV] ");
  chMIdProb_el->GetXaxis()->SetLabelSize(0.04);
  chMIdProb_el->GetXaxis()->SetLabelOffset(-0.009);
  chMIdProb_el->GetYaxis()->SetLabelSize(0.05);
  chMIdProb_el->GetXaxis()->SetTitleOffset(0.90);
  chMIdProb_el->GetXaxis()->SetTitleSize(0.05);
  chMIdProb_el->GetYaxis()->SetTitleSize(0.06);
  chMIdProb_el->GetYaxis()->SetTitleOffset(0.80);
  chMIdProb_el->GetYaxis()->SetTitle("#eta  ");
  chMIdProb_el->GetZaxis()->SetRangeUser(0.000001,0.01);
  chMIdProb_el->SetContour(100);
  chMIdProb_el->Draw("colz texte");
  cEl->SetLogx(1);
  chMIdProb_el->GetXaxis()->SetMoreLogLabels();
  

  TLine* l1=new TLine(20,0 ,20,2.5);       l1->SetLineWidth(2);l1->SetLineStyle(7);l1->SetLineColor(kGray+2);
  TLine* l2=new TLine(50,0 ,50,2.5);       l2->SetLineWidth(2);l2->SetLineStyle(7);l2->SetLineColor(kGray+2);
  TLine* l3=new TLine(100,0 ,100,2.5);     l3->SetLineWidth(2);l3->SetLineStyle(7);l3->SetLineColor(kGray+2);
  TLine* l4=new TLine(200,0 ,200,2.5);     l4->SetLineWidth(2);l4->SetLineStyle(7);l4->SetLineColor(kGray+2);
  TLine* l6=new TLine(0,1.479 ,300,1.479); l6->SetLineWidth(2);l6->SetLineStyle(7);l6->SetLineColor(kGray+2);
  TLine* l7=new TLine(0,0.8 ,300,0.8);     l7->SetLineWidth(2);l7->SetLineStyle(7);l7->SetLineColor(kGray+2);

  l1->Draw();
  l2->Draw();
  l3->Draw();
  l4->Draw();
  l6->Draw();
  l7->Draw();

  
}
