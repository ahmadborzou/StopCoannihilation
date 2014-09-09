#include <cstdio>
using namespace std;

drawPlot(string cutname="RA2nocut", string histname="NJet", string scale="log"){

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////Some cosmetic work for official documents. 
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi_v2.C");

  writeExtraText = false;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"

  int iPeriod = 14;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 14= PU=140,14TeV 
  int iPos = 33;       // left-aligned
  
  int W = 600;
  int H = 600;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  //
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //
  TFile * BJfile  =new TFile("../PhaseII4_BJ_14TEV_140PileUp.root","R");
  TFile * TTfile  =new TFile("../PhaseII4_TT_14TEV_140PileUp.root","R");
  TFile * Bfile   =new TFile("../PhaseII4_B_14TEV_140PileUp.root","R");  
  TFile * BBfile  =new TFile("../PhaseII4_BB_14TEV_140PileUp.root","R");
  TFile * BBBfile =new TFile("../PhaseII4_BBB_14TEV_140PileUp.root","R");
  TFile * BJJfile =new TFile("../PhaseII4_BJJ_14TEV_140PileUp.root","R");
  TFile * Hfile   =new TFile("../PhaseII4_H_14TEV_140PileUp.root","R");
  TFile * LLfile  =new TFile("../PhaseII4_LL_14TEV_140PileUp.root","R");
  TFile * LLBfile =new TFile("../PhaseII4_LLB_14TEV_140PileUp.root","R");
  TFile * TBfile  =new TFile("../PhaseII4_TB_14TEV_140PileUp.root","R");
  TFile * TJfile  =new TFile("../PhaseII4_TJ_14TEV_140PileUp.root","R");
  TFile * TTBfile =new TFile("../PhaseII4_TTB_14TEV_140PileUp.root","R");
  //
  THStack * tempstack;
  THStack * finalstack = new THStack("finalstack","");
  TH1D * temphist, * VJhist, * Thist, * Otherhist, * TThist, * VVhist;
  char tempname[200];
  char xtitlename[200];
  char ytitlename[200];

  //
  // Define canvas
  // TCanvas * c1 = new TCanvas("c1", "c1", 800, 800); 
  // KH-added  
  char canvName[200];
  sprintf(canvName,"%s_%s",cutname.c_str(), histname.c_str());
  TCanvas* c1 = new TCanvas(canvName,canvName,10,10,W,H);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin( L/W );
  c1->SetRightMargin( R/W );
  c1->SetTopMargin( T/H );
  c1->SetBottomMargin( B/H );
  c1->SetTickx(0);
  c1->SetTicky(0);
  // KH-added

  //
  // Define TLegend
  Float_t legendX1 = .60;  //.50;
  Float_t legendX2 = .90;  //.70;
  Float_t legendY1 = .45;  //.65;
  Float_t legendY2 = .80;
  TLegend* catLeg1 = new TLegend(legendX1,legendY1,legendX2,legendY2);
  catLeg1->SetTextSize(0.032);
  catLeg1->SetTextFont(42);

  //
  // Define rebinning
  int rebin_MHT=2;
  int rebin_HT=2;

  /*
  ///add TTB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TTBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(7);
  catLeg1->AddEntry(temphist,"t#bar{t} + V","f");
  finalstack->Add(temphist);
  */

  /////////All Other Backgrounds
  //
  // Otherhist: TTB,BBB,H,LL,LLB
  // Thist:  TJ+TB
  // VVhist: BB
  // VJhist: B+BJ+BJJ
  // TThist: TT
  // 

  ///add BBB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BBBfile->Get(tempname)->Clone();
  Otherhist = (TH1D *) tempstack->GetStack()->Last();

  ///add TTB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TTBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Otherhist->Add(temphist);

  ///add H to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)Hfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Otherhist->Add(temphist);

  ///add LL to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)LLfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Otherhist->Add(temphist);

  ///add LLB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)LLBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Otherhist->Add(temphist);
  Otherhist->SetFillColor(4);  
  if(histname=="MHT" && rebin_MHT!=1) Otherhist->Rebin(rebin_MHT);
  if(histname=="HT"  && rebin_HT!=1)  Otherhist->Rebin(rebin_HT);
  finalstack->Add(Otherhist);

  /////endl of other backgrounds

  ///add BB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BBfile->Get(tempname)->Clone();
  VVhist = (TH1D *) tempstack->GetStack()->Last();
  VVhist->SetFillColor(6);
  if(histname=="MHT" && rebin_MHT!=1) VVhist->Rebin(rebin_MHT);
  if(histname=="HT"  && rebin_HT!=1)  VVhist->Rebin(rebin_HT);
  finalstack->Add(VVhist);

  ///add TB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TBfile->Get(tempname)->Clone();
  Thist = (TH1D *) tempstack->GetStack()->Last();

  ///add TJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  Thist->Add(temphist);
  Thist->SetFillColor(9);
  if(histname=="MHT" && rebin_MHT!=1) Thist->Rebin(rebin_MHT);
  if(histname=="HT"  && rebin_HT!=1)  Thist->Rebin(rebin_HT);
  finalstack->Add(Thist);

  ///add B to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)Bfile->Get(tempname)->Clone();
  VJhist = (TH1D *) tempstack->GetStack()->Last();

  ///add BJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  VJhist->Add(temphist);

  ///add BJJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  VJhist->Add(temphist);
  VJhist->SetFillColor(46);
  if(histname=="MHT" && rebin_MHT!=1) VJhist->Rebin(rebin_MHT);
  if(histname=="HT"  && rebin_HT!=1)  VJhist->Rebin(rebin_HT);
  finalstack->Add(VJhist);

  ///add TTbar to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TTfile->Get(tempname)->Clone();
  TThist = (TH1D *) tempstack->GetStack()->Last();
  TThist->SetFillColor(8);
  if(histname=="MHT" && rebin_MHT!=1) TThist->Rebin(rebin_MHT);
  if(histname=="HT"  && rebin_HT!=1)  TThist->Rebin(rebin_HT);
  finalstack->Add(TThist);

  ///
  catLeg1->AddEntry(TThist,"t#bar{t}","f");
  catLeg1->AddEntry(VJhist,"V + jets","f");
  catLeg1->AddEntry(Thist,"Single Top","f");
  catLeg1->AddEntry(VVhist,"VV","f");
  catLeg1->AddEntry(Otherhist,"Other BG","f");
  //

  ///Draw the background stack
  finalstack->Draw();
  //  finalstack->SetMinimum(200.);
  //  finalstack->GetHistogram()->GetXaxis()->SetTitle("p_{T}(jet1) [GeV]");

  //
  if(histname=="NBtagLoose"){
    sprintf(xtitlename,"Number of b-tags");
    sprintf(ytitlename,"Events");
    finalstack->GetXaxis()->SetLimits(0., 10.);
    if (scale=="log"){
      gPad->SetLogy();
      finalstack->SetMaximum(1000.);    
      finalstack->SetMinimum(0.1);    
    }
  }
  //
  else if(histname=="NJet"){
    //gPad->SetLogy();
    sprintf(xtitlename,"Number of jets");
    sprintf(ytitlename,"Events");
    finalstack->GetXaxis()->SetLimits(0., 20.);
    finalstack->SetMaximum(7.);
    if (cutname=="RA2noleptoncutht2500mht1300") finalstack->SetMaximum(100.);
  }
  //
  else if(histname=="MHT"){
    //gPad->SetLogy();
    sprintf(xtitlename,"#slash{H}_{T} (GeV)");
    sprintf(ytitlename,"Events / 100 GeV");
    finalstack->GetXaxis()->SetLimits(500., 3000.);
    if (scale=="log"){
      gPad->SetLogy();
      finalstack->SetMaximum(1000.);    
      finalstack->SetMinimum(0.3);    
    }
  }
  //
  else if(histname=="HT"){
    sprintf(xtitlename,"H_{T} (GeV)");
    sprintf(ytitlename,"Events / 200 GeV");
    finalstack->GetXaxis()->SetLimits(1000., 5000.);
    finalstack->GetXaxis()->SetNdivisions(505);
    finalstack->SetMaximum(15.);
    if (scale=="log"){
      gPad->SetLogy();
      finalstack->SetMaximum(1000.);    
      finalstack->SetMinimum(0.3);    
    }
  }
  //
  else{
    sprintf(xtitlename,"%s",histname.c_str());
    sprintf(ytitlename,"Events");
  }

  finalstack->GetHistogram()->GetXaxis()->SetTitle(xtitlename);
  finalstack->GetHistogram()->GetYaxis()->SetTitle(ytitlename);
  if(histname=="MHT"){
  }
  else if(histname=="HT"){
  }
  else{
  }
  //  c1->Modified();
 
  //------------------------------
  // Signal
  TFile * Sigfile =new TFile("../PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp.root","R");

  //Draw the signal on the same canvas
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());

  tempstack = (THStack *) Sigfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(4);
  temphist->SetLineWidth(2);
  //temphist->SetLineStyle(2);
  temphist->SetFillStyle(0);
  if(histname=="MHT" && rebin_MHT!=1) temphist->Rebin(rebin_MHT);
  if(histname=="HT"  && rebin_HT!=1)  temphist->Rebin(rebin_HT);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"STOC","l");

  gPad->RedrawAxis();

  //
  catLeg1->SetFillColor(kWhite);
  catLeg1->SetBorderSize(0);
  catLeg1->Draw();


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //For official documents
  // writing the lumi information and the CMS "logo"
  //
  {
    CMS_lumi_v2( c1, iPeriod, iPos );
  }

  c1->Update();
  c1->RedrawAxis();
  c1->GetFrame()->Draw();
  sprintf(tempname,"%s_%s.pdf",cutname.c_str(), histname.c_str());
  c1->Print(tempname,".pdf");
  sprintf(tempname,"%s_%s.png",cutname.c_str(), histname.c_str());
  //  c1->Print(tempname,".png");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
