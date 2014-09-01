#include <cstdio>
using namespace std;

drawPlot(string cutname="nocut", string histname="NJet"){

  //
  TFile * BJfile =new TFile("../PhaseII4_BJ_14TEV_140PileUp.root","R");
  TFile * TTfile =new TFile("../PhaseII4_TT_14TEV_140PileUp.root","R");
  TFile * Bfile =new TFile("../PhaseII4_B_14TEV_140PileUp.root","R");  
  TFile * BBfile =new TFile("../PhaseII4_BB_14TEV_140PileUp.root","R");
  TFile * BBBfile =new TFile("../PhaseII4_BBB_14TEV_140PileUp.root","R");
  TFile * BJJfile =new TFile("../PhaseII4_BJJ_14TEV_140PileUp.root","R");
  TFile * Hfile =new TFile("../PhaseII4_H_14TEV_140PileUp.root","R");
  TFile * LLfile =new TFile("../PhaseII4_LL_14TEV_140PileUp.root","R");
  TFile * LLBfile =new TFile("../PhaseII4_LLB_14TEV_140PileUp.root","R");
  TFile * TBfile =new TFile("../PhaseII4_TB_14TEV_140PileUp.root","R");
  TFile * TJfile =new TFile("../PhaseII4_TJ_14TEV_140PileUp.root","R");
//  TFile * TTBfile =new TFile("../PhaseII4_TTB_14TEV_140PileUp.root","R");
//
  THStack * tempstack;
  THStack * finalstack = new THStack("finalstack","Final Plot");
  TH1D * temphist;
  char tempname[200];
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
  gPad->SetLogy();

  Float_t legendX1 = .60;  //.50;
  Float_t legendX2 = .90;  //.70;
  Float_t legendY1 = .45;  //.65;
  Float_t legendY2 = .85;
  TLegend* catLeg1 = new TLegend(legendX1,legendY1,legendX2,legendY2);
  catLeg1->SetTextSize(0.032);
  catLeg1->SetTextFont(42);

  ///add Wlv to the finalstack
  sprintf(tempname,"Wlv/%s/%s_%s_Wlv",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(2);
  catLeg1->AddEntry(temphist,"W(lv)","f");
  //catLeg1->AddEntry(temphist,"Wlv","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;

  ///add TTbar to the finalstack
  sprintf(tempname,"TTbar/%s/%s_%s_TTbar",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TTfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(3);
  catLeg1->AddEntry(temphist,"t#bar{t}","f");
  //catLeg1->AddEntry(temphist,"TTbar","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;

  ///add Zvv to the finalstack
  sprintf(tempname,"Zvv/%s/%s_%s_Zvv",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(4);
  catLeg1->AddEntry(temphist,"Z(#nu#bar{#nu})","f");
  //catLeg1->AddEntry(temphist,"Zvv","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;

  ///add B to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)Bfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(5);
  catLeg1->AddEntry(temphist,"B","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;


  ///add BB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(6);
  catLeg1->AddEntry(temphist,"BB","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;


  ///add BBB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BBBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(7);
  catLeg1->AddEntry(temphist,"BBB","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;


  ///add BJJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)BJJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(40);
  catLeg1->AddEntry(temphist,"BJJ","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;

  ///add H to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)Hfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(28);
  catLeg1->AddEntry(temphist,"H","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;


  ///add LL to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)LLfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(12);
  catLeg1->AddEntry(temphist,"LL","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;


  ///add LLB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)LLBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(41);
  catLeg1->AddEntry(temphist,"LLB","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;


  ///add TB to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TBfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(8);
  catLeg1->AddEntry(temphist,"TB","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;


  ///add TJ to the finalstack
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
  tempstack = (THStack *)TJfile->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetFillColor(1);
  catLeg1->AddEntry(temphist,"TJ","f");
  finalstack->Add(temphist);
  std::cout << temphist->GetSumOfWeights() << std::endl;





  ///Draw the background stack
  finalstack->SetMinimum(200.);
  finalstack->SetMaximum(1000000.);
  finalstack->Draw();
//  finalstack->GetHistogram()->GetXaxis()->SetTitle("p_{T}(jet1) [GeV]");
sprintf(tempname,"%s",histname.c_str());
  finalstack->GetHistogram()->GetXaxis()->SetTitle(tempname);
  finalstack->GetHistogram()->GetYaxis()->SetTitle("Number of Events / 100 GeV");
//  finalstack->GetXaxis()->SetLimits(0., 20.);
//  c1->Modified();

  //------------------------------
  // Signal

/*TFile * Sigfile1 =new TFile("../PhaseII4_Stop_CharmLSP_14TEV_140PileUp_00.root","R");
  TFile * Sigfile2 =new TFile("../PhaseII4_Stop_CharmLSPv2_14TEV_140PileUp_00.root","R");
  TFile * Sigfile3 =new TFile("../PhaseII4_Stop_CharmLSPv3_14TEV_140PileUp_00.root","R");

  TFile * Sigfile4 =new TFile("../PhaseII4_t2cc450410_14TEV_140PileUp_00.root","R");
  TFile * Sigfile5 =new TFile("../PhaseII4_t2cc450440_14TEV_140PileUp_00.root","R");
  TFile * Sigfile6 =new TFile("../PhaseII4_t2cc400360_14TEV_140PileUp_00.root","R");
  TFile * Sigfile7 =new TFile("../PhaseII4_t2cc400390_14TEV_140PileUp_00.root","R");
  TFile * Sigfile8 =new TFile("../PhaseII4_t2cc350310_14TEV_140PileUp_00.root","R");
  TFile * Sigfile9 =new TFile("../PhaseII4_t2cc350340_14TEV_140PileUp_00.root","R");

  //----------
  //Draw the signal on the same canvas
  sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());

  
  tempstack = (THStack *) Sigfile1->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(4);
  temphist->SetLineWidth(2);
  //temphist->SetLineStyle(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"T2cc(400,360)","l");

  tempstack = (THStack *) Sigfile2->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(4);
  temphist->SetLineWidth(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"T2cc(400,390)","l");
  

  tempstack = (THStack *) Sigfile3->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(4);
  temphist->SetLineWidth(2);
  //temphist->SetLineStyle(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"STOC (449,410)","l");
  
  //-----
  tempstack = (THStack *) Sigfile4->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(6);
  temphist->SetLineWidth(2);
  temphist->SetLineStyle(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"T2cc(450,410)","l");

  tempstack = (THStack *) Sigfile5->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(6);
  temphist->SetLineWidth(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"T2cc(450,440)","l");

  
  tempstack = (THStack *) Sigfile6->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(4);
  temphist->SetLineWidth(2);
  temphist->SetLineStyle(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"T2cc(400,360)","l");

  tempstack = (THStack *) Sigfile7->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(4);
  temphist->SetLineWidth(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"T2cc(400,390)","l");

  tempstack = (THStack *) Sigfile8->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(2);
  temphist->SetLineWidth(2);
  temphist->SetLineStyle(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"T2cc(350,310)","l");

  tempstack = (THStack *) Sigfile9->Get(tempname)->Clone();
  temphist = (TH1D *) tempstack->GetStack()->Last();
  temphist->SetLineColor(2);
  temphist->SetLineWidth(2);
  temphist->SetFillStyle(0);
  temphist->Draw("SAME");
  std::cout << temphist->GetSumOfWeights() << std::endl;
  catLeg1->AddEntry(temphist,"T2cc(350,340)","l");
  */

  gPad->RedrawAxis();

  //
  catLeg1->SetFillColor(kWhite);
  catLeg1->SetBorderSize(0);
  catLeg1->Draw();

  sprintf(tempname,"%s_%s.pdf",cutname.c_str(), histname.c_str());
  c1->SaveAs(tempname);

}
