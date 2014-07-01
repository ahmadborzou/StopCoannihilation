#include <cstdio>
using namespace std;

drawPlot(string cutname="MET", string histname="j1Pt"){
TFile * BJfile =new TFile("../PhaseI_BJ_14TEV_NoPileUp.root","R");
TFile * TTfile =new TFile("../PhaseI_TT_14TEV_NoPileUp.root","R");
TFile * Sigfile =new TFile("../PhaseI_t2cc350310_14TEV_NoPileUp_00.root","R");
THStack * tempstack;
THStack * finalstack = new THStack("finalstack","Final Plot");
TH1D * temphist;
char tempname[200];
TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);

Float_t legendX1 = .50;
   Float_t legendX2 = .70;
   Float_t legendY1 = .65;
   Float_t legendY2 = .85;
   TLegend* catLeg1 = new TLegend(legendX1,legendY1,legendX2,legendY2);
catLeg1->SetTextSize(0.032);
   catLeg1->SetTextFont(42);


///add Wlv to the finalstack
sprintf(tempname,"Wlv/%s/%s_%s_Wlv",cutname.c_str(),histname.c_str(),cutname.c_str());
tempstack = (THStack *)BJfile->Get(tempname)->Clone();
temphist = (TH1D *) tempstack->GetStack()->Last();
temphist->SetFillColor(2);
catLeg1->AddEntry(temphist,"Wlv","f");
finalstack->Add(temphist);

///add TTbar to the finalstack
sprintf(tempname,"TTbar/%s/%s_%s_TTbar",cutname.c_str(),histname.c_str(),cutname.c_str());
tempstack = (THStack *)TTfile->Get(tempname)->Clone();
temphist = (TH1D *) tempstack->GetStack()->Last();
temphist->SetFillColor(8);
catLeg1->AddEntry(temphist,"TTbar","f");
finalstack->Add(temphist);



///add Zvv to the finalstack
sprintf(tempname,"Zvv/%s/%s_%s_Zvv",cutname.c_str(),histname.c_str(),cutname.c_str());
tempstack = (THStack *)BJfile->Get(tempname)->Clone();
temphist = (TH1D *) tempstack->GetStack()->Last();
temphist->SetFillColor(4);
catLeg1->AddEntry(temphist,"Zvv","f");
finalstack->Add(temphist);

///Draw the background stack
finalstack->Draw();

//Draw the signal on the same canvas
sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
tempstack = (THStack *) Sigfile->Get(tempname)->Clone();
temphist = (TH1D *) tempstack->GetStack()->Last();
temphist->SetLineColor(1);
temphist->SetLineWidth(2);
temphist->Draw("SAME");
catLeg1->AddEntry(temphist,"Signal","l");





catLeg1->SetFillColor(kWhite);
   catLeg1->SetBorderSize(0);
   catLeg1->Draw();





}
