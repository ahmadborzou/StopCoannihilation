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


///add Wlv to the finalstack
sprintf(tempname,"Wlv/%s/%s_%s_Wlv",cutname.c_str(),histname.c_str(),cutname.c_str());
tempstack = (THStack *)BJfile->Get(tempname)->Clone();
temphist = (TH1D *) tempstack->GetStack()->Last();
temphist->SetFillColor(2);
finalstack->Add(temphist);

///add TTbar to the finalstack
sprintf(tempname,"TTbar/%s/%s_%s_TTbar",cutname.c_str(),histname.c_str(),cutname.c_str());
tempstack = (THStack *)TTfile->Get(tempname)->Clone();
temphist = (TH1D *) tempstack->GetStack()->Last();
temphist->SetFillColor(8);
finalstack->Add(temphist);



///add Zvv to the finalstack
sprintf(tempname,"Zvv/%s/%s_%s_Zvv",cutname.c_str(),histname.c_str(),cutname.c_str());
tempstack = (THStack *)BJfile->Get(tempname)->Clone();
temphist = (TH1D *) tempstack->GetStack()->Last();
temphist->SetFillColor(4);
finalstack->Add(temphist);

///Draw the background stack
finalstack->Draw();

//Draw the signal on the same canvas
sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
tempstack = (THStack *) Sigfile->Get(tempname)->Clone();
temphist = (TH1D *) tempstack->GetStack()->Last();
temphist->SetLineColor(6);
temphist->Draw("SAME");



}
