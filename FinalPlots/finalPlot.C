#include <vector>
#include <cstdio>

using namespace std;


finalPlot(int luminosity=3000000,string cutname="pt250", string histname="j1Pt"){

vector<double> Sig_xs_vec, Sig_scalevec;
char tempname[200];
double tempvalue;
THStack * tempstack;
TH1D *temphist;
vector<TFile *> Sig_inputfilevec;


// Create canvas
TCanvas *c1 = new TCanvas("c1","c1",600,800);
gStyle->SetOptStat(0);
c1->SetLogy();


  //build a vector of scale factors
  //first load the cross sections into a vector
Sig_xs_vec.push_back(0.757); /// v1
//Sig_xs_vec.push_back(1.12); // v2
//Sig_xs_vec.push_back(1.15); // v3
//Sig_xs_vec.push_back(1.14); // M(Stop,LSP)=(450,410) and also M(Stop,LSP)=(450,440)
//Sig_xs_vec.push_back(2.18); // M(Stop,LSP)=(400,390) and also M(Stop,LSP)=(400,360)
//Sig_xs_vec.push_back(4.41); // M(Stop,LSP)=(350,340) and also M(Stop,LSP)=(350,310)

///Compute the Scale factors and load them to a vector
sprintf(tempname,"../Results/results_PhaseI_Stop_CharmLSP_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_Stop_CharmLSPv2_14TEV_NoPileUp_02.root");
//sprintf(tempname,"../Results/results_PhaseI_Stop_CharmLSPv3_14TEV_NoPileUp_03.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc450410_14TEV_NoPileUp_00.root");  
//sprintf(tempname,"../Results/results_PhaseI_t2cc450440_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc400390_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc400360_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc350340_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc350310_14TEV_NoPileUp_00.root");
  file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*Sig_xs_vec[0])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    Sig_scalevec.push_back(tempvalue);
  std::cout << "normalization scale factor determination done" << std::endl;

///Open the input files
sprintf(tempname,"../Results/results_PhaseI_Stop_CharmLSP_14TEV_NoPileUp_00.root");
Sig_inputfilevec.push_back(TFile::Open(tempname,"R"));


sprintf(tempname,"allEvents/%s/%s_%s_allEvents",cutname.c_str(),histname.c_str(),cutname.c_str());
temphist = (TH1D *) Sig_inputfilevec.at(0)->Get(tempname)->Clone();
temphist->Scale(Sig_scalevec[0]);
temphist->SetLineColor(kBlack);
temphist->Draw();







///Get a the background plot
TFile *file=new TFile("plot.root","R");
sprintf(tempname,"%s_%s",cutname.c_str(),histname.c_str());
tempstack = (THStack *)file->Get(tempname)->Clone();
tempstack->Draw("SAME");

















}
