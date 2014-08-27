//There some points here that should be cared about
//First, cutname must be initialized the same as the one in main.cpp. 
//The same thing is true about type. Also,  
//there is a loop in the part that loads histograms into vec_map_map. 
//Number of loops should be the same as size of the vector of histograms
//in main.cpp, i.e. number of histograms in each branch. Currently,
//we have 4 histograms Weight,HT,MHT, NJets 

#include <cassert>
#include "TH1.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include "TSystem.h"
#include "THStack.h"
#include "TFile.h"
#include <vector>
#include <map>
#include "TDirectory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"

using namespace std;

////Define a class to calculate the significance
double zbi(double n_on, double mu_b_hat, double sigma_b){

  //double n_on     = 140.;                         // total events in signal region (S+B)
  //double mu_b_hat = 83.33;                        // mean num of BG events expected in sig. region
  //double sigma_b  = 8.333;                        // uncertainty of mu_b_hat

  double tau      = mu_b_hat / (sigma_b*sigma_b); // scale factor to corresp. Noff/Non              
  double n_off    = tau*mu_b_hat;
  double P_Bi     = TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1);
  double Z_Bi     = sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);

  cout  <<"  total events in signal region (S+B)               - n_on     " <<n_on      <<endl
        <<"  mean num of BG events expected in sig. region     - mu_b_hat " <<mu_b_hat  <<endl
        <<"  uncertainty of mu_b_hat                           - sigma_b  " <<sigma_b   <<endl
        <<"  scale factor to corresp. Noff/Non                 - tau      " <<tau       <<endl
        <<"  tau*mu_b_hat                                      - n_off    " <<n_off     <<endl
        <<"  TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1) - P_Bi     " <<P_Bi      <<endl
        <<"  sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi)           - Z_Bi     " <<Z_Bi      <<endl;

  return Z_Bi;
}


class mainClass{
double tempvalue;
vector<double> Sig_scalevec, BJ_scalevec,TT_scalevec ,B_scalevec,BJJ_scalevec,BB_scalevec,BBB_scalevec,H_scalevec,LL_scalevec,LLB_scalevec,TB_scalevec,TJ_scalevec;
char tempname[200];
char newhistname[200];
bool stackingswitch;
string yesno;
vector<double> Sig_xs_vec,BJ_xs_vec,TT_xs_vec,B_xs_vec,BJJ_xs_vec,BB_xs_vec,BBB_xs_vec,H_xs_vec,LL_xs_vec,LLB_xs_vec,TB_xs_vec,TJ_xs_vec;
map<int, string> cutname; 
map<int, string> sigtype,BJtype,TTtype,OtherBGtype;
map<int, string> histname;
map<int, vector<double> > yieldmap; //the first component must correspond to the components of the cutname map. For instance yieldmap[0]=0 would correspond to nocut while yieldmap[0]=1 to Asys. The other components represnt the # of events of Wlv, Zvv, TTbar, and Signal.
map<string, map<string , vector<TH1D> > > map_map; //here the string shows the type of events and map is between cutname and a vector of histograms like HT, MHT,...
vector< map<string, map<string , vector<TH1D> > >  > vec_map_map; //this is a vector of previous map. Each component of the vector corresponds to one HTbin. 
                                                                  //In BJ case there are 7 of them.
TFile *file;
vector<TFile *> Sig_inputfilevec, BJ_inputfilevec, TT_inputfilevec,B_inputfilevec,BJJ_inputfilevec,BB_inputfilevec,BBB_inputfilevec,H_inputfilevec,LL_inputfilevec,LLB_inputfilevec,TB_inputfilevec,TJ_inputfilevec;
TH1D *temphist;
//KKH TH1D temphist;
vector<TH1D > vtemphist;
vector<TH1D > temphistvec;
THStack * tempstack;
TDirectory *cdtoitt;
TDirectory *cdtoit;


public:
mainClass(int luminosity){//constructor
//Importnat
//make sure this initialization of the 
//maps is the same as that in main.cpp
 
    cutname[0]="nocut";
    cutname[1]="Asys";
    cutname[2]="MET200";
    cutname[3]="jetone";
    cutname[4]="jettwo";
    cutname[5]="3jet";
    cutname[6]="dphi";
    cutname[7]="nolep";
    cutname[8]="MET";
    cutname[9]="pt250";
    cutname[10]="pt300";
    cutname[11]="pt350";
    cutname[12]="pt400";
    cutname[13]="pt450";
    cutname[14]="pt500";
    cutname[15]="pt600";
    cutname[16]="pt700";
    cutname[17]="pt800";
    cutname[18]="pt900";
    cutname[19]="pt1000";
    cutname[20]="pt1100";
    cutname[21]="pt1200";
    cutname[22]="pt1300";
    cutname[23]="pt1400";
    cutname[24]="pt1500";
/*
    cutname[25]="METBtag2";
    cutname[26]="pt250Btag2";
    cutname[27]="pt300Btag2";
    cutname[28]="pt350Btag2";
    cutname[29]="pt400Btag2";
    cutname[30]="pt450Btag2";
    cutname[31]="pt500Btag2";
    cutname[32]="pt600Btag2";
    cutname[33]="pt700Btag2";
    cutname[34]="pt800Btag2";
    cutname[35]="pt900Btag2";
    cutname[36]="pt1000Btag2";
    cutname[37]="pt1100Btag2";
    cutname[38]="pt1200Btag2";
    cutname[39]="pt1300Btag2";
    cutname[40]="pt1400Btag2";
    cutname[41]="pt1500Btag2";*/


    sigtype[0]="allEvents";
    OtherBGtype[0]="allEvents"; 
    BJtype[0]="allEvents";
    BJtype[1]="W";
    BJtype[2]="Wlv";
    BJtype[3]="Wjj";
    BJtype[4]="Z";
    BJtype[5]="Zll";
    BJtype[6]="Zvv";
    BJtype[7]="Zjj";
    BJtype[8]="photon";
    BJtype[9]="H";
 
   TTtype[0]="allEvents";
   TTtype[1]="TTbar";
   TTtype[2]="TTSingLep";
   TTtype[3]="TTdiLep";
   TTtype[4]="TThadronic";


  //KH
  histname[0]="weight";
  histname[1]="METAsys";
  histname[2]="MET";
  histname[3]="NJet";
  histname[4]="j1Pt";
  histname[5]="Jet1Eta";
  histname[6]="Jet1Phi";
histname[7]="j2Pt";
histname[8]="Jet2Eta";
histname[9]="Jet2Phi";
histname[10]="j3Pt";
histname[11]="Jet3Eta";
histname[12]="Jet3Phi";
histname[13]="DelPhij1j2";
histname[14]="NLep";
histname[15]="NElec";
histname[16]="NMuon";
histname[17]="NTau";
/*histname[18]="NBtagLoose";
histname[19]="BtagLoose1Pt";
histname[20]="BtagLoose1Eta";
histname[21]="BtagLoose1Phi";
histname[22]="BtagLoose2Pt";
histname[23]="BtagLoose2Eta";
histname[24]="BtagLoose2Phi";*/


  ///end of initialization of the maps

yieldmap.clear();



//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section

  //build a vector of scale factors
  //first load the cross sections into a vector
//Sig_xs_vec.push_back(0.757); /// v1
//Sig_xs_vec.push_back(1.12); // v2
//Sig_xs_vec.push_back(1.15); // v3
//Sig_xs_vec.push_back(1.14); // M(Stop,LSP)=(450,410) and also M(Stop,LSP)=(450,440)
//Sig_xs_vec.push_back(2.18); // M(Stop,LSP)=(400,390) and also M(Stop,LSP)=(400,360)
//Sig_xs_vec.push_back(4.41); // M(Stop,LSP)=(350,340) and also M(Stop,LSP)=(350,310)
Sig_xs_vec.push_back(2.11); // v4
//Sig_xs_vec.push_back(0.009635); //FullExceptStopv4


  double Sig_numberofevents =0;//this will use GetSumOfWeights() 
  const int Sig_nHT = 1;   // Total number of HT bin samples
  const int nHist = 18; // Number of histograms in each TDirectory


  for(int i=1; i<=Sig_nHT ; i++){
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSP_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv2_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv3_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc450410_14TEV_140PileUp_00.root");  
//sprintf(tempname,"../Results/results_PhaseII4_t2cc450440_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc400390_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc400360_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc350340_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc350310_14TEV_140PileUp_00.root");
sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp.root");
//sprintf(tempname,"../Results/results_PhaseII4_FullExceptStopv4_14TEV_140PileUp.root");

  file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*Sig_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    Sig_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=Sig_nHT; i++){
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSP_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv2_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv3_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc450410_14TEV_140PileUp_00.root");  
//sprintf(tempname,"../Results/results_PhaseII4_t2cc450440_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc400390_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc400360_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc350340_14TEV_140PileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseII4_t2cc350310_14TEV_140PileUp_00.root");
sprintf(tempname,"../Results/results_PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp.root");
//sprintf(tempname,"../Results/results_PhaseII4_FullExceptStopv4_14TEV_140PileUp.root");

Sig_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

tempstack = new THStack("stack","Binned Sample Stack");
//sprintf(tempname,"PhaseII4_Stop_CharmLSP_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_Stop_CharmLSPv2_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_Stop_CharmLSPv3_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc450410_14TEV_140PileUp_00.root");  
//sprintf(tempname,"PhaseII4_t2cc450440_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc400390_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc400360_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc350340_14TEV_140PileUp_00.root");
//sprintf(tempname,"PhaseII4_t2cc350310_14TEV_140PileUp_00.root");
sprintf(tempname,"PhaseII4_Stop_CharmLSPv4_14TEV_140PileUp.root");
//sprintf(tempname,"PhaseII4_FullExceptStopv4_14TEV_140PileUp.root");

file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=sigtype.begin(); itt!=sigtype.end();itt++){        // loop over different event types
    cdtoitt = file->mkdir((itt->second).c_str());
    cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
      cdtoit =  cdtoitt->mkdir((it->second).c_str());
      cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<Sig_nHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) Sig_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(Sig_scalevec[i]);
if(histname[j]=="MET"){
Sig_numberofevents+=(double)temphist->GetSumOfWeights();
} 
temphist->SetFillColor(i+2);
tempstack->Add(temphist);


               }//end of loop over HTbins 1..7
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(Sig_numberofevents);
}

}
Sig_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
        tempstack->Write(tempname);
        delete tempstack;
        tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
 c+=1;   }//end of loop over cutnames
  }//end of loop over event types
file->Close();
///////////////////////////////////////////////////////////////End of signal
cout << "##\n##\n##\n " << endl;
cout << "Do you need the stacked root files or the cutflow is enough? [y/n?] " << endl;
cin >> yesno;
while(yesno!="y" && yesno!="n" && yesno!="Y" && yesno!="N"){
cout << "please type y or n: " << endl;
cin >> yesno;
}
if(yesno=="Y" || yesno=="y"){stackingswitch=true;}
else{stackingswitch=false;}

//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section//BJ Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  BJ_xs_vec.push_back(34409.92339);
  BJ_xs_vec.push_back(2642.85309);
  BJ_xs_vec.push_back(294.12311);
  BJ_xs_vec.push_back(25.95000);
  BJ_xs_vec.push_back(2.42111);
  BJ_xs_vec.push_back(0.22690);
  BJ_xs_vec.push_back(0.02767);

  double BJ_numberofevents =0;
  const int bjnHT = 7;   // Total number of HT bin samples


  for(int i=1; i<=bjnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_BJ_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*BJ_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    BJ_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=bjnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_BJ_14TEV_HT%d_140PileUp.root",i);
BJ_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_BJ_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=BJtype.begin(); itt!=BJtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<bjnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) BJ_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(BJ_scalevec[i]);
if(histname[j]=="MET"){
BJ_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..7
if(histname[j]=="MET"){
if(itt->second=="Wlv" || itt->second=="Zvv"){
yieldmap[c].push_back(BJ_numberofevents);
}

}
BJ_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ//EndOfBJ


//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section

  //build a vector of scale factors
  //first load the cross sections into a vector
  TT_xs_vec.push_back(530.89358);
  TT_xs_vec.push_back(42.55351);
  TT_xs_vec.push_back(4.48209);
  TT_xs_vec.push_back(0.52795);
  TT_xs_vec.push_back(0.05449);

  double TT_numberofevents =0;
  const int ttnHT = 5;   // Total number of HT bin samples


  for(int i=1; i<=ttnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_TT_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*TT_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    TT_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=ttnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_TT_14TEV_HT%d_140PileUp.root",i);
TT_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_TT_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=TTtype.begin(); itt!=TTtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<ttnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) TT_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(TT_scalevec[i]);
if(histname[j]=="MET"){
TT_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(TT_numberofevents);
}

}
TT_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;   }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
///////////////end of TTbar


//B Section//B Section//B Section//B Section//B Section//B Section//B Section//B Section//B Section//B Section//B Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  B_xs_vec.push_back(200944.68129);
 
  double B_numberofevents =0;
  const int bnHT = 1;   // Total number of HT bin samples

  for(int i=1; i<=bnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_B_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*B_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    B_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=bnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_B_14TEV_HT%d_140PileUp.root",i);
B_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_B_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<bnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) B_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(B_scalevec[i]);
if(histname[j]=="MET"){
B_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..7
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(B_numberofevents);
}

}
B_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB//EndOfB


//BJJ Section//BJJ Section//BJJ Section//BJJ Section//BJJ Section//BJJ Section//BJJ Section//BJJ Section//BJJ Section//BJJ Section//BJJ Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  BJJ_xs_vec.push_back(86.45604);
  BJJ_xs_vec.push_back(4.34869);
  BJJ_xs_vec.push_back(0.32465);
  BJJ_xs_vec.push_back(0.03032);
  BJJ_xs_vec.push_back(0.00313);

  double BJJ_numberofevents =0;
  const int bjjnHT = 5;   // Total number of HT bin samples


  for(int i=1; i<=bjjnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_BJJ_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*BJJ_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    BJJ_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=bjjnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_BJJ_14TEV_HT%d_140PileUp.root",i);
BJJ_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_BJJ_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<bjjnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) BJJ_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(BJJ_scalevec[i]);
if(histname[j]=="MET"){
BJJ_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(BJJ_numberofevents);
}
}
BJJ_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ//EndOfBJJ

//BB Section//BB Section//BB Section//BB Section//BB Section//BB Section//BB Section//BB Section//BB Section//BB Section//BB Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  BB_xs_vec.push_back(249.97710);
  BB_xs_vec.push_back(35.23062);
  BB_xs_vec.push_back(4.13743);
  BB_xs_vec.push_back(0.41702);
  BB_xs_vec.push_back(0.04770);

  double BB_numberofevents =0;
  const int bbnHT = 5;   // Total number of HT bin samples


  for(int i=1; i<=bbnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_BB_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*BB_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    BB_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=bbnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_BB_14TEV_HT%d_140PileUp.root",i);
BB_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_BB_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<bbnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) BB_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(BB_scalevec[i]);
if(histname[j]=="MET"){
BB_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(BB_numberofevents);
}
}
BB_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB//EndOfBB


//BBB Section//BBB Section//BBB Section//BBB Section//BBB Section//BBB Section//BBB Section//BBB Section//BBB Section//BBB Section//BBB Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  BBB_xs_vec.push_back(2.57304);
  BBB_xs_vec.push_back(0.14935);
  BBB_xs_vec.push_back(0.01274);

  double BBB_numberofevents =0;
  const int bbbnHT = 3;   // Total number of HT bin samples


  for(int i=1; i<=bbbnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_BBB_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*BBB_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    BBB_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=bbbnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_BBB_14TEV_HT%d_140PileUp.root",i);
BBB_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_BBB_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<bbbnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) BBB_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(BBB_scalevec[i]);
if(histname[j]=="MET"){
BBB_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(BBB_numberofevents);
}
}
BBB_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB//EndOfBBB


//H Section//H Section//H Section//H Section//H Section//H Section//H Section//H Section//H Section//H Section//H Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  H_xs_vec.push_back(21.55990);
  H_xs_vec.push_back(1.11282);
  H_xs_vec.push_back(0.09188);
  H_xs_vec.push_back(0.01009);

  double H_numberofevents =0;
  const int hnHT = 4;   // Total number of HT bin samples


  for(int i=1; i<=hnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_H_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*H_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    H_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=hnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_H_14TEV_HT%d_140PileUp.root",i);
H_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_H_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<hnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) H_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(H_scalevec[i]);
if(histname[j]=="MET"){
H_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(H_numberofevents);
}
}
H_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH//EndOfH


//LL Section//LL Section//LL Section//LL Section//LL Section//LL Section//LL Section//LL Section//LL Section//LL Section//LL Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  LL_xs_vec.push_back(1341.36923);
  LL_xs_vec.push_back(156.29534);
  LL_xs_vec.push_back(42.40132);
  LL_xs_vec.push_back(2.84373);
  LL_xs_vec.push_back(0.20914);
  LL_xs_vec.push_back(0.02891);

  double LL_numberofevents =0;
  const int llnHT = 6;   // Total number of HT bin samples


  for(int i=1; i<=llnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_LL_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*LL_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    LL_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=llnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_LL_14TEV_HT%d_140PileUp.root",i);
LL_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_LL_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<llnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) LL_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(LL_scalevec[i]);
if(histname[j]=="MET"){
LL_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(LL_numberofevents);
}
}
LL_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL//EndOfLL



//LLB Section//LLB Section//LLB Section//LLB Section//LLB Section//LLB Section//LLB Section//LLB Section//LLB Section//LLB Section//LLB Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  LLB_xs_vec.push_back(2.97380);
  LLB_xs_vec.push_back(0.22854);
  LLB_xs_vec.push_back(0.02080);

  double LLB_numberofevents =0;
  const int llbnHT = 3;   // Total number of HT bin samples


  for(int i=1; i<=llbnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_LLB_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*LLB_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    LLB_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=llbnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_LLB_14TEV_HT%d_140PileUp.root",i);
LLB_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_LLB_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<llbnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) LLB_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(LLB_scalevec[i]);
if(histname[j]=="MET"){
LLB_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(LLB_numberofevents);
}
}
LLB_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB//EndOfLLB



//TB Section//TB Section//TB Section//TB Section//TB Section//TB Section//TB Section//TB Section//TB Section//TB Section//TB Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  TB_xs_vec.push_back(63.88923);
  TB_xs_vec.push_back(7.12172);
  TB_xs_vec.push_back(0.98030);
  TB_xs_vec.push_back(0.08391);
  TB_xs_vec.push_back(0.00953);


  double TB_numberofevents =0;
  const int tbnHT = 5;   // Total number of HT bin samples


  for(int i=1; i<=tbnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_TB_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*TB_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    TB_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=tbnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_TB_14TEV_HT%d_140PileUp.root",i);
TB_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_TB_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<tbnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) TB_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(TB_scalevec[i]);
if(histname[j]=="MET"){
TB_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(TB_numberofevents);
}
}
TB_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB//EndOfTB




//TJ Section//TJ Section//TJ Section//TJ Section//TJ Section//TJ Section//TJ Section//TJ Section//TJ Section//TJ Section//TJ Section

 //build a vector of scale factors
  //first load the cross sections into a vector
  TJ_xs_vec.push_back(109.73602);
  TJ_xs_vec.push_back(5.99325);
  TJ_xs_vec.push_back(0.37680);
  TJ_xs_vec.push_back(0.03462);
  TJ_xs_vec.push_back(0.00312);


  double TJ_numberofevents =0;
  const int tjnHT = 5;   // Total number of HT bin samples


  for(int i=1; i<=tjnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseII4_TJ_14TEV_HT%d_140PileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*TJ_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    TJ_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=tjnHT; i++){
sprintf(tempname,"../Results/results_PhaseII4_TJ_14TEV_HT%d_140PileUp.root",i);
TJ_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
sprintf(tempname,"PhaseII4_TJ_14TEV_140PileUp.root");
if(stackingswitch==true)file = new TFile(tempname,"RECREATE");
 for(map<int , string >::iterator itt=OtherBGtype.begin(); itt!=OtherBGtype.end();itt++){        // loop over different event types
if(stackingswitch==true)cdtoitt = file->mkdir((itt->second).c_str());
if(stackingswitch==true)cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
if(stackingswitch==true)cdtoit =  cdtoitt->mkdir((it->second).c_str());
if(stackingswitch==true)cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<tjnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) TJ_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(TJ_scalevec[i]);
if(histname[j]=="MET"){
TJ_numberofevents+=(double)temphist->GetSumOfWeights();
}
temphist->SetFillColor(i+2);
if(stackingswitch==true)tempstack->Add(temphist);
               }//end of loop over HTbins 1..5
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(TJ_numberofevents);
}
}
TJ_numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
if(stackingswitch==true)tempstack->Write(tempname);
if(stackingswitch==true)delete tempstack;
if(stackingswitch==true)tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
if(stackingswitch==true)file->Close();
//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ//EndOfTJ


///write the output in a file
fstream ff;
ff.open("CutFlow.txt", std::fstream::out);
ff << " Cut Name,    " << "  Signal,      " << "  Wlv,      " << "  Zvv,     " << "  TTbar,      "<< "    Total BG,   " << " % Signal/Background,   "  <<  "    Significance, " << "  ZBI" << endl; 
double totalBG=0, delWlv=0, delZvv=0, delTT=0, delB=0, delBsquare=0;
for(int i=0; i<yieldmap.size(); i++){
totalBG=(double) (yieldmap[i].at(1)+yieldmap[i].at(2)+yieldmap[i].at(3));
delWlv= 0.08*yieldmap[i].at(1);///uncrtainty for Wlv is 8%
delZvv= 0.05*yieldmap[i].at(2);
delTT= 0.5*yieldmap[i].at(3);///uncrtainty for TTbar is 50%
delBsquare=pow(delWlv,2)+pow(delZvv,2)+pow(delTT,2);///delta_background = sqrt(delWlv^2+delZvv^2+delTT^2)
ff << "  " <<cutname[i]<<",     " << yieldmap[i].at(0) << ",     " << yieldmap[i].at(1) <<",     " << yieldmap[i].at(2) <<",     " <<yieldmap[i].at(3) << ",      "<< totalBG << ",      " << yieldmap[i].at(0)/totalBG*100  <<  ",       " << yieldmap[i].at(0)/sqrt(delBsquare+totalBG+yieldmap[i].at(0))<< ",    " << zbi((totalBG+yieldmap[i].at(0)),totalBG ,sqrt(delBsquare))  << endl;  
}
ff.close();


}//end of the constructor
};
int main(){
mainClass mainObj(3000000);
//mainClass mainObj(19700);
cout << " done :) " << endl;

}

