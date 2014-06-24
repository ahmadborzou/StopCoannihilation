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

using namespace std;

class mainClass{
double tempvalue;
vector<double> Sig_scalevec, BJ_scalevec, TT_scalevec;
char tempname[200];
char newhistname[200];
vector<double> Sig_xs_vec,BJ_xs_vec,TT_xs_vec;
map<int, string> cutname; 
map<int, string> type;
map<int, string> histname;
map<int, vector<double> > yieldmap; //the first component must correspond to the components of the cutname map. For instance yieldmap[0]=0 would correspond to nocut while yieldmap[0]=1 to Asys. The other components represnt the # of events of Wlv, Zvv, TTbar, and Signal.
map<int, vector<double> > yield_map; //the first component must correspond to the components of the cutname map. For instance yieldmap[0]=0 would correspond to nocut while yieldmap[0]=1 to Asys. The other components represnt the # of events of Wlv, Zvv, TTbar, and Signal.

map<string, map<string , vector<TH1D> > > map_map; //here the string shows the type of events and map is between cutname and a vector of histograms like HT, MHT,...
vector< map<string, map<string , vector<TH1D> > >  > vec_map_map; //this is a vector of previous map. Each component of the vector corresponds to one HTbin. 
                                                                  //In BJ case there are 7 of them.
TFile *file;
vector<TFile *> Sig_inputfilevec, BJ_inputfilevec, TT_inputfilevec;
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


    type[0]="allEvents";
    type[1]="W";
    type[2]="Wlv";
    type[3]="Wjj";
    type[4]="Z";
    type[5]="Zll";
    type[6]="Zvv";
    type[7]="Zjj";
    type[8]="photon";
    type[9]="H";
   type[10]="TTbar";
   type[11]="TTSingLep";
   type[12]="TTdiLep";
   type[13]="TThadronic";



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
  ///end of initialization of the maps

int metcut = 0;
//for(int metcut=250; metcut<1500 ; metcut+=50){

yieldmap.clear();
yield_map.clear();



//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section//Signal Section

  //build a vector of scale factors
  //first load the cross sections into a vector
//Sig_xs_vec.push_back(0.757); /// v1
//Sig_xs_vec.push_back(1.12); // v2
//Sig_xs_vec.push_back(1.15); // v3
//Sig_xs_vec.push_back(1.14); // M(Stop,LSP)=(450,410) and also M(Stop,LSP)=(450,440)
//Sig_xs_vec.push_back(2.18); // M(Stop,LSP)=(400,390) and also M(Stop,LSP)=(400,360)
Sig_xs_vec.push_back(4.41); // M(Stop,LSP)=(350,340) and also M(Stop,LSP)=(350,310)



  double Sig_numberofevents =0;//this will use GetSumOfWeights() 
  double Sig_numofevents =0; //this will use Integral() to study tighter cuts
  const int Sig_nHT = 1;   // Total number of HT bin samples
  const int nHist = 18; // Number of histograms in each TDirectory


  for(int i=1; i<=Sig_nHT ; i++){
//sprintf(tempname,"../Results/results_PhaseI_Stop_CharmLSP_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_Stop_CharmLSPv2_14TEV_NoPileUp_02.root");
//sprintf(tempname,"../Results/results_PhaseI_Stop_CharmLSPv3_14TEV_NoPileUp_03.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc450410_14TEV_NoPileUp_00.root");  
//sprintf(tempname,"../Results/results_PhaseI_t2cc450440_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc400390_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc400360_14TEV_NoPileUp_00.root");
//sprintf(tempname,"../Results/results_PhaseI_t2cc350340_14TEV_NoPileUp_00.root");
sprintf(tempname,"../Results/results_PhaseI_t2cc350310_14TEV_NoPileUp_00.root");




  file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*Sig_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    Sig_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=Sig_nHT; i++){
sprintf(tempname,"../Results/results_PhaseI_Stop_CharmLSP_14TEV_NoPileUp_00.root");
Sig_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

tempstack = new THStack("stack","Binned Sample Stack");
//file = new TFile("Sig_stack.root","RECREATE");
 for(map<int , string >::iterator itt=type.begin(); itt!=type.end();itt++){        // loop over different event types
//    cdtoitt = file->mkdir((itt->second).c_str());
//    cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
//      cdtoit =  cdtoitt->mkdir((it->second).c_str());
//      cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<Sig_nHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) Sig_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(Sig_scalevec[i]);
if(histname[j]=="MET"){
Sig_numberofevents+=(double)temphist->GetSumOfWeights();
Sig_numofevents+=(double)temphist->Integral(temphist->GetXaxis()->FindBin(metcut),temphist->GetXaxis()->FindBin(5000));
} 
temphist->SetFillColor(i+2);
tempstack->Add(temphist);


               }//end of loop over HTbins 1..7
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(Sig_numberofevents);
yield_map[c].push_back(Sig_numofevents);
}
/*cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
cout << "type: " << (itt->second).c_str() << ",  cutname: " << (it->second).c_str() << endl;
cout << "Number of events:  " << Sig_numberofevents << endl;
cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
*/
}
Sig_numberofevents=0;
Sig_numofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
//        tempstack->Write(tempname);
        delete tempstack;
        tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
 c+=1;   }//end of loop over cutnames
  }//end of loop over event types
//file->Close();

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
  double BJ_numofevents =0;
  const int bjnHT = 7;   // Total number of HT bin samples


  for(int i=1; i<=bjnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*BJ_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    BJ_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=bjnHT; i++){
sprintf(tempname,"../Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp.root",i);
BJ_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

tempstack = new THStack("stack","Binned Sample Stack");
//file = new TFile("stack.root","RECREATE");
 for(map<int , string >::iterator itt=type.begin(); itt!=type.end();itt++){        // loop over different event types
 //   cdtoitt = file->mkdir((itt->second).c_str());
 //   cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
  //    cdtoit =  cdtoitt->mkdir((it->second).c_str());
  //    cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<bjnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) BJ_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(BJ_scalevec[i]);
if(histname[j]=="MET"){
BJ_numberofevents+=(double)temphist->GetSumOfWeights();
BJ_numofevents+=(double)temphist->Integral(temphist->GetXaxis()->FindBin(metcut),temphist->GetXaxis()->FindBin(5000));
}
temphist->SetFillColor(i+2);
tempstack->Add(temphist);
               }//end of loop over HTbins 1..7
if(histname[j]=="MET"){
if(itt->second=="Wlv" || itt->second=="Zvv"){
yieldmap[c].push_back(BJ_numberofevents);
yield_map[c].push_back(BJ_numofevents);
}
/*cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
cout << "type: " << (itt->second).c_str() << ",  cutname: " << (it->second).c_str() << endl;
cout << "Number of events:  " << BJ_numberofevents << endl;
cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
*/
}
BJ_numberofevents=0;
BJ_numofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
//        tempstack->Write(tempname);
        delete tempstack;
        tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;    }//end of loop over cutnames
  }//end of loop over event types
//file->Close();

//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section//TTbar Section

  //build a vector of scale factors
  //first load the cross sections into a vector
  TT_xs_vec.push_back(530.89358);
  TT_xs_vec.push_back(42.55351);
  TT_xs_vec.push_back(4.48209);
  TT_xs_vec.push_back(0.52795);
  TT_xs_vec.push_back(0.05449);

  double TT_numberofevents =0;
  double TT_numofevents =0;
  const int ttnHT = 5;   // Total number of HT bin samples


  for(int i=1; i<=ttnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseI_TT_14TEV_HT%d_NoPileUp_00.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*TT_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    TT_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=ttnHT; i++){
sprintf(tempname,"../Results/results_PhaseI_TT_14TEV_HT%d_NoPileUp_00.root",i);
TT_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

tempstack = new THStack("stack","Binned Sample Stack");
//file = new TFile("stack.root","RECREATE");
 for(map<int , string >::iterator itt=type.begin(); itt!=type.end();itt++){        // loop over different event types
  //  cdtoitt = file->mkdir((itt->second).c_str());
  //  cdtoitt->cd();
int c=0;
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
   //   cdtoit =  cdtoitt->mkdir((it->second).c_str());
   //   cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<ttnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) TT_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(TT_scalevec[i]);
if(histname[j]=="MET"){
TT_numberofevents+=(double)temphist->GetSumOfWeights();
TT_numofevents+=(double)temphist->Integral(temphist->GetXaxis()->FindBin(metcut),temphist->GetXaxis()->FindBin(5000));
}
temphist->SetFillColor(i+2);
          tempstack->Add(temphist);
               }//end of loop over HTbins 1..7
if(histname[j]=="MET"){
if(itt->second=="allEvents"){
yieldmap[c].push_back(TT_numberofevents);
yield_map[c].push_back(TT_numofevents);
}
/*cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
cout << "type: " << (itt->second).c_str() << ",  cutname: " << (it->second).c_str() << endl;
cout << "Number of events:  " << TT_numberofevents << endl;
cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
*/
}
TT_numberofevents=0;
TT_numofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
//        tempstack->Write(tempname);
        delete tempstack;
        tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
c+=1;   }//end of loop over cutnames
  }//end of loop over event types
//file->Close();

///write the output in a file
fstream ff;
ff.open("CutFlow.txt", std::fstream::out);
ff << " Cut Name,    " << "  Signal,      " << "  Wlv,      " << "  Zvv,     " << "  TTbar,      "<< "    Total BG,   " << " % Signal/Background,   "  <<  "    Significance " << endl; 
double totalBG=0, delWlv=0, delZvv=0, delTT=0, delB=0, delBsquare=0;
for(int i=0; i<yieldmap.size(); i++){
totalBG=(double) (yieldmap[i].at(1)+yieldmap[i].at(2)+yieldmap[i].at(3));
delWlv= 0.08*yieldmap[i].at(1);///uncrtainty for Wlv is 8%
delZvv= 0.05*yieldmap[i].at(2);
delTT= 0.5*yieldmap[i].at(3);///uncrtainty for TTbar is 50%
delBsquare=pow(delWlv,2)+pow(delZvv,2)+pow(delTT,2);///delta_background = sqrt(delWlv^2+delZvv^2+delTT^2)
ff << "  " <<cutname[i]<<",     " << yieldmap[i].at(0) << ",     " << yieldmap[i].at(1) <<",     " << yieldmap[i].at(2) <<",     " <<yieldmap[i].at(3) << ",      "<< totalBG << ",      " << yieldmap[i].at(0)/totalBG*100  <<  ",       " << yieldmap[i].at(0)/sqrt(delBsquare+totalBG+yieldmap[i].at(0))  <<endl;  
}
ff.close();
/*
fstream fff;
fff.open("TighterCutFlow.txt", std::fstream::out | std::fstream::app);
fff << " Cut Name    " << "  Signal        " << "  Wlv          " << "  Zvv           " << "  TTbar          "<< "    Total BG  " << " % Signal/Background"  <<  "          Significance " << endl; 
double totBG=0, delWlv=0, delZvv=0, delTT=0, delB=0, delBsquare=0;
for(int i=0; i<yield_map.size(); i++){
totBG=(double) (yield_map[i].at(1)+yield_map[i].at(2)+yield_map[i].at(3));
delWlv= 0.08*yield_map[i].at(1);///uncrtainty for Wlv is 8%
delZvv= 0.05*yield_map[i].at(2);
delTT= 0.5*yield_map[i].at(3);///uncrtainty for TTbar is 50%
delBsquare=pow(delWlv,2)+pow(delZvv,2)+pow(delTT,2);///delta_background = sqrt(delWlv^2+delZvv^2+delTT^2)
fff << "  " <<cutname[i]<<"       " << yield_map[i].at(0) << "       " << yield_map[i].at(1) <<"       " << yield_map[i].at(2) <<"       " <<yield_map[i].at(3) << "         "<< totBG << "           " << yield_map[i].at(0)/totBG*100  <<  "               " << yield_map[i].at(0)/sqrt(delBsquare+totBG+yield_map[i].at(0))  <<endl;  
}
fff.close();
*/
/*
///report high significances and S/B ratios
double totBG=0;
for(int l=0; l<yield_map.size(); l++){
totBG=(double) (yield_map[l].at(1)+yield_map[l].at(2)+yield_map[l].at(3));
if( (double) yield_map[l].at(0)/totBG*100 > 2.0 && (double) yield_map[l].at(0)/sqrt(totBG+yield_map[l].at(0)) >5.0){
cout << "metcut: " << metcut << ", cutname: " << cutname[l] << ", %S/B: " << yield_map[l].at(0)/totBG*100 << ", significance: " << yield_map[l].at(0)/sqrt(totBG+yield_map[l].at(0)) << endl; 
}
}
*/
//}//end of loop over metcut

}//end of the constructor
};
int main(){
mainClass mainObj(3000000);
//mainClass mainObj(19700);
cout << " done :) " << endl;

}

