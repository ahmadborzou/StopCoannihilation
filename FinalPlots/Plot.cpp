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
map<int, string> BGtype;
map<int, string> type;
map<int, string> histname;
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


//We only need 3 types of backgrounds 1. Wlv 2. Zvv 3. TTbar
BGtype[0]="Wlv";BGtype[1]="Zvv";BGtype[2]="TTbar";

  const int bjnHT = 7;   // Total number of HT bin samples
  const int ttnHT = 5;

 //build a vector of scale factors
  //first load the cross sections into a vector
  BJ_xs_vec.push_back(34409.92339);
  BJ_xs_vec.push_back(2642.85309);
  BJ_xs_vec.push_back(294.12311);
  BJ_xs_vec.push_back(25.95000);
  BJ_xs_vec.push_back(2.42111);
  BJ_xs_vec.push_back(0.22690);
  BJ_xs_vec.push_back(0.02767);

  //build a vector of scale factors
  //first load the cross sections into a vector
  TT_xs_vec.push_back(530.89358);
  TT_xs_vec.push_back(42.55351);
  TT_xs_vec.push_back(4.48209);
  TT_xs_vec.push_back(0.52795);
  TT_xs_vec.push_back(0.05449);


///Compute the Scale factors and load them to a vector
  for(int i=1; i<=bjnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*BJ_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    BJ_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;

///Open the input files
for(int i=1; i<=bjnHT; i++){
sprintf(tempname,"../Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp.root",i);
BJ_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

///Compute the Scale factors and load them to a vector
  for(int i=1; i<=ttnHT ; i++){
    sprintf(tempname,"../Results/results_PhaseI_TT_14TEV_HT%d_NoPileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*TT_xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    TT_scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;

///Open the input files
for(int i=1; i<=ttnHT; i++){
sprintf(tempname,"../Results/results_PhaseI_TT_14TEV_HT%d_NoPileUp.root",i);
TT_inputfilevec.push_back(TFile::Open(tempname,"R"));
}

tempstack = new THStack("stack","Binned Sample Stack");
file = new TFile("plot.root","RECREATE");
for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames

for(int j=0; j<histname.size(); j++){                                        // loop over different histograms

for(map<int , string >::iterator itt=BGtype.begin(); itt!=BGtype.end();itt++){        // loop over different event types

if(itt->second=="Wlv" || itt->second=="Zvv"){
for(int i=0; i<bjnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) BJ_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(BJ_scalevec[i]);
//temphist->SetFillColor(((int)(itt->first))+2);
tempstack->Add(temphist);
}//end of loop over bj ht bins
}else{
for(int i=0; i<ttnHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) TT_inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(TT_scalevec[i]);
//temphist->SetFillColor(((int)(itt->first))+2);
tempstack->Add(temphist);
}//end of loop over ttbar ht bins
}
}//end of loop over event types
sprintf(tempname,"%s_%s",(it->second).c_str(),histname[j].c_str());
tempstack->Write(tempname);
delete tempstack;
tempstack = new THStack("stack","Binned Sample Stack");

}//end of loop over histograms

}//end of loop over cutnames

file->Close();


}//end of the constructor
};
int main(){
mainClass mainObj(3000000);
cout << " done :) " << endl;

}

