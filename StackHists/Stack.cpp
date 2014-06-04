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

using namespace std;

class mainClass{
double tempvalue;
vector<double> scalevec;
char tempname[200];
char newhistname[200];
vector<double> xs_vec;
map<int, string> cutname; 
map<int, string> type;
map<int, string> histname;
map<string, map<string , vector<TH1D> > > map_map; //here the string shows the type of events and map is between cutname and a vector of histograms like HT, MHT,...
vector< map<string, map<string , vector<TH1D> > >  > vec_map_map; //this is a vector of previous map. Each component of the vector corresponds to one HTbin. 
                                                                  //In BJ case there are 7 of them.
TFile *file;
vector<TFile *> inputfilevec;
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
//    type[9]="H";
 //   type[10]="TTbar";
 //   type[11]="TTSingLep";
 //   type[12]="TTdiLep";
 //   type[13]="TThadronic";



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

  //build a vector of scale factors
  //first load the cross sections into a vector
  xs_vec.push_back(34409.92339);
  xs_vec.push_back(2642.85309);
  xs_vec.push_back(294.12311);
  xs_vec.push_back(25.95000);
  xs_vec.push_back(2.42111);
  xs_vec.push_back(0.22690);
  xs_vec.push_back(0.02767);

  double numberofevents =0;
  const int nHT = 7;   // Total number of HT bin samples
  const int nHist = 18; // Number of histograms in each TDirectory


  for(int i=1; i<=nHT ; i++){
    sprintf(tempname,"../Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp.root",i);
    file = new TFile(tempname, "R");
    sprintf(tempname,"allEvents/nocut/MET_nocut_allEvents");
    tempvalue = (luminosity*xs_vec[i-1])/((* (TH1D* ) file->Get(tempname)).GetEntries());
    scalevec.push_back(tempvalue);
  }//end of loop over HTbins 
  std::cout << "normalization scale factor determination done" << std::endl;
for(int i=1; i<=nHT; i++){
sprintf(tempname,"../Results/results_PhaseI_BJ_14TEV_HT%d_NoPileUp.root",i);
inputfilevec.push_back(TFile::Open(tempname,"R"));
}

tempstack = new THStack("stack","Binned Sample Stack");
file = new TFile("stack.root","RECREATE");
 for(map<int , string >::iterator itt=type.begin(); itt!=type.end();itt++){        // loop over different event types
    cdtoitt = file->mkdir((itt->second).c_str());
    cdtoitt->cd();
    for(map<int , string >::iterator it=cutname.begin(); it!=cutname.end();it++){   // loop over different cutnames
      cdtoit =  cdtoitt->mkdir((it->second).c_str());
      cdtoit->cd();
      for(int j=0; j<histname.size(); j++){                                        // loop over different histograms
        for(int i=0; i<nHT ; i++){                                                  // loop over different HT bins
sprintf(tempname,"%s/%s/%s_%s_%s",(itt->second).c_str(),(it->second).c_str(),(histname[j]).c_str(),(it->second).c_str(),(itt->second).c_str());
temphist = (TH1D *) inputfilevec.at(i)->Get(tempname)->Clone();
temphist->Scale(scalevec[i]);
if(histname[j]=="MET"){numberofevents+=(double)temphist->GetSumOfWeights();} //all the histograms in one directory have the same number of events
//if(histname[j]=="MET"){cout << " temphist->GetSumOfWeights() " << temphist->GetSumOfWeights() << endl;}
/*if(i==0){
cout << "" << endl;
cout << "type: " << (itt->second).c_str() << ",  cutname: " << (it->second).c_str()<< ", histname: " << histname[j].c_str() << ", bin#: " << i << endl;
cout << "temphist->GetEntries():  " << temphist->GetEntries() << endl;
cout << "temphist->GetSumOfWeights(): " << temphist->GetSumOfWeights() << endl;
cout << " ===============================================================" << endl;
}
*/
temphist->SetFillColor(i+2);
          tempstack->Add(temphist);
               }//end of loop over HTbins 1..7
if(histname[j]=="MET"){
cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
cout << "type: " << (itt->second).c_str() << ",  cutname: " << (it->second).c_str() << endl;
cout << "Number of events:  " << numberofevents << endl;
cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}
numberofevents=0;
        sprintf(tempname,"%s_%s_%s",histname[j].c_str(),(it->second).c_str(),(itt->second).c_str());
        tempstack->Write(tempname);
        delete tempstack;
        tempstack = new THStack("stack","Binned Sample Stack");
      }//end of loop over histograms
    }//end of loop over cutnames
  }//end of loop over event types
file->Close();


}//end of the constructor
};

int main(){
//mainClass mainObj(3000000);
mainClass mainObj(19700);
cout << " done :) " << endl;

}

