{
#include <vector>

using namespace std;
int j=0;
char tempname[200];
TH1D * temphist;
vector<TFile*> filevec;



for(int i=0; i <=16; i++){//loop over all the results_PhaseII4_B_14TEV_HT1_140PileUp_ in the Results directory
if(i<10)sprintf(tempname,"../Results/results_PhaseII4_B_14TEV_HT1_140PileUp_0%d.root",i);
else sprintf(tempname,"../Results/results_PhaseII4_B_14TEV_HT1_140PileUp_%d.root",i);
TFile *File = new TFile(tempname, "R"); 

sprintf(tempname,"allEvents/pt600/MET_pt600_allEvents");
temphist = (TH1D *)File->Get(tempname)->Clone();

if(i<10)sprintf(tempname,"results_PhaseII4_B_14TEV_HT1_140PileUp_0%d.root",i);
else sprintf(tempname,"results_PhaseII4_B_14TEV_HT1_140PileUp_%d.root",i);
if(temphist->GetEntries() > 0 ) cout << " Suspicious Files: " << tempname <<  endl;


}//end of loop over files 



































}
