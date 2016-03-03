#define Jet_cxx
#include "/home/hep/caines/ly247/QA/eff/Jet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Jet::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Jet.C
//      Root > Jet t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}



double getphi(double px, double py) {

        double phi = ((px==0)?0:atan(py/px)) ;

        if(px<0&&py<0) phi-=TMath::Pi();
        if(px<0&&py>0) phi+=TMath::Pi();

        return phi;
}


double geteta(double px, double py, double pz) {

        double p = sqrt(px*px+py*py+pz*pz);

        if(fabs(p-pz)<1e-12) return 1000;

        double e = 0.5*log((p+pz)/(p-pz));

        return e;
}

double getpt(double px, double py) {
        return sqrt(px*px+py*py);
}


