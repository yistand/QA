//====================================================================================================
//
//	2016.01.06	Li YI
//	Tracks, hits dependence on ZDC to invest pile-up effect
//
//====================================================================================================

#include <iostream>

#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TMath.h"

#include "/home/hep/caines/ly247/QA/eff/Jet.C"

using namespace std;

void ZdcDependence() {


	TChain *chain = new TChain("JetTree");
	chain->Add("/home/hep/caines/ly247/Scratch/pp12MBPico_151207/sum*.root");	

	Jet *t = new Jet(chain);


	double ptmax = 10;
	double ptmin = 0.15;
	int ptbins = 100; 
	int Ntrkmax = 100;

	TProfile *hrefmult = new TProfile("hrefmult","refmult vs ZDC",1000,1000,50000);
	TProfile *hglobal = new TProfile("hglobal","global vs ZDC",1000,1000,50000);
	TProfile *hprimary = new TProfile("hprimary","primary vs ZDC",1000,1000,50000);
	TProfile *htowermatched = new TProfile("htowermatched","towermatched vs ZDC",1000,1000,50000);
	TProfile *htower = new TProfile("htower","tower vs ZDC",1000,1000,50000);
	TProfile *hprimarymatched = new TProfile("hprimarymatched","primarymatched vs ZDC",1000,1000,50000);
	TProfile *hemcpoints = new TProfile("hemcpoints","emcpoints vs ZDC",1000,1000,50000);
	TProfile *hbbc = new TProfile("hbbc","bbc vs ZDC",1000,1000,50000);
	TProfile *hgoodtrk = new TProfile("hgoodtrk","goodtrk vs ZDC",1000,1000,50000);
	TProfile *hgoodtrkbemcmatch = new TProfile("hgoodtrkbemcmatch","goodtrkbemcmatch vs ZDC",1000,1000,50000);
	TProfile *hgoodtrktofmatch = new TProfile("hgoodtrktofmatch","goodtrktofmatch vs ZDC",1000,1000,50000);

	hrefmult->Sumw2();
	hglobal->Sumw2();
	hprimary->Sumw2();
	htowermatched->Sumw2();
	htower->Sumw2();
	hprimarymatched->Sumw2();
	hemcpoints->Sumw2();
	hbbc->Sumw2();
	hgoodtrk->Sumw2();
	hgoodtrkbemcmatch->Sumw2();
	hgoodtrktofmatch->Sumw2();

	for(int ievt = 0 ; ievt<chain->GetEntries() ; ievt++) {
                t->GetEntry(ievt);
                if(ievt%100000==0) cout<<ievt<<endl;


		if(fabs(t->fEventHeader_fPVz)>5) continue;
		if(fabs(t->fEventHeader_fPVz - t->fEventHeader_fvpdVz)>3) continue;
		if(t->fEventHeader_fRank<0) continue;
		int triggered = 0;
		for(int itrg = 0; itrg<t->fEventHeader_fNOfTriggerIds; itrg++) {
			if(t->fEventHeader_fTriggerIdArray[itrg]==370011) triggered = 1;		// VPDMB_nobsmd
			//if(t->mTrigId[itrg]==370341) triggered = 1;		// TOFMult4
			//if(t->mTrigId[itrg]==370361) triggered = 1;		// TOFMult3*VPD
		}
		if(triggered==0) continue;	

		double zdc = t->fEventHeader_fZdcCoincidenceRate;
		
		hrefmult->Fill(zdc,t->fEventHeader_fRefMult);
		hglobal->Fill(zdc,t->fEventHeader_fNOfGlobalTracks);
		hprimary->Fill(zdc,t->fEventHeader_fNOfPrimaryTracks);
		htowermatched->Fill(zdc,t->fEventHeader_fNOfTowerTrackMatched);
		htower->Fill(zdc,t->fEventHeader_fNOfTowers);
		hprimarymatched->Fill(zdc,t->fEventHeader_fNOfMatchedTracks);
		hemcpoints->Fill(zdc,t->fEventHeader_fNOfEMCPoints);
		hbbc->Fill(zdc,t->fEventHeader_fBbcCoincidenceRate);

		double goodtrk = 0;
		double goodtrkbemcmatch = 0;
		double goodtrktofmatch = 0;

		for(int j = 0; j<t->fPrimaryTracks_; j++) {
			// track cuts
			double ieta = geteta(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j],t->fPrimaryTracks_fPz[j]);
			double ipt = getpt(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j]);
			double iphi = getphi(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j]);
			double idca = fabs(t->fPrimaryTracks_fDCA[j]);
			

			if((t->fPrimaryTracks_fFlag[j]<=0)) continue;
			if(fabs(ieta)>0.5) continue;			
			if(ipt<ptmin && ipt>ptmax ) continue;
			if(idca>1) continue;
			if(t->fPrimaryTracks_fNFittedHits[j]<25) continue;
			if(1.*t->fPrimaryTracks_fNFittedHits[j]/t->fPrimaryTracks_fNHitsPoss[j]<0.52) continue;
			if(1.*t->fPrimaryTracks_fNFittedHits[j]/t->fPrimaryTracks_fNHitsPoss[j]>1.02) continue;

			goodtrk++;
			if(t->fPrimaryTracks_fBemcMatchFlag[j]>0) goodtrkbemcmatch++;
			if(t->fPrimaryTracks_fTofMatchFlag[j]>0) goodtrktofmatch++;

		}
		hgoodtrk->Fill(zdc,goodtrk);
		hgoodtrkbemcmatch->Fill(zdc,goodtrkbemcmatch);
		hgoodtrktofmatch->Fill(zdc,goodtrktofmatch);


	}

	TFile *fout = new TFile("ppTracksVsZdc_0.root","recreate");
	hrefmult->Write();
	hglobal->Write();
	hprimary->Write();
	htowermatched->Write();
	htower->Write();
	hprimarymatched->Write();
	hemcpoints->Write();
	hbbc->Write();
	hgoodtrk->Write();
	hgoodtrkbemcmatch->Write();
	hgoodtrktofmatch->Write();

	fout->Close();


}
