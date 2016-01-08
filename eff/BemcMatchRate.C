//=============================================================================================
//
//		2015.12.16	Li Yi
//		Fraction of tracks matched to TOF / total tracks in TPC
//
//=============================================================================================
#include <iostream>

#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TMath.h"

#include "/home/hep/caines/ly247/QA/eff/Jet.C"


void BemcMatchRate() {

	TChain *chain = new TChain("JetTree");
	chain->Add("/home/hep/caines/ly247/Scratch/pp12MBPico_151207/sum*.root");	

	Jet *t = new Jet(chain);


	double ptmax = 10;
	double ptmin = 0.15;
	int ptbins = 100; 
	int Ntrkmax = 100;

	const int NoZdc = 10;
	int ZdcCut[NoZdc+1] = {0,1000,2000,3000,4000,5000,6000,7000,8000,10000,50000};
	TH1D *htpc[NoZdc];
	TH1D *hemc[NoZdc];
	for(int i = 0; i<NoZdc; i++) {
		htpc[i] = new TH1D(Form("htpc%d",i), Form("Number of TPC tracks vs pt for %d<zdc<%d",ZdcCut[i],ZdcCut[i+1]),ptbins,ptmin,ptmax);
		hemc[i] = new TH1D(Form("hemc%d",i), Form("Number of BEMC tracks vs pt for %d<zdc<%d",ZdcCut[i],ZdcCut[i+1]),ptbins,ptmin,ptmax);
		htpc[i]->Sumw2();
		hemc[i]->Sumw2();
	}

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

		int zdc = 0;
		while(t->fEventHeader_fZdcCoincidenceRate>ZdcCut[zdc]) zdc++;
		if(zdc>=NoZdc) zdc = NoZdc-1;

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

			htpc[zdc]->Fill(ipt);
			if(t->fPrimaryTracks_fBemcMatchFlag[j]>0) hemc[zdc]->Fill(ipt);

		}
	}

	TH1D *hr[NoZdc];
	for(int i = 0; i<NoZdc; i++) {
		hr[i] = (TH1D*)hemc[i]->Clone(Form("hr%d",i));
		hr[i]->Divide(htpc[i]);
	}

	TFile *fout = new TFile("ppBEMCeff.root","recreate");
	for(int i = 0; i<NoZdc; i++) {
		htpc[i]->Write();
		hemc[i]->Write();
		hr[i]->Write();
	}
	fout->Close();

}


