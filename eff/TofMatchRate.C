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


void TofMatchRate() {

	TChain *chain = new TChain("JetTree");
	chain->Add("/home/hep/caines/ly247/Scratch/pp12MBPico_151207/sum*.root");	

	Jet *t = new Jet(chain);


	double ptmax = 10;
	double ptmin = 0.15;
	int ptbins = 100; 

	double etamax = 1;
	double etamin = -1;
	int etabins = 100; 

	int Ntrkmax = 100;

	const int NoZdc = 10;
	int ZdcCut[NoZdc+1] = {0,1000,2000,3000,4000,5000,6000,7000,8000,10000,50000};
	TH1D *htpc[NoZdc];
	TH1D *htof[NoZdc];
	TH1D *htpcVseta[NoZdc];
	TH1D *htofVseta[NoZdc];
	for(int i = 0; i<NoZdc; i++) {
		htpc[i] = new TH1D(Form("htpc%d",i), Form("Number of TPC tracks vs pt for %d<zdc<%d",ZdcCut[i],ZdcCut[i+1]),ptbins,ptmin,ptmax);
		htof[i] = new TH1D(Form("htof%d",i), Form("Number of TOF tracks vs pt for %d<zdc<%d",ZdcCut[i],ZdcCut[i+1]),ptbins,ptmin,ptmax);
		htpc[i]->Sumw2();
		htof[i]->Sumw2();

		htpcVseta[i] = new TH1D(Form("htpcVseta%d",i), Form("Number of TPC tracks vs #eta for %d<zdc<%d",ZdcCut[i],ZdcCut[i+1]),etabins,etamin,etamax);
		htofVseta[i] = new TH1D(Form("htofVseta%d",i), Form("Number of TOF tracks vs #eta for %d<zdc<%d",ZdcCut[i],ZdcCut[i+1]),etabins,etamin,etamax);
		htpcVseta[i]->Sumw2();
		htofVseta[i]->Sumw2();
	}

	const int NoVz = 12;
	int VzCut[NoVz+1] = {-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30};
	TH1D *htpcVseta_vzdep[NoVz];
	TH1D *htofVseta_vzdep[NoVz];
	int zdcmax_vzdep = 3000;
	for(int i = 0; i<NoVz; i++) {
		htpcVseta_vzdep[i] = new TH1D(Form("htpcVseta_vzdep%d",i), Form("Number of TPC tracks vs #eta for %d<zdc<%d %d<vz<%d",0,zdcmax_vzdep,VzCut[i],VzCut[i+1]),etabins,etamin,etamax);			
		htofVseta_vzdep[i] = new TH1D(Form("htofVseta_vzdep%d",i), Form("Number of TOF tracks vs #eta for %d<zdc<%d %d<vz<%d",0,zdcmax_vzdep,VzCut[i],VzCut[i+1]),etabins,etamin,etamax);
		htpcVseta_vzdep[i]->Sumw2();
		htofVseta_vzdep[i]->Sumw2();
	}

	cout<<"Total Events " <<chain->GetEntries()<<endl;
	for(int ievt = 0 ; ievt<chain->GetEntries() ; ievt++) {
                t->GetEntry(ievt);
                if(ievt%100000==0) cout<<ievt<<endl;


		//if(fabs(t->fEventHeader_fPVz)>5) continue;
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


		int vz = 0;
		while(t->fEventHeader_fPVz>VzCut[vz]) vz++;
		if(vz>=NoVz) vz = NoVz-1;

		for(int j = 0; j<t->fPrimaryTracks_; j++) {

			// track cuts
			double ieta = geteta(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j],t->fPrimaryTracks_fPz[j]);
			double ipt = getpt(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j]);
			double iphi = getphi(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j]);
			double idca = fabs(t->fPrimaryTracks_fDCA[j]);
			

			if((t->fPrimaryTracks_fFlag[j]<=0)) continue;
			if(ipt<ptmin && ipt>ptmax ) continue;
			if(idca>1) continue;
			if(t->fPrimaryTracks_fNFittedHits[j]<25) continue;
			if(1.*t->fPrimaryTracks_fNFittedHits[j]/t->fPrimaryTracks_fNHitsPoss[j]<0.52) continue;
			if(1.*t->fPrimaryTracks_fNFittedHits[j]/t->fPrimaryTracks_fNHitsPoss[j]>1.02) continue;

			if(t->fEventHeader_fZdcCoincidenceRate<zdcmax_vzdep) {
				htpcVseta_vzdep[vz]->Fill(ieta);
				if(t->fPrimaryTracks_fTofMatchFlag[j]>0) htofVseta_vzdep[vz]->Fill(ieta);
			}

			if(fabs(t->fEventHeader_fPVz)>5) continue;
			htpcVseta[zdc]->Fill(ieta);
			if(t->fPrimaryTracks_fTofMatchFlag[j]>0) htofVseta[zdc]->Fill(ieta);

			if(fabs(ieta)>0.5) continue;			
			htpc[zdc]->Fill(ipt);
			if(t->fPrimaryTracks_fTofMatchFlag[j]>0) htof[zdc]->Fill(ipt);

		}
	}

	TH1D *hr[NoZdc];
	TH1D *hrVseta[NoZdc];
	for(int i = 0; i<NoZdc; i++) {
		hr[i] = (TH1D*)htof[i]->Clone(Form("hr%d",i));
		hr[i]->Divide(htpc[i]);
		hrVseta[i] = (TH1D*)htofVseta[i]->Clone(Form("hrVseta%d",i));
		hrVseta[i]->Divide(htpcVseta[i]);
	}
	TH1D *hrVseta_vzdep[NoVz];
	for(int i = 0; i<NoVz; i++) {
		hrVseta_vzdep[i] = (TH1D*)htofVseta_vzdep[i]->Clone(Form("hrVseta_vzdep%d",i));
		hrVseta_vzdep[i]->Divide(htpcVseta_vzdep[i]);
	}

	// as the total effect has already been taken into account, the normalization is 1 here
	for(int i = 0; i<NoZdc; i++) {
		double sumE = hrVseta[i]->Integral();
		if(sumE) {
			hrVseta[i]->Scale(hrVseta[i]->GetNbinsX()/sumE);			
		}
	}

	TFile *fout = new TFile("ppTOFeff.root","recreate");
	for(int i = 0; i<NoZdc; i++) {
		htpc[i]->Write();
		htof[i]->Write();
		hr[i]->Write();
		htpcVseta[i]->Write();
		htofVseta[i]->Write();
		hrVseta[i]->Write();
	}
	for(int i = 0; i<NoVz; i++) {
		htpcVseta_vzdep[i]->Write();
		htofVseta_vzdep[i]->Write();
		hrVseta_vzdep[i]->Write();
	}
	fout->Close();

}


