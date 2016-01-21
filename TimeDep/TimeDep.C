//====================================================================================================
//
//	2016.01.08	Li Yi
//	Variables Vs runid
//
//====================================================================================================


#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <iomanip>	// std::setprecision


#include "TFile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TChain.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TClonesArray.h"

#include "/home/fas/caines/ly247/Software/PicoCode/eventStructuredAu/TStarJetPicoEvent.h"
#include "/home/fas/caines/ly247/Software/PicoCode/eventStructuredAu/TStarJetPicoEventCuts.h"
#include "/home/fas/caines/ly247/Software/PicoCode/eventStructuredAu/TStarJetPicoTriggerInfo.h"
#include "/home/fas/caines/ly247/Software/PicoCode/eventStructuredAu/TStarJetPicoTower.h"
#include "/home/fas/caines/ly247/Software/PicoCode/eventStructuredAu/TStarJetPicoPrimaryTrack.h"

using namespace std;

int convertphi(double phi, int offset=0) {		// convert phi into section. offset is used to set east & west section.
// West section is 1 to 12. offset should be 0
// East section is 13 to 24. offset should be 12
// https://drupal.star.bnl.gov/STAR/files/starnotes/csn0229a.pdf
	int sec = 0;
	double phip = TMath::Pi()*5./12.-phi;

	sec = TMath::Ceil(phip/(TMath::Pi()/6.));
	if(sec<=0) sec+=12;
	if(sec>12) sec-=12;
	sec+=offset;
	return sec;
}

int convertsetwest(int sec) {			// input 1-12, output 1-12 no change
	return sec;
}

int convertseteast(int sec) {			// input 1-12, output 13-24
// https://drupal.star.bnl.gov/STAR/files/starnotes/csn0229a.pdf
	int eastsec = 24 - sec;
	if(eastsec==12) eastsec = 24;
	return eastsec;
}

void convertrunindex(TProfile *oldp, TProfile* newp)	{ 		// runid -> runindex, use hp to figure out runindex

	vector<int> runindex;
	for(int oi = 0; oi<oldp->GetNbinsX(); oi++) {
		if(oldp->GetBinContent(oi+1)>0) {
			runindex.push_back(oi+1);
			//cout<<"+ "<<oi+1<<endl;
		}
	}

// copy TProfile needs to be done in the ways:
// - sum of bin y values 
// - sum of bin y*y values 
// - sum of bin weight = number of bin entries if profile is filled with weights =1 
// - sum of bin weight ^2 if the profile is filled with weights different than 1 
// all been copied for each bin
// SetBinContent & SetBinError as TH1D method not working for TProfile properly
	for(unsigned int ir = 0; ir<runindex.size(); ir++) {
		int newbin=ir+1;
		int oldbin=runindex.at(ir);
		//cout<<oldbin<<" --> "<<newbin<<endl;
		(*newp)[newbin]=(*oldp)[oldbin];		// copy y
		(*newp->GetSumw2())[newbin] = (*oldp->GetSumw2())[oldbin];		// copy y*y
		newp->SetBinEntries(newbin,oldp->GetBinEntries(oldbin));		// copy entries
		if( oldp->GetBinSumw2()->fN > oldbin) {					// copy (if needed) bin sum of weight square
			newp->Sumw2();
			(*newp->GetBinSumw2())[newbin] = (*oldp->GetBinSumw2())[oldbin];
		}
	}
	//cout<<"Done runid --> runindex"<<endl;
}


void convertrunindex(TProfile *hp, TH1D* h)	{ 		// runid -> runindex, use hp to figure out runindex
	vector<int> runindex;
	for(int i = 0; i<hp->GetNbinsX(); i++) {
		if(hp->GetBinContent(i+1)>0) {
			runindex.push_back(i+1);
		}
	}

	for(unsigned int ir = 0; ir<runindex.size(); ir++) {
		h->SetBinContent(ir+1,hp->GetBinContent(runindex.at(ir)));
		h->SetBinError(ir+1,hp->GetBinError(runindex.at(ir)));
	}
}


void convertrunindex(TH1D *hp, TH1D* h)	{ 		// runid -> runindex, use hp to figure out runindex
	vector<int> runindex;
	for(int i = 0; i<hp->GetNbinsX(); i++) {
		if(hp->GetBinContent(i+1)>0) {
			runindex.push_back(i+1);
		}
	}

	for(unsigned int ir = 0; ir<runindex.size(); ir++) {
		h->SetBinContent(ir+1,hp->GetBinContent(runindex.at(ir)));
		h->SetBinError(ir+1,hp->GetBinError(runindex.at(ir)));
	}
}


void convertrunindex(const vector<int>&  runindex, TProfile *hp, TH1D* h)	{ 		// runid -> runindex
	for(unsigned int ir = 0; ir<runindex.size(); ir++) {
		h->SetBinContent(ir+1,hp->GetBinContent(runindex.at(ir)));
		h->SetBinError(ir+1,hp->GetBinError(runindex.at(ir)));
	}
}

void convertrunindex(const vector<int>&  runindex, TProfile *oldp, TProfile* newp)	{ 		// runid -> runindex
// copy TProfile needs to be done in the ways:
// - sum of bin y values 
// - sum of bin y*y values 
// - sum of bin weight = number of bin entries if profile is filled with weights =1 
// - sum of bin weight ^2 if the profile is filled with weights different than 1 
// all been copied for each bin
// SetBinContent & SetBinError as TH1D method not working for TProfile properly
	for(unsigned int ir = 0; ir<runindex.size(); ir++) {
		int newbin=ir+1;
		int oldbin=runindex.at(ir);
		(*newp)[newbin]=(*oldp)[oldbin];		// copy y
		(*newp->GetSumw2())[newbin] = (*oldp->GetSumw2())[oldbin];		// copy y*y
		newp->SetBinEntries(newbin,oldp->GetBinEntries(oldbin));		// copy entries
		if( oldp->GetBinSumw2()->fN > oldbin) {					// copy (if needed) bin sum of weight square
			newp->Sumw2();
			(*newp->GetBinSumw2())[newbin] = (*oldp->GetBinSumw2())[oldbin];
		}
	}
}

bool AveSig(TH1D *h, double &ave, double &sigma) {		// average and sigma 
	ave = 0;
	sigma = 0;
	if(!h || h->GetNbinsX()==0) {
		cout<<"NULL input histogram in AveSig()"<<endl;
		return false;
	}
	if(h->GetNbinsX()==1) {
		ave = h->GetBinContent(1);
		sigma = 0;
		return false;
	}
	for(int i = 0; i<h->GetNbinsX(); i++) {
		ave+=h->GetBinContent(i+1);
	}
	ave/=h->GetNbinsX();

	for(int i = 0; i<h->GetNbinsX(); i++) {
		sigma+=pow((h->GetBinContent(i+1)-ave),2);
	}
	sigma/=h->GetNbinsX()-1;
	sigma=sqrt(sigma);
		
	return true;
}	


bool AveSig(TProfile *h, double &ave, double &sigma) {		// average and sigma 
	ave = 0;
	sigma = 0;
	if(!h || h->GetNbinsX()==0) {
		cout<<"NULL input histogram in AveSig()"<<endl;
		return false;
	}
	if(h->GetNbinsX()==1) {
		ave = h->GetBinContent(1);
		sigma = 0;
		return false;
	}
	for(int i = 0; i<h->GetNbinsX(); i++) {
		ave+=h->GetBinContent(i+1);
	}
	ave/=h->GetNbinsX();

	for(int i = 0; i<h->GetNbinsX(); i++) {
		sigma+=pow((h->GetBinContent(i+1)-ave),2);
	}
	sigma/=h->GetNbinsX()-1;
	sigma=sqrt(sigma);
		
	return true;
}	

void BadRun(TH1D *h, TH1D *hmap, vector<int> & badlist, const char *tag_trig = "NPE25", TString outfiletag="", const char *dir = "$HOME/Scratch/mapBEMCauau11Pico/", double sigmacut = 3, int CutOnOneSide = 0, double absolutecut = 0, double inputave = 0) {		
 // abnormal (>sigmacut-sigma) runid for variable in h (vs run index) 
 // If CutOnOneSide == 0:  symmetric cut on each side 
 // If CutOnOneSide == 1:  cut off the one exceed upper limit  
 // If CutOnOneSide == -1: cut off the one lower than low limit 
 // If absolutecut != 0 : use absolutecut as mean+/-absolutecut

	double ave = 0, sig = 0;
	AveSig(h,ave,sig);
	badlist.clear();

	if( fabs(inputave) > 1e-6) {
		ave = inputave;
	}
	absolutecut = fabs(absolutecut);
	if( fabs(absolutecut) < 1e-6) {		// if absolutecut==0, use relative cut
		absolutecut = sig*sigmacut;
	}

	cout<<"cut "<<ave<<"+/-"<<absolutecut<<endl;
	for(int i = 0; i<h->GetNbinsX(); i++) {
		if( CutOnOneSide==0 && (h->GetBinContent(i+1)>ave+absolutecut || h->GetBinContent(i+1)<ave-absolutecut) ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
		if( CutOnOneSide==1 && h->GetBinContent(i+1)>ave+absolutecut ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
		if( CutOnOneSide==-1 && h->GetBinContent(i+1)<ave-absolutecut ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
	}

	
	cout<<"Bad run for "<<h->GetTitle()<<endl<<"badrunlist["<<badlist.size()<<"] = {";
	for(list<int>::iterator it = badlist.begin(); it!=badlist.end(); ++it) {
		cout<<*it<<", ";
	}	
	cout<<'\b'; // cursor moves 1 position backwards
	cout<<'\b'; // cursor moves 1 position backwards
	cout<<"};"<<endl;

	TCanvas *c = new TCanvas(Form("c%s",h->GetName()),h->GetTitle());
	gStyle->SetOptStat(0);
	//h->SetMaximum(ave+sig*10);
	//h->SetMinimum(ave-sig*10);
	//h->Draw("pe");
	TH1D *htemp = (TH1D*)hrunid(h,hmap);
	htemp->SetMarkerStyle(h->GetMarkerStyle());
	htemp->SetMarkerColor(h->GetMarkerColor());
	htemp->SetLineStyle(h->GetLineStyle());
	htemp->SetLineColor(h->GetLineColor());
	htemp->Draw("p");
	h->Draw("pesame");
	double xmin = h->GetBinCenter(1);
	double xmax = h->GetBinCenter(h->GetNbinsX()+1);
	TLine *l = new TLine();
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->DrawLine(xmin,ave,xmax,ave);
	l->SetLineStyle(1);
	if( CutOnOneSide>=0 ) {
		l->DrawLine(xmin,ave+absolutecut,xmax,ave+absolutecut);
	}
	if( CutOnOneSide<=0 ) {
		l->DrawLine(xmin,ave-absolutecut,xmax,ave-absolutecut);
	}

	c->SaveAs(Form("%sQA_%s_%s%s.png",dir,h->GetName(),tag_trig,outfiletag.Data()));
}


void BadRun(TH1D *h, TH1D *hmap, list<int> & badlist, const char *tag_trig = "NPE25", TString outfiletag="", const char *dir = "$HOME/scratch/mapBEMCauau11Pico/", double sigmacut = 3, int CutOnOneSide = 0, double absolutecut = 0, double inputave = 0) {			
 // abnormal (>sigmacut-sigma) runid for variable in h (vs run index) 
 // If CutOnOneSide == 0:  symmetric cut on each side 
 // If CutOnOneSide == 1:  cut off the one exceed upper limit  
 // If CutOnOneSide == -1: cut off the one lower than low limit 
 // If absolutecut != 0 : use absolutecut as mean+/-absolutecut

	double ave = 0, sig = 0;
	AveSig(h,ave,sig);
	badlist.clear();

	if( fabs(inputave) > 1e-6) {
		ave = inputave;
	}
	absolutecut = fabs(absolutecut);
	if( fabs(absolutecut) < 1e-6) {		// if absolutecut==0, use relative cut
		absolutecut = sig*sigmacut;
	}

	for(int i = 0; i<h->GetNbinsX(); i++) {
		if( CutOnOneSide==0 && (h->GetBinContent(i+1)>ave+absolutecut || h->GetBinContent(i+1)<ave-absolutecut) ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
		if( CutOnOneSide==1 && h->GetBinContent(i+1)>ave+absolutecut ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
		if( CutOnOneSide==-1 && h->GetBinContent(i+1)<ave-absolutecut ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
	}

	
	cout<<"Bad run for "<<h->GetTitle()<<endl<<"badrunlist["<<badlist.size()<<"] = {";
	for(list<int>::iterator it = badlist.begin(); it!=badlist.end(); ++it) {
		cout<<*it<<", ";
	}	
	cout<<'\b'; // cursor moves 1 position backwards
	cout<<'\b'; // cursor moves 1 position backwards
	cout<<"};"<<endl;

	TCanvas *c = new TCanvas(Form("c%s",h->GetName()),h->GetTitle());
	gStyle->SetOptStat(0);
	//h->SetMaximum(ave+sig*10);
	//h->SetMinimum(ave-sig*10);
	//h->Draw("pe");
	TH1D *htemp = (TH1D*)hrunid(h,hmap);
	htemp->SetMarkerStyle(h->GetMarkerStyle());
	htemp->SetMarkerColor(h->GetMarkerColor());
	htemp->SetLineStyle(h->GetLineStyle());
	htemp->SetLineColor(h->GetLineColor());
	htemp->Draw("p");
	h->Draw("pesame");
	double xmin = h->GetBinCenter(1);
	double xmax = h->GetBinCenter(h->GetNbinsX()+1);
	TLine *l = new TLine();
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->DrawLine(xmin,ave,xmax,ave);
	l->SetLineStyle(1);
	if( CutOnOneSide>=0 ) {
		l->DrawLine(xmin,ave+absolutecut,xmax,ave+absolutecut);
	}
	if( CutOnOneSide<=0 ) {
		l->DrawLine(xmin,ave-absolutecut,xmax,ave-absolutecut);
	}

	c->SaveAs(Form("%sQA_%s_%s%s.png",dir,h->GetName(),tag_trig,outfiletag.Data()));
}


void BadRun(TProfile *h, TH1D *hmap, list<int> & badlist, const char *tag_trig = "NPE25", TString outfiletag="", const char *dir = "$HOME/scratch/mapBEMCauau11Pico/", double sigmacut = 3, int CutOnOneSide = 0, double absolutecut = 0, double inputave = 0) {			
 // abnormal (>sigmacut-sigma) runid for variable in h (vs run index) 
 // If CutOnOneSide == 0:  symmetric cut on each side 
 // If CutOnOneSide == 1:  cut off the one exceed upper limit  
 // If CutOnOneSide == -1: cut off the one lower than low limit 
 // If absolutecut != 0 : use absolutecut as mean+/-absolutecut

	double ave = 0, sig = 0;
	AveSig(h,ave,sig);
	badlist.clear();

	if( fabs(inputave) > 1e-6) {
		ave = inputave;
	}
	absolutecut = fabs(absolutecut);
	if( fabs(absolutecut) < 1e-6) {		// if absolutecut==0, use relative cut
		absolutecut = sig*sigmacut;
	}

	for(int i = 0; i<h->GetNbinsX(); i++) {
		if( CutOnOneSide==0 && (h->GetBinContent(i+1)>ave+absolutecut || h->GetBinContent(i+1)<ave-absolutecut) ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
		if( CutOnOneSide==1 && h->GetBinContent(i+1)>ave+absolutecut ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
		if( CutOnOneSide==-1 && h->GetBinContent(i+1)<ave-absolutecut ) {
			badlist.push_back(hmap->GetBinContent(i+1));
		}
	}
	
	cout<<"Bad run for "<<h->GetTitle()<<endl<<"badrunlist["<<badlist.size()<<"] = {";
	for(list<int>::iterator it = badlist.begin(); it!=badlist.end(); ++it) {
		cout<<*it<<", ";
	}	
	cout<<'\b'; // cursor moves 1 position backwards
	cout<<'\b'; // cursor moves 1 position backwards
	cout<<"};"<<endl;

	TCanvas *c = new TCanvas(Form("c%s",h->GetName()),h->GetTitle());
	gStyle->SetOptStat(0);
	//h->SetMaximum(ave+sig*10);
	//h->SetMinimum(ave-sig*10);
	//h->Draw("pe");
	h->SetMarkerStyle(4);
	TH1D *htemp = (TH1D*)hrunid(h,hmap);
	htemp->SetMarkerStyle(h->GetMarkerStyle());
	htemp->SetMarkerColor(h->GetMarkerColor());
	htemp->SetLineStyle(h->GetLineStyle());
	htemp->SetLineColor(h->GetLineColor());
	htemp->Draw("p");
	h->Draw("pesame");
	double xmin = h->GetBinCenter(1);
	double xmax = h->GetBinCenter(h->GetNbinsX()+1);
	TLine *l = new TLine();
	l->SetLineColor(2);
	l->SetLineStyle(2);
	l->DrawLine(xmin,ave,xmax,ave);
	l->SetLineStyle(1);
	if( CutOnOneSide>=0 ) {
		l->DrawLine(xmin,ave+absolutecut,xmax,ave+absolutecut);
	}
	if( CutOnOneSide<=0 ) {
		l->DrawLine(xmin,ave-absolutecut,xmax,ave-absolutecut);
	}

	c->SaveAs(Form("%sQA_%s_%s%s.png",dir,h->GetName(),tag_trig,outfiletag.Data()));
}

void MergeBadRunList(vector<int>& a, const vector<int>& b) {		// merge b into a in order, assume both a and b are in ascending order 
	std::vector<int>::iterator ia;
	std::vector<int>::const_iterator ib;
	for(ib = b.begin(); ib!=b.end(); ib++) {
		for(ia = a.begin(); ia!=a.end(); ia++) {
			if( *ib>*ia ) continue;
			if( *ib==*ia ) break;
			if( *ib<*ia ) {
				ia=a.insert(ia,*ib);
			}
		}
	}
	
}


//void MergeBadRunList(list<int>& a, const list<int>& b) {		// merge b into a in order, assume both a and b is in ascending order , b will keep the same
//	std::list<int>::iterator ia;
//	std::list<int>::const_iterator ib;
//	for(ib = b.begin(); ib!=b.end(); ib++) {
//		for(ia = a.begin(); ia!=a.end(); ia++) {
//			if( *ib>*ia ) continue;
//			if( *ib==*ia ) break;
//			if( *ib<*ia ) {
//				ia=a.insert(ia,*ib);
//			}
//		}
//	}
//	
//}


void MergeBadRunList(list<int>& a, list<int>& b) {		// merge b into a in order, assume both a and b are in ascending order , b will be empty
	a.merge(b);	
	a.unique();
}

TH1D *hrunid(TProfile *h, TH1D *hmap) {	// h is the one vs runindex, not the real run number
        string name = h->GetName();
        TH1D *hout = new TH1D(Form("%s_perrun",name.data()),Form("%s",h->GetTitle()),h->GetNbinsX(),0,h->GetNbinsX());
        hout->GetYaxis()->SetTitle(Form("%s",name.data()));
        hout->GetXaxis()->SetTitle(Form("day"));
        hout->SetMarkerStyle(4);
	int currentday = 0;
	int year = floor(hmap->GetBinContent(1)/1000000);
	//cout<<"year 20"<<year-1<<endl;
        for(int i = 0; i<h->GetNbinsX(); i++) {
                hout->SetBinContent(i+1,h->GetBinContent(i+1));
                double runnumber = hmap->GetBinContent(i+1);
                if(fabs(currentday-(floor(runnumber/1000)-year*1000))>1e-6) {
                        char day[3];
                        sprintf(day,"%d",floor(runnumber/1000-year*1000));
			//cout<<day<<endl;
                        hout->GetXaxis()->SetBinLabel(i+1,day);
                        hout->LabelsOption("hd");
			currentday= floor(runnumber/1000)-year*1000;
                }
        }
        double ave = 0, sig = 0;
        AveSig(hout,ave,sig);
        hout->SetMaximum(ave+10*sig);
        hout->SetMinimum(ave-10*sig);

        return hout;
}



TH1D *hrunid(TProfile *h, TH1D *hmap, TCanvas *c, const char *tag_trig = "NPE25", TString outfiletag="", const char *dir = "$HOME/Scratch/mapBEMCauau11Pico/") {	// h is the one vs runindex, not the real run number
        c->SetLogy(0);
        gStyle->SetOptStat(0);
        string name = h->GetName();
        TH1D *hout = new TH1D(Form("%s_perrun",name.data()),Form("%s",h->GetTitle()),h->GetNbinsX(),0,h->GetNbinsX());
        hout->GetYaxis()->SetTitle(Form("%s",name.data()));
        hout->GetXaxis()->SetTitle(Form("runid"));
        hout->SetMarkerStyle(4);
	int currentday = 0;
	int year = floor(hmap->GetBinContent(1)/1000000);
	//cout<<"year 20"<<year-1<<endl;
        for(int i = 0; i<h->GetNbinsX(); i++) {
                hout->SetBinContent(i+1,h->GetBinContent(i+1));
                double runnumber = hmap->GetBinContent(i+1);
                if(fabs(currentday-(floor(runnumber/1000)-year*1000))>1e-6) {
                        char day[3];
                        sprintf(day,"%d",floor(runnumber/1000-year*1000));
			//cout<<day<<endl;
                        hout->GetXaxis()->SetBinLabel(i+1,day);
                        hout->LabelsOption("hd");
			currentday= floor(runnumber/1000)-year*1000;
                }
        }
        double ave = 0, sig = 0;
        AveSig(hout,ave,sig);
        hout->SetMaximum(ave+5*sig);
        hout->SetMinimum(ave-5*sig);

        hout->Draw("p");


        double ymax = hout->GetMaximum();
        double ymin = hout->GetMinimum();
        double ylat = (ymax-ymin)*0.8+ymin;
        double xline = 0;
        TLine *l = new TLine();
        l->SetLineStyle(2);
        l->SetLineColor(kGreen);
        TLatex *lat = new TLatex();
        lat->SetTextAngle(45);
        lat->SetTextColor(kGreen);
        lat->SetTextFont(62);
        lat->SetTextSize(0.03);
        //xline = hout->GetBinCenter(h->FindBin(12138080.5));
        //l->DrawLine(xline,ymin,xline,ymax);     //add FTPC trigger
        //lat->DrawLatex(xline,ylat,"add FTPC trigger");

        //xline = hout->GetBinCenter(h->FindBin(12140029.5));
        //l->DrawLine(xline,ymin,xline,ymax);     //add future_guadian
        //lat->DrawLatex(xline,ylat,"add future_guadian");

        //xline = hout->GetBinCenter(h->FindBin(12144030.5));
        //l->DrawLine(xline,ymin,xline,ymax);     //add future_guadian again
        //lat->DrawLatex(xline,ylat,"future_guadian again");

        //xline = hout->GetBinCenter(h->FindBin(12145020.5));
        //l->DrawLine(xline,ymin,xline,ymax);     //add new trigger file
        //lat->DrawLatex(xline,ylat,"new trigger file");

	c->SaveAs(Form("%sQA_%s_%s%s.png",dir,h->GetName(),tag_trig,outfiletag.Data()));
        return hout;

}






//==============================================================
//================== Main Function =============================
//==============================================================
//void TimeDep(TString fin="/home/ly247/code/BEMCHTFinder/AuAu200Run11NPE18Central.list") {
void TimeDep(TString fin="/home/fas/caines/ly247/scratch/run12ppQA/sum*.root",TString outfiletag="",const char *tag_trig="MB") {
//void TimeDep(int fileid = 0) {
//	gSystem->Load("libPhysics");		// needed to TLorentVector
//	gSystem->Load("libHist");		// needed to TLorentVector
//	gSystem->Load("~/Software/PicoCode/eventStructuredAu/libTStarJetPico.so");

//	TString fin="/home/hep/caines/ly247/BEMCHTFinder/AuAu200Run11NPE15.list";  // This is a text file with a list of PicoDst Files, one file per line
//	TString fin="/home/fas/caines/ly247/scratch/run12ppQA/sum";
//	fin+=fileid;
//	fin+=".root";
	TChain *J=new TChain("JetTree");
	TString cName="";

	// Following load root files
	cout<<fin<<endl;
	int fcount = 0;
	if (fin.Contains("root")){
		cName=fin;
		cout<<"Add "<<cName<<endl;	
		J->Add(cName);
		fcount++;
	}
	else{
		if(fin.Contains(".txt")||fin.Contains(".list")){
			std::ifstream in(fin);
			std::string str;
			while (std::getline(in,str))
			{
				cout<<str<<endl;
				TString* strin=new TString(str);
				const char *name;
				name=strin->Data();
				cName=name;
				J->Add(strin->Data());
				fcount++;
				//if(fcount>100) break;	// test
			}
		}
		else{
			cName=fin;
			cName +="pico*.root";
			J->Add(cName);
			fcount++;
		}
		cout<<"read in "<<fcount<<" files"<<endl;
	}
	cout<<"reading finished"<<endl;

	// Creat output root file and histogram
	// NPE-15
	//TFile *fout = new TFile(Form("$HOME/scratch/mapBEMCauau11Pico/BEMCTowerHits_NPE15.root"),"RECREATE");
	//char datadescription[100] = "NPE15 Run11 AuAu 200 GeV";
	//int startrun = 12126079;
	//int endrun = 12171016;
	
	// NPE-25, NPE-18
	//const char *tag_trig="MB";   // "NPE25";
	const char *outdir="/home/fas/caines/ly247/scratch/run12ppQA/out/";
	char datadescription[100];sprintf(datadescription, "%s Run12 pp 200 GeV",tag_trig);
	int startrun = 13039166;
	int endrun = 13104049+1;		// make sure all possible run included
	int runrange = endrun-startrun;		// how many histogram bins

	// per run
	TH1D *hrunid4Run = new TH1D("hrunid4Run",Form("runid vs runid (to filter out void runid) %s",datadescription),runrange,startrun,endrun);
	TProfile *hrefmult4Run = new TProfile("hrefmult4Run",Form("refmult vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hvx4Run = new TProfile("hvx4Run",Form("vx vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hvy4Run = new TProfile("hvy4Run",Form("vy vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hvz4Run = new TProfile("hvz4Run",Form("vz vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hNPVertex4Run = new TProfile("hNPVertex4Run",Form("NumberOfPrimaryVertices vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hratioP2G4Run = new TProfile("hratioP2G4Run",Form("NumberOfPrimaryTrack/GlobalTrack vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hratioTrackMatch4Run = new TProfile("hratioTrackMatch4Run",Form("NOfMatchedTracks/NOfPTracks vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hzdccoinrate4Run = new TProfile("hzdccoinrate4Run",Form("ZDC coincidence Rate vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hbbccoinrate4Run = new TProfile("hbbccoinrate4Run",Form("BBC coincidence Rate vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hbemcNOfTower4Run = new TProfile("hbemcNOfTower4Run",Form("BEMC Number Of Towers vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hbemcratioNOfTowerMatch4Run = new TProfile("hbemcratioNOfTowerMatch4Run",Form("BEMC Number Of Matched Towers/Number Of Towers vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hbemcadcsum4Run = new TProfile("hbemcadcsum4Run",Form("BBEMC ADC sum vs runid %s",datadescription),runrange,startrun,endrun);

	TProfile *hglobal4Run = new TProfile("hglobal4Run",Form("global vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hgoodtrk4Run = new TProfile("hgoodtrk4Run",Form("goodtrk vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hgoodtrkbemcmatch4Run = new TProfile("hgoodtrkbemcmatch4Run",Form("goodtrkbemcmatch vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hgoodtrktofmatch4Run = new TProfile("hgoodtrktofmatch4Run",Form("goodtrktofmatch vs runid %s",datadescription),runrange,startrun,endrun);

	TProfile *hdca4Run = new TProfile("hdca4Run",Form("dca vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hpt4Run = new TProfile("hpt4Run",Form("pt vs runid %s",datadescription),runrange,startrun,endrun);
	TProfile *hEt4Run = new TProfile("hEt4Run",Form("Et vs runid %s",datadescription),runrange,startrun,endrun);
	
	// all run
	TH1D *hvz = new TH1D("hvz",Form("Vz for all run %s",datadescription),10000,-100,100);
	TH2D *hvxy = new TH2D("hvxy",Form("Vx - Vy for all run %s",datadescription),1000,-4,4,1000,-4,4);
	TH1D *hzdc = new TH1D("hzdc",Form("Zdc Coin. Rate for all run %s",datadescription),1000,0,20000);
	TH1D *hbbc = new TH1D("hbbc",Form("Bbc Coin. Rate for all run %s",datadescription),1000,0,1500000);
	TH2D *hBbcVsZdc = new TH2D("hBbcVsZdc",Form("hBbcVsZdc  for all run %s",datadescription),500,0,20000,500,0,1500000);
	TH2D *hNOfGlobalVsZdc = new TH2D("hNOfGlobalVsZdc",Form("hNOfGlobalVsZdc for all run %s",datadescription),500,0,20000,500,0,1500);
	TH2D *hgoodtrktofmatchVsNOfGlobal = new TH2D("hgoodtrktofmatchVsNOfGlobal",Form("hgoodtrktofmatchVsNOfGlobal for all run %s",datadescription),500,0,1500,10,0,10);
	TH1D *hbemc = new TH1D("hbemc",Form("BEMC Tower ADC for all run %s",datadescription),4800,1,4801);
	TH2D *hbemc2d = new TH2D("hbemc2d",Form("BEMC Tower ADC for all run %s",datadescription),4800,1,4801,1000,0,4096);	// 2^12, ADC 12 bits
	TH2D *hbemc2d_trig = new TH2D("hbemc2d_trig",Form("Triggered BEMC Tower ADC for all run %s",datadescription),4800,1,4801,1000,0,4096);	// hbemc2d for trig id tower only 

	//TH2D *hbemc4Run = new TH2D("hbemc4Run",Form("bemc vs runid %s",datadescription),runrange,startrun,endrun,10000,-100,100);


	//  Read in Tree
	cout<<"read in tree"<<endl;
	TStarJetPicoEvent *mEv = new TStarJetPicoEvent();
	J->SetBranchAddress("PicoJetTree",&mEv);

	TStarJetPicoEventCuts *EvCut = new TStarJetPicoEventCuts();
	EvCut->SetTriggerSelection(Form("pp%s",tag_trig));

	Int_t ievents = J->GetEntries();

	cout<<"---> Number of events = "<<ievents<<endl;

	// Loop over all events
	for (Int_t i=0; i<ievents; i++){
		//for (Int_t i=0; i<10; i++){
		if( i%10000 == 0){
			cout << "On event " << i << endl;
		}
		//if( i>2 ) break; 	// test
		J->GetEvent(i);

		// Event trigcut
		if(!EvCut->IsTriggerIdOK(mEv)) continue;


		// Event information
		int runid = mEv->GetHeader()->GetRunId();
		hvz->Fill(mEv->GetHeader()->GetPrimaryVertexZ());
		hvxy->Fill(mEv->GetHeader()->GetPrimaryVertexX(),mEv->GetHeader()->GetPrimaryVertexY());
		hzdc->Fill(mEv->GetHeader()->GetZdcCoincidenceRate());
		hBbcVsZdc->Fill(mEv->GetHeader()->GetZdcCoincidenceRate(),mEv->GetHeader()->GetBbcCoincidenceRate());
		hNOfGlobalVsZdc->Fill(mEv->GetHeader()->GetZdcCoincidenceRate(),mEv->GetHeader()->GetNGlobalTracks());
		hbbc->Fill(mEv->GetHeader()->GetBbcCoincidenceRate());

		hrunid4Run->Fill(runid);
		hrefmult4Run->Fill(runid,mEv->GetHeader()->GetReferenceMultiplicity());
		hvx4Run->Fill(runid,mEv->GetHeader()->GetPrimaryVertexX());
		hvy4Run->Fill(runid,mEv->GetHeader()->GetPrimaryVertexY());
		hvz4Run->Fill(runid,mEv->GetHeader()->GetPrimaryVertexZ());
		hNPVertex4Run->Fill(runid,mEv->GetHeader()->GetNumberOfVertices());
		hratioP2G4Run->Fill(runid,(mEv->GetHeader()->GetNGlobalTracks()==0)?0:(1.*mEv->GetHeader()->GetNOfPrimaryTracks()/mEv->GetHeader()->GetNGlobalTracks()));
		hratioTrackMatch4Run->Fill(runid,(mEv->GetHeader()->GetNOfPrimaryTracks()==0)?0:(1.*mEv->GetHeader()->GetNOfMatchedTracks()/mEv->GetHeader()->GetNOfPrimaryTracks()));
		hzdccoinrate4Run->Fill(runid,mEv->GetHeader()->GetZdcCoincidenceRate());
		hbbccoinrate4Run->Fill(runid,mEv->GetHeader()->GetBbcCoincidenceRate());
		hbemcNOfTower4Run->Fill(runid,mEv->GetHeader()->GetNOfTowers());
		hbemcratioNOfTowerMatch4Run->Fill(runid,(mEv->GetHeader()->GetNOfTowers()==0)?0:1.*mEv->GetHeader()->GetNOfMatchedTowers()/mEv->GetHeader()->GetNOfTowers());
		hglobal4Run->Fill(runid,mEv->GetHeader()->GetNGlobalTracks());
		

		// Get the Tower trigger the HT event
		for(int itrg = 0; itrg<mEv->GetHeader()->GetNOfTrigObjs(); itrg ++) {
			hbemc2d_trig->Fill(mEv->GetTrigObj(itrg)->GetADC(),mEv->GetTrigObj(itrg)->GetId());
		}


		// Read BEMC Tower
		for (Int_t itwr=0; itwr<mEv->GetHeader()->GetNOfTowers() ; itwr++) {
			TStarJetPicoTower *mTwr = mEv->GetTower(itwr);
			hbemc->Fill(mTwr->GetId(),mTwr->GetEnergy());
			hbemc2d->Fill(mTwr->GetId(),mTwr->GetEnergy());
			//cout<<"tower "<<itwr<<"\t"<<mTwr->GetId()<<"\t"<<mTwr->GetADC()<<endl; // test
			hbemcadcsum4Run->Fill(runid,mTwr->GetEnergy());
			hEt4Run->Fill(runid,mTwr->GetEnergy());
		}

		// Primary track loop
		int goodtrk = 0;
		int goodtrktofmatch = 0;
		int goodtrkbemcmatch = 0;
		int flagEvt = 0;
		if(fabs(mEv->GetHeader()->GetPrimaryVertexZ())<=5 && fabs(mEv->GetHeader()->GetPrimaryVertexZ()-mEv->GetHeader()->GetVpdVz())<=3 && mEv->GetHeader()->GetPrimaryVertexRanking()>0) flagEvt = 1;		// good

		for (Int_t itrk=0; itrk<mEv->GetHeader()->GetNOfPrimaryTracks(); itrk++) {
			TStarJetPicoPrimaryTrack *mTrk = mEv->GetPrimaryTrack(itrk); 
			double dca = mTrk->GetDCA();
			double px = mTrk->GetPx();
			double py = mTrk->GetPy();
			double pz = mTrk->GetPz();
			double pt = sqrt(px*px+py*py+pz*pz);
			
			hdca4Run->Fill(runid,dca);
			hpt4Run->Fill(runid,pt);

			if(flagEvt>0) {
				if(mTrk->GetFlag()<=0) continue;
				if(fabs(mTrk->GetEta())>0.5) continue;
				if(pt<0.15&&pt>10) continue;
				if(dca>1) continue;
				float nfit = mTrk->GetNOfFittedHits();
				float nposs = mTrk->GetNOfPossHits();
				if(nfit<25) continue;
				if(1.0*nfit/nposs<0.52) continue;
				if(1.0*nfit/nposs>1.02) continue;

				goodtrk++;
				if(mTrk->GetBemcMatchFlag()) goodtrkbemcmatch++;
				if(mTrk->GetTofMatchFlag()) goodtrktofmatch++;

			}
			hgoodtrk4Run->Fill(runid,goodtrk);
			hgoodtrkbemcmatch4Run->Fill(runid,goodtrkbemcmatch);
			hgoodtrktofmatch4Run->Fill(runid,goodtrktofmatch);

			hgoodtrktofmatchVsNOfGlobal->Fill(mEv->GetHeader()->GetNGlobalTracks(),goodtrktofmatch);
		}

	}// End of Evt Loop

// Convert run number into runindex array, so that we don't have gap from not used runnumber in the histogram x
	vector<int> runindex;
	vector<int> runindex4map;
	for(int ir = 0; ir<runrange; ir++) {
		if(hrunid4Run->GetBinContent(ir+1)>0) {		// valid runnumber
			runindex.push_back(ir+1);
			runindex4map.push_back(hrunid4Run->GetBinLowEdge(ir+1));
		}
	}
	int runindexsize = runindex.size();
	TH1D *hmaprunindex2runid = new TH1D("hmaprunindex2runid","Map runindex to runid",runindexsize,0,runindexsize);
	for(unsigned int idx=0;idx<runindex4map.size(); idx++) {
		hmaprunindex2runid->SetBinContent(idx+1,runindex4map.at(idx));
	}
	// per run rearranged as runindex (0...number of runs) to avoid the gap from the run number not in use
	TProfile *hrefmult4runindex = new TProfile("hrefmult4runindex",Form("refmult vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hvx4runindex = new TProfile("hvx4runindex",Form("vx vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hvy4runindex = new TProfile("hvy4runindex",Form("vy vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hvz4runindex = new TProfile("hvz4runindex",Form("vz vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hNPVertex4runindex = new TProfile("hNPVertex4runindex",Form("NumberOfPrimaryVertices vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hratioP2G4runindex = new TProfile("hratioP2G4runindex",Form("NumberOfPrimaryTrack/GlobalTrack vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hratioTrackMatch4runindex = new TProfile("hratioTrackMatch4runindex",Form("NOfMatchedTracks/NOfPTracks vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hzdccoinrate4runindex = new TProfile("hzdccoinrate4runindex",Form("ZDC coincidence Rate vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hbbccoinrate4runindex = new TProfile("hbbccoinrate4runindex",Form("BBC coincidence Rate vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hbemcNOfTower4runindex = new TProfile("hbemcNOfTower4runindex",Form("BEMC Number Of Towers vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hbemcratioNOfTowerMatch4runindex = new TProfile("hbemcratioNOfTowerMatch4runindex",Form("BEMC Number Of Matched Towers/Number Of Towers vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hbemcadcsum4runindex = new TProfile("hbemcadcsum4runindex",Form("BBEMC ADC sum vs runid %s",datadescription),runindexsize,0,runindexsize);
	
	TProfile *hglobal4runindex = new TProfile("hglobal4runindex",Form("global vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hgoodtrk4runindex = new TProfile("hgoodtrk4runindex",Form("goodtrk vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hgoodtrkbemcmatch4runindex = new TProfile("hgoodtrkbemcmatch4runindex",Form("goodtrkbemcmatch vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hgoodtrktofmatch4runindex = new TProfile("hgoodtrktofmatch4runindex",Form("goodtrktofmatch vs runid %s",datadescription),runindexsize,0,runindexsize);

	TProfile *hdca4runindex = new TProfile("hdca4runindex",Form("dca vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hpt4runindex = new TProfile("hpt4runindex",Form("pt vs runid %s",datadescription),runindexsize,0,runindexsize);
	TProfile *hEt4runindex = new TProfile("hEt4runindex",Form("Et vs runid %s",datadescription),runindexsize,0,runindexsize);

	convertrunindex(runindex,hrefmult4Run,hrefmult4runindex);
	convertrunindex(runindex,hvx4Run,hvx4runindex);
	convertrunindex(runindex,hvy4Run,hvy4runindex);
	convertrunindex(runindex,hvz4Run,hvz4runindex);
	convertrunindex(runindex,hNPVertex4Run,hNPVertex4runindex);
	convertrunindex(runindex,hratioP2G4Run,hratioP2G4runindex);
	convertrunindex(runindex,hratioTrackMatch4Run,hratioTrackMatch4runindex);
	convertrunindex(runindex,hzdccoinrate4Run,hzdccoinrate4runindex);
	convertrunindex(runindex,hbbccoinrate4Run,hbbccoinrate4runindex);
	convertrunindex(runindex,hbemcNOfTower4Run,hbemcNOfTower4runindex);
	convertrunindex(runindex,hbemcratioNOfTowerMatch4Run,hbemcratioNOfTowerMatch4runindex);
	convertrunindex(runindex,hbemcadcsum4Run,hbemcadcsum4runindex);
	convertrunindex(runindex,hglobal4Run,hglobal4runindex);

	convertrunindex(runindex,hgoodtrk4Run,hgoodtrk4runindex);
	convertrunindex(runindex,hgoodtrkbemcmatch4Run,hgoodtrkbemcmatch4runindex);
	convertrunindex(runindex,hgoodtrktofmatch4Run,hgoodtrktofmatch4runindex);

	convertrunindex(runindex,hdca4Run,hdca4runindex);
	convertrunindex(runindex,hpt4Run,hpt4runindex);
	convertrunindex(runindex,hEt4Run,hEt4runindex);

// Write out hot tower id into txt file
	ofstream ftxt;
	ftxt.open(Form("list4BEMCTowerHits_%s.txt",tag_trig));

	double bemcave=0;
	for(int i = 1; i<=4800; i++) {
		bemcave+=hbemc->GetBinContent(i);	
	}
	bemcave/=4800;
	double bemcsigma=0;
	for(int i = 1; i<=4800; i++) {
		bemcsigma+=pow((hbemc->GetBinContent(i)-bemcave),2);
	}
	bemcsigma/=4800-1;
	bemcsigma=sqrt(bemcsigma);
	double times=3;			// times-bemcsigma deviation
	//double ratio = 0.2;		// +- ratio in the good tower range
	cout<<"Hot Tower Id: "<<endl;
	for(int i = 1; i<=4800; i++) {
		//if(hbemc->GetBinContent(i)>bemcave*(1+ratio)) {
		if(hbemc->GetBinContent(i)>bemcave+times*bemcsigma) {
			cout<<i<<endl;
			ftxt<<i<<endl;	
		}
	}
	cout<<endl;
	ftxt.close();

// Draw bemc
	TCanvas *c = new TCanvas();
	hbemc2d->Draw("col");
	c->SetLogz();
	TLine *l = new TLine();
	l->SetLineColor(2);
	l->SetLineStyle(2);
	//l->DrawLine(1,ave*(1+ratio),4801,ave*(1+ratio));
	//l->DrawLine(1,bemcave*(1-ratio),4801,bemcave*(1-ratio));
	l->DrawLine(1,bemcave+times*bemcsigma,4801,bemcave+times*bemcsigma);
	l->DrawLine(1,bemcave-times*bemcsigma,4801,bemcave-times*bemcsigma);
	c->SaveAs(Form("%sQA_%s_%s.png",outdir,"bemc2d",tag_trig));

	c->Clear();
	hbemc2d_trig->Draw("col");
	c->SaveAs(Form("%sQA_%s_%s.png",outdir,"bemc2d_trg",tag_trig));

	c->Clear();
	TH2D *hbemc2d_notrig = (TH2D*)hbemc2d->Clone("hbemc2d_notrig");
	hbemc2d_notrig->Add(hbemc2d_trig,-1);
	hbemc2d_notrig->SetTitle(Form("Non-triggered BEMC Tower ADC for all run %s",datadescription));
	hbemc2d_notrig->Draw("col");
	l->DrawLine(1,bemcave+times*bemcsigma,4801,bemcave+times*bemcsigma);
	l->DrawLine(1,bemcave-times*bemcsigma,4801,bemcave-times*bemcsigma);
	c->SaveAs(Form("%sQA_%s_%s.png",outdir,"bemc2d_nontrg",tag_trig));
	
	
// Search for abnormal run
	//vector<int> refmultlist, vxlist, vylist, vzlist, ratioP2Glist, ratioTrackMatchlist, bemcratioNOfTowerMatchlist; 
	list<int> refmultlist, vxlist, vylist, vzlist, ratioP2Glist, ratioTrackMatchlist, bemcratioNOfTowerMatchlist; 
	BadRun(hrefmult4runindex,hmaprunindex2runid,refmultlist,tag_trig,outfiletag,outdir);
	BadRun(hvx4runindex,hmaprunindex2runid,vxlist,tag_trig,outfiletag,outdir);
	BadRun(hvy4runindex,hmaprunindex2runid,vylist,tag_trig,outfiletag,outdir);
	BadRun(hvz4runindex,hmaprunindex2runid,vzlist,tag_trig,outfiletag,outdir);
	BadRun(hratioP2G4runindex,hmaprunindex2runid,ratioP2Glist,tag_trig,outfiletag,outdir);
	BadRun(hratioTrackMatch4runindex,hmaprunindex2runid,ratioTrackMatchlist,tag_trig,outfiletag,outdir);
	BadRun(hbemcratioNOfTowerMatch4runindex,hmaprunindex2runid,bemcratioNOfTowerMatchlist,tag_trig,outfiletag,outdir);
// Merge bad run list from all variables
	//vector<int> badrunlist;
	list<int> badrunlist;
	MergeBadRunList(badrunlist,refmultlist);
	MergeBadRunList(badrunlist,vxlist);
	MergeBadRunList(badrunlist,vylist);
	MergeBadRunList(badrunlist,vzlist);
	MergeBadRunList(badrunlist,ratioP2Glist);
	MergeBadRunList(badrunlist,ratioTrackMatchlist);
	MergeBadRunList(badrunlist,bemcratioNOfTowerMatchlist);
	
	cout<<endl<<"Potential bad run list: "<<endl;
	list<int>::iterator ilist;
	for(ilist = badrunlist.begin(); ilist!=badrunlist.end(); ilist++) {
		cout<<*ilist<<endl;
	}
	cout<<endl;


// Write histograms into root file
	//TFile *fout = new TFile(Form("%s%s_%d.root",outdir,tag_trig,fileid),"RECREATE");
	TFile *fout = new TFile(Form("%s%s.root",outdir,tag_trig),"RECREATE");
	fout->cd();
	hrunid4Run->Write();
	hmaprunindex2runid->Write();
	hvz->Write();
	hvxy->Write();
	hzdc->Write();
	hbbc->Write();
	hBbcVsZdc->Write();
	hNOfGlobalVsZdc->Write();
	hgoodtrktofmatchVsNOfGlobal->Write();
	hbemc->Write();
	hbemc2d->Write();
	hbemc2d_trig->Write();
	hrefmult4Run->Write();
	hvx4Run->Write();
	hvy4Run->Write();
	hvz4Run->Write();
	hNPVertex4Run->Write();
	hratioP2G4Run->Write();
	hratioTrackMatch4Run->Write();
	hzdccoinrate4Run->Write();
	hbbccoinrate4Run->Write();
	hbemcNOfTower4Run->Write();
	hbemcratioNOfTowerMatch4Run->Write();
	hbemcadcsum4Run->Write();
	hglobal4Run->Write();

	hgoodtrk4Run->Write();
	hgoodtrkbemcmatch4Run->Write();
	hgoodtrktofmatch4Run->Write();

	hdca4Run->Write();
	hpt4Run->Write();
	hEt4Run->Write();	

	hrefmult4runindex->Write();
	hvx4runindex->Write();
	hvy4runindex->Write();
	hvz4runindex->Write();
	hNPVertex4runindex->Write();
	hratioP2G4runindex->Write();
	hratioTrackMatch4runindex->Write();
	hzdccoinrate4runindex->Write();
	hbbccoinrate4runindex->Write();
	hbemcNOfTower4runindex->Write();
	hbemcratioNOfTowerMatch4runindex->Write();
	hbemcadcsum4runindex->Write();
	hglobal4runindex->Write();

	hgoodtrk4runindex->Write();
	hgoodtrkbemcmatch4runindex->Write();
	hgoodtrktofmatch4runindex->Write();

	hdca4runindex->Write();
	hpt4runindex->Write();
	hEt4runindex->Write();	


	fout->Close();

}




