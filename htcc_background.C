#include <iostream>
#include <cmath>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLine.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TLegend.h>
#include "TFile.h"
#include "TApplication.h"
#include <TROOT.h>

using namespace std; 
TApplication gui("GUI",0,NULL);

int main() {

   cout << "Setting Style" << endl;
   
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasColor(10);

   gStyle->SetPadBorderMode(0);
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadRightMargin(0.15);
   gStyle->SetPadTopMargin(0.1);
   gStyle->SetPadBottomMargin(0.13);
   gStyle->SetPadColor(10);


   gStyle->SetTitleFont(72,"X");
   gStyle->SetTitleFont(72,"Y");
   gStyle->SetTitleOffset(0.8,"X");
   gStyle->SetTitleOffset(0.9,"Y");
   gStyle->SetTitleSize(0.08,"X");
   gStyle->SetTitleSize(0.08,"Y");

   gStyle->SetLabelFont(72,"X");
   gStyle->SetLabelFont(72,"Y");
   gStyle->SetLabelFont(72,"Z");
   gStyle->SetLabelSize(0.06,"X");
   gStyle->SetLabelSize(0.06,"Y");
   gStyle->SetLabelSize(0.06,"Z");
   gStyle->SetPalette(1);
   gStyle->SetOptFit(111);
   gStyle->SetOptStat("nemriou");
         gStyle->SetOptStat("");

   
// htcc
   int nhtcchit;
   vector<double>   *htcc_sector=new vector<double>;
   vector<double>   *htcc_ring=new vector<double>;
   vector<double>   *htcc_half=new vector<double>;
   vector<double>   *htcc_nphe=new vector<double>;
   vector<double> *htcc_Edep=new vector<double>;
   vector<double> *htcc_E=new vector<double>;
   vector<double> *htcc_x=new vector<double>;
   vector<double> *htcc_y=new vector<double>;
   vector<double> *htcc_z=new vector<double>;
   vector<double> *htcc_lx=new vector<double>;
   vector<double> *htcc_ly=new vector<double>;
   vector<double> *htcc_lz=new vector<double>;
   vector<double> *htcc_t=new vector<double>;
   vector<double>   *htcc_pid=new vector<double>;
   vector<double>   *htcc_mpid=new vector<double>;
   vector<double> *htcc_vx=new vector<double>;
   vector<double> *htcc_vy=new vector<double>;
   vector<double> *htcc_vz=new vector<double>;
   vector<double> *htcc_mvx=new vector<double>;
   vector<double> *htcc_mvy=new vector<double>;
   vector<double> *htcc_mvz=new vector<double>;
   
   cout << "Creating Tree chains" << endl;
   TChain *g= new TChain("generated");
   g->Add("out*.root");
   TChain *htcc= new TChain("htcc");
   htcc->Add("out*.root");
   
// GENERATED
//   g->SetBranchAddress("evn" ,&evn);
//   g->SetBranchAddress("ngen",&ngen);

   // htcc
   htcc->SetBranchAddress("sector" ,&htcc_sector);
   htcc->SetBranchAddress("ring" ,  &htcc_ring);
   htcc->SetBranchAddress("half" ,  &htcc_half);
   htcc->SetBranchAddress("nphe"   ,&htcc_nphe);
   htcc->SetBranchAddress("trackE"      ,&htcc_E);
   htcc->SetBranchAddress("totEdep"   ,&htcc_Edep);
   htcc->SetBranchAddress("avg_t"      ,&htcc_t);
   htcc->SetBranchAddress("pid"    ,&htcc_pid);
   htcc->SetBranchAddress("avg_x"      ,&htcc_x);
   htcc->SetBranchAddress("avg_y"      ,&htcc_y);
   htcc->SetBranchAddress("avg_z"      ,&htcc_z);
   htcc->SetBranchAddress("avg_lx"     ,&htcc_lx);
   htcc->SetBranchAddress("avg_ly"     ,&htcc_ly);
   htcc->SetBranchAddress("avg_lz"     ,&htcc_lz);
   htcc->SetBranchAddress("vx"     ,&htcc_vx);
   htcc->SetBranchAddress("vy"     ,&htcc_vy);
   htcc->SetBranchAddress("vz"     ,&htcc_vz);
//   htcc->SetBranchAddress("mpid"   ,&htcc_mpid);
//   htcc->SetBranchAddress("mvx"    ,&htcc_mvx);
//   htcc->SetBranchAddress("mvy"    ,&htcc_mvy);
//   htcc->SetBranchAddress("mvz"    ,&htcc_mvz);
   
   Long64_t nentries = htcc->GetEntries();
   cout << "N. entries htcc:" << nentries << endl;


// Create histos
   TH2F *hi_htcc_occ = new TH2F("hi_htcc_occ(%)", "",6, 0.,6.,50, 0.,50);
   hi_htcc_occ->GetXaxis()->SetTitle("Ring");
   hi_htcc_occ->GetYaxis()->SetTitle("n. p.e.");

   TH2F *hi_htcc_occ_rate = new TH2F("hi_htcc_occ_rate(kHZ)", "",6, 0.,6.,50, 0.,50);
   hi_htcc_occ_rate->GetXaxis()->SetTitle("Ring");
   hi_htcc_occ_rate->GetYaxis()->SetTitle("n. p.e.");
    
    TH1F *hi_htcc_rate = new TH1F("hi_htcc_rate(kHZ)", "",6, 0.,6.);
    hi_htcc_rate->GetXaxis()->SetTitle("Ring");
    hi_htcc_rate->GetYaxis()->SetTitle("Rate (kHz)");

    TH1F *hi_htcc_rate_cut = new TH1F("hi_htcc_rate_cut(kHZ)", "",6, 0.,6.);
    hi_htcc_rate_cut->GetXaxis()->SetTitle("Ring");
    hi_htcc_rate_cut->GetYaxis()->SetTitle("Rate (kHz)");

    TH1F *hi_htcc_nphe = new TH1F("hi_htcc_nphe", "",100, 0.,100.);
    hi_htcc_nphe->GetXaxis()->SetTitle("n. p.e.");
    hi_htcc_nphe->GetYaxis()->SetTitle("Rate (kHz)");
    
    TH1F *hi_htcc_nphe_ring4 = new TH1F("hi_htcc_nphe_ring4", "",100, 0.,100.);
    hi_htcc_nphe_ring4->GetXaxis()->SetTitle("n. p.e.");
    hi_htcc_nphe_ring4->GetYaxis()->SetTitle("Rate (kHz)");
    hi_htcc_nphe_ring4->SetLineColor(2);

   TH1F *hi_htcc_vz_all = new TH1F("hi_htcc_vz_all", "",200, -50.,150.);
   hi_htcc_vz_all->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_htcc_vz_all->GetYaxis()->SetTitle("Rate (MHz/cm)");

   TH1F *hi_htcc_vz_e = new TH1F("hi_htcc_vz_e", "",200, -50.,150.);
   hi_htcc_vz_e->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_htcc_vz_e->GetYaxis()->SetTitle("Rate (MHz/cm)");

   TH1F *hi_htcc_vz_g = new TH1F("hi_htcc_vz_g", "",200, -50.,150.);
   hi_htcc_vz_g->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_htcc_vz_g->GetYaxis()->SetTitle("Rate (MHz/cm)");

   TH1F *hi_htcc_vz_h = new TH1F("hi_htcc_vz_h", "",200, -50.,150.);
   hi_htcc_vz_h->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_htcc_vz_h->GetYaxis()->SetTitle("Rate (MHz/cm)");

   TH1F *hi_htcc_vz_n = new TH1F("hi_htcc_vz_n", "",200, -50.,150.);
   hi_htcc_vz_n->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_htcc_vz_n->GetYaxis()->SetTitle("Rate (MHz/cm)");

   TH2F *hi_htcc_origin_all = new TH2F("hi_htcc_origin_all", "",200, 0.,200.,200, 0.,100.);
   hi_htcc_origin_all->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_htcc_origin_all->GetYaxis()->SetTitle("R_{vertex} (cm)");

   TH2F *hi_htcc_azimuthal_all = new TH2F("hi_htcc_azimuthal_all", "",200, -100.,100.,200, -100.,100.);
   hi_htcc_azimuthal_all->GetXaxis()->SetTitle("X_{vertex} (cm)");
   hi_htcc_azimuthal_all->GetYaxis()->SetTitle("Y_{vertex} (cm)");
   
   TH1F *hi_htcc_occ_nphe = new TH1F("hi_htcc_occ_nphe (uA)", "",6, 0.,6.);
   hi_htcc_occ_nphe->GetXaxis()->SetTitle("ring");
   hi_htcc_occ_nphe->GetYaxis()->SetTitle("<n. p.e.>");
   




// in standard run mode (1 electron at a time)
//   float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (50*2300 or 118600 electrons in 250ns window)
   float norm=1.0/(nentries);
   float lumi=nentries*250/10; // (mbarn^-1 or 10^27 cm^-2)
   float time=250/norm/1000;

   cout << "normalization factor          = " << norm  << endl;
   cout << "Run time                      = " << time  << " us" << endl;


   float Nthr=3;            // energy thresholds in nphe
   int ngoodentries=0;
   for (Long64_t jentry=0; jentry < nentries; jentry++) {

       htcc_sector->clear();
       htcc_ring->clear();
       htcc_half->clear();
       htcc_nphe->clear();
       htcc_E->clear();
       htcc_Edep->clear();
       htcc_t->clear();
       htcc_pid->clear();
       htcc_x->clear();
       htcc_y->clear();
       htcc_z->clear();
       htcc_lx->clear();
       htcc_ly->clear();
       htcc_lz->clear(); 
       htcc_vx->clear();
       htcc_vy->clear();
       htcc_vz->clear();
       
       int nb = g->GetEntry(jentry); 
       htcc->GetEntry(jentry);
       ngoodentries++;
       if(int(jentry/1000)*1000==jentry) cout << "Analyzed " << jentry << " events of " << nentries << endl;
	
       nhtcchit=htcc_pid->size();
       for(int i=0; i<nhtcchit; i++) {
	 if((*htcc_sector)[i]>=0 && (*htcc_nphe)[i]>0){
	   hi_htcc_occ->Fill((*htcc_ring)[i]*1.,(*htcc_nphe)[i]*1.,1./12.);
           hi_htcc_occ_rate->Fill((*htcc_ring)[i]*1.,(*htcc_nphe)[i]*1.,1./12.);
           hi_htcc_rate->Fill((*htcc_ring)[i]*1.,1./12.);
           hi_htcc_nphe->Fill((*htcc_nphe)[i],1./12.);
           if((*htcc_ring)[i]==4) hi_htcc_nphe_ring4->Fill((*htcc_nphe)[i],1./12.);
           
           hi_htcc_occ_nphe->Fill((*htcc_ring)[i]*1.,(*htcc_nphe)[i]*1/12.);
           hi_htcc_vz_all->Fill((*htcc_vz)[i]/10.);
           if(abs((*htcc_pid)[i])==11) hi_htcc_vz_e->Fill((*htcc_vz)[i]/10.);
           else if((*htcc_pid)[i]==22) hi_htcc_vz_g->Fill((*htcc_vz)[i]/10.);
           else                        hi_htcc_vz_h->Fill((*htcc_vz)[i]/10.);
           if((*htcc_pid)[i]==2112)    hi_htcc_vz_n->Fill((*htcc_vz)[i]/10.);
           hi_htcc_origin_all->Fill((*htcc_vz)[i]/10.,sqrt((*htcc_vx)[i]*(*htcc_vx)[i]/100.+(*htcc_vy)[i]*(*htcc_vy)[i]/100.));
           hi_htcc_azimuthal_all->Fill((*htcc_vx)[i]/10.,(*htcc_vy)[i]/10.);
           if((*htcc_nphe)[i]>Nthr) {
               hi_htcc_rate_cut->Fill((*htcc_ring)[i]*1.,1./12.);
           }
	 }
       }
   }

// in standard run mode (1 electron at a time)
//   float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (50*2300 or 118600 electrons in 250ns window)
   norm=1.0/ngoodentries;
   lumi=ngoodentries*250/10; // (mbarn^-1 or 10^27 cm^-2)
   time=250/norm/1000;

   cout << "normalization factor          = " << norm  << endl;
   cout << "Run time                      = " << time  << " us" << endl;
   
   cout << ngoodentries << " events analyzed" << endl;

   // normalizing rate histogram to 10^35 luminosity
   hi_htcc_occ->Scale(100./ngoodentries);
   hi_htcc_occ_rate->Scale(1E3/time);
   hi_htcc_rate->Scale(1E3/time);
   hi_htcc_rate_cut->Scale(1E3/time);
   hi_htcc_nphe->Scale(1E3/time);
   hi_htcc_nphe_ring4->Scale(1E3/time);
   hi_htcc_occ_nphe->Scale(1/time);
   hi_htcc_vz_all->Scale(1/time);
   hi_htcc_vz_e->Scale(1/time);
   hi_htcc_vz_g->Scale(1/time);
   hi_htcc_vz_h->Scale(1/time);
   hi_htcc_vz_n->Scale(1/time);
   hi_htcc_origin_all->Scale(1/time);
   hi_htcc_azimuthal_all->Scale(1/time);
   

   TCanvas *c_occ=new TCanvas("c_occ","Occupancy",750,1000);
   c_occ->Divide(1,2);
   c_occ->cd(1);
   hi_htcc_occ->Draw("COLZ");
   c_occ->cd(2);
   hi_htcc_occ_rate->Draw("COLZ");
   c_occ->Print("htcc_occupancy.pdf");

    TCanvas *c_rate=new TCanvas("c_rate","Rate",750,500);
    TLegend *l_occ=new TLegend(0.72,0.75,0.96,0.96);
    gPad->SetLogy();
    hi_htcc_rate->SetMarkerColor(kBlue+2);
    hi_htcc_rate->SetMarkerStyle(20);
    hi_htcc_rate->SetLineColor(kBlue+2);
    hi_htcc_rate->SetMinimum(0.5);
    hi_htcc_rate->Draw("PE");
    l_occ->AddEntry(hi_htcc_rate,"N. PhE>0");
    hi_htcc_rate_cut->SetMarkerColor(kRed+2);
    hi_htcc_rate_cut->SetMarkerStyle(20);
    hi_htcc_rate_cut->SetLineColor(kRed+2);
    hi_htcc_rate_cut->Draw("SAMEPE");
    l_occ->AddEntry(hi_htcc_rate_cut,Form("N. PhE>%.0f",Nthr));
    l_occ->SetTextFont(72);
    l_occ->Draw();
    c_rate->Print("htcc_rate.pdf");
    
    TCanvas *c_nphe=new TCanvas("c_nphe","nphe",750,500);
    TLegend *l_nphe=new TLegend(0.72,0.75,0.96,0.96);
    gPad->SetLogy();
    hi_htcc_nphe->Draw();
    hi_htcc_nphe_ring4->Draw("SAME");
    l_nphe->AddEntry(hi_htcc_nphe,"All PMTs");
    l_nphe->AddEntry(hi_htcc_nphe_ring4,"4th ring");
    l_nphe->SetTextFont(72);
    l_nphe->Draw();
    c_nphe->Print("htcc_nphe.pdf");


   TCanvas *c_origin=new TCanvas("c_origin","Origin",750,1000);
   c_origin->Divide(1,2);
   c_origin->cd(1);
   gPad->SetLogz();
   hi_htcc_origin_all->Draw("COLZ");
   c_origin->cd(2);
   gPad->SetLogz();
   hi_htcc_azimuthal_all->Draw("COLZ");
   c_origin->Print("htcc_origin.pdf");

   TCanvas *c_vertex=new TCanvas("c_vertex","vertex",750,500);
   gPad->SetLogz();   
   hi_htcc_vz_all->SetMinimum(0.001);
   hi_htcc_vz_all->SetLineColor(1);
   hi_htcc_vz_all->Draw();
   hi_htcc_vz_e->SetLineColor(2);
   hi_htcc_vz_e->Draw("SAME");
   hi_htcc_vz_g->SetLineColor(4);
   hi_htcc_vz_g->Draw("SAME");
   hi_htcc_vz_h->SetLineColor(kViolet);
   hi_htcc_vz_h->Draw("SAME");
   hi_htcc_vz_n->SetLineColor(kGreen+2);
   hi_htcc_vz_n->Draw("SAME");
   c_vertex->Print("htcc_vertex.pdf");
   
   TCanvas *c_occ_nphe=new TCanvas("c_occ_nphe","Occ_nphe",750,500);
   hi_htcc_occ_nphe->Draw("");
   c_occ_nphe->Print("htcc_occ_nphe.pdf");
   
   gui.Run(1);

}   


