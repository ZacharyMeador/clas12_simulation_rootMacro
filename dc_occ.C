#include <cstdlib>
#include <iostream>
#include <cmath>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLine.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include "TFile.h"
#include <TLegend.h>
#include "TApplication.h"
#include <TROOT.h>

using namespace std;
TApplication gui("GUI",0,NULL);

int main() {

    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(10);

    gStyle->SetPadBorderMode(0);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadColor(10);


    gStyle->SetTitleFont(72,"X");
    gStyle->SetTitleFont(72,"Y");
    gStyle->SetTitleOffset(0.9,"X");
    gStyle->SetTitleOffset(1.2,"Y");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleSize(0.05,"Y");
    
    gStyle->SetLabelFont(72,"X");
    gStyle->SetLabelFont(72,"Y");
    gStyle->SetLabelFont(72,"Z");
//    gStyle->SetOptFit(111);
    gStyle->SetOptStat("nemriou");
    gStyle->SetOptStat("");

    gStyle->SetPalette(1);

    // generated
    int   ngen;
    vector<double>  *pid_gen = new vector<double>;
    // dc
    const int NDCHITMAX = 100000;
    int   ndchit;
    vector<double> *sector = new vector<double>;
    vector<double>  *layer = new vector<double>;
    vector<double>   *wire = new vector<double>;
    vector<double>   *hitn = new vector<double>;
    vector<double>   *pid  = new vector<double>;
    vector<double>   *mpid = new vector<double>;
    vector<double>   *Edep = new vector<double>;
    vector<double>      *E = new vector<double>;
    vector<double>     *vx = new vector<double>;
    vector<double>     *vy = new vector<double>;
    vector<double>     *vz = new vector<double>;
    vector<double>    *mvx = new vector<double>;
    vector<double>    *mvy = new vector<double>;
    vector<double>    *mvz = new vector<double>;
	

    
    // GENERATED
    TChain *g= new TChain("generated");
    g->Add("out*.root");
    g->SetBranchAddress("pid",&pid_gen);
    
    
    // DC
    TChain *d= new TChain("dc");
    d->Add("out*.root");
    d->SetBranchAddress("wire",&wire);
    d->SetBranchAddress("layer",&layer);
    d->SetBranchAddress("sector",&sector);
    d->SetBranchAddress("hitn",&hitn);
    d->SetBranchAddress("trackE",&E);
    d->SetBranchAddress("totEdep",&Edep);
    d->SetBranchAddress("pid",&pid);
    d->SetBranchAddress("mpid",&mpid);
    d->SetBranchAddress("mvx",&mvx);
    d->SetBranchAddress("mvy",&mvy);
    d->SetBranchAddress("mvz",&mvz);
    d->SetBranchAddress("vx",&vx);
    d->SetBranchAddress("vy",&vy);
    d->SetBranchAddress("vz",&vz);
   
    
    Long64_t nentries = d->GetEntries();
    cout << nentries << endl;




    // Create histos
    TH2F *hi_dcocc_all   = new TH2F("DC Occ. All", "DC Occ. All",112, 1.,113., 36, 1.,37.);
    hi_dcocc_all->GetXaxis()->SetTitle("wire");
    hi_dcocc_all->GetYaxis()->SetTitle("layer");

    TH2F *hi_dcocc_ecut   = new TH2F("DC Occ. E > 50 eV", "DC Occ. E > 50 eV",112, 1.,113., 36, 1.,37.);
    hi_dcocc_ecut->GetXaxis()->SetTitle("wire");
    hi_dcocc_ecut->GetYaxis()->SetTitle("layer");
    
    char *histname = new char[50];
    TH1F *hi_dcocc_region[3];
    for(int i=0; i<3; i++) {
      //  sprintf(histname,"hi_dcocc_region%i",i);
	hi_dcocc_region[i]= new TH1F(Form("hi_dcocc_region%i",i+1),"",6,0.5,6.5);
	hi_dcocc_region[i]->GetXaxis()->SetTitle("Sector");
	hi_dcocc_region[i]->GetYaxis()->SetTitle("Occupancy (%)");
    }

    TH2F *hi_bg_origin   = new TH2F("Origin of Bg", "",160, -500.,3500., 160, 0.,800.);
    hi_bg_origin->SetTitle("Origin of Bg");
    hi_bg_origin->GetXaxis()->SetTitle("z(mm)");                                        
    hi_bg_origin->GetYaxis()->SetTitle("r(mm)");                                        

    TH2F *hi_bg_r_vs_z_reg[3];
    TH2F *hi_bg_r_vs_z_vs_ene_reg[3];
    TH2F *hi_bg_r_vs_z_vs_ene_reg_temp[3];
    TH2F *hi_bg_y_vs_x_reg[3];
    
    
    for(int i=0; i<3; i++) {
       // sprintf(histname,"hi_bg_r_vs_z_region%i",i);
        hi_bg_r_vs_z_reg[i]= new TH2F(Form("hi_bg_r_vs_z_region%i",i+1),Form("hi_bg_r_vs_z_region%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_reg[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_reg[i]->GetYaxis()->SetTitle("r(mm)");
      //  sprintf(histname,"hi_bg_y_vs_x_region%i",i);
        hi_bg_y_vs_x_reg[i]= new TH2F(Form("hi_bg_y_vs_x_region%i",i+1), Form("hi_bg_y_vs_x_region%i",i+1),100, -1000.,1000., 100, -1000.,1000.);
        hi_bg_y_vs_x_reg[i]->GetXaxis()->SetTitle("x(mm)");
        hi_bg_y_vs_x_reg[i]->GetYaxis()->SetTitle("y(mm)");
        
        hi_bg_r_vs_z_vs_ene_reg[i]= new TH2F(Form("hi_bg_r_vs_z_vs_ene_reg%i",i+1),Form("hi_bg_r_vs_z_vs_ene_reg%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_vs_ene_reg[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_vs_ene_reg[i]->GetYaxis()->SetTitle("r(mm)");
        hi_bg_r_vs_z_vs_ene_reg[i]->GetZaxis()->SetTitle("Energy(MeV)");
        

        hi_bg_r_vs_z_vs_ene_reg_temp[i]= new TH2F(Form("hi_bg_r_vs_z_vs_ene_reg_temp%i",i+1),Form("hi_bg_r_vs_z_vs_ene_reg_temp%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_vs_ene_reg_temp[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_vs_ene_reg_temp[i]->GetYaxis()->SetTitle("r(mm)");
        hi_bg_r_vs_z_vs_ene_reg_temp[i]->GetZaxis()->SetTitle("Energy(MeV)");
        
    }

    TH2F *hi_bg_energy_tmp   = new TH2F("Energy of Bg tmp", "",160, -500.,3500., 160, 0.,800.);
    hi_bg_energy_tmp->GetXaxis()->SetTitle("z(mm)");
    hi_bg_energy_tmp->GetYaxis()->SetTitle("r(mm)");

    TH2F *hi_bg_energy   = new TH2F("Energy of Bg", "Energy of Bg",160, -500.,3500., 160, 0.,800.);
    hi_bg_energy->GetXaxis()->SetTitle("z(mm)");
    hi_bg_energy->GetYaxis()->SetTitle("r(mm)");
    hi_bg_energy->GetZaxis()->SetTitle("Energy(MeV)");
    
    TH1F *hi_bg_z   = new TH1F("Vz of Bg", "Vz of Bg",400, -500.,3500.);
    //hi_bg_z->SetTitle("Vz of Bg");
    hi_bg_z->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z->GetYaxis()->SetTitle("Rate (MHz)");

    TH1F *hi_bg_z_reg[3];
    TH1F *hi_bg_z_e_reg[3];
    TH1F *hi_bg_z_g_reg[3];
    TH1F *hi_bg_z_p_reg[3];
    TH1F *hi_bg_z_pi_reg[3];
    TH1F *hi_bg_z_n_reg[3];
    TH1F *hi_bg_z_o_reg[3];
    
    for(int i=0; i<3; i++) {
    //    sprintf(histname,"hi_bg_z_region%i",i);
    hi_bg_z_reg[i]= new TH1F(Form("hi_bg_z_region%i",i+1), Form("hi_bg_z_region%i",i+1),100, -200.,6500.);
	hi_bg_z_reg[i]->GetXaxis()->SetTitle("z(mm)");
	hi_bg_z_reg[i]->SetTitle("Rate (MHz)");

        
    hi_bg_z_e_reg[i]= new TH1F(Form("hi_bg_z_e_reg%i",i+1), Form("hi_bg_z_e_reg%i",i+1),100, -200.,6500.);
    hi_bg_z_e_reg[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_e_reg[i]->SetTitle("Rate (MHz)");
        
    hi_bg_z_g_reg[i]= new TH1F(Form("hi_bg_z_g_reg%i",i+1), Form("hi_bg_z_g_reg%i",i+1),100, -200.,6500.);
    hi_bg_z_g_reg[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_g_reg[i]->SetTitle("Rate (MHz)");
    
    hi_bg_z_n_reg[i]= new TH1F(Form("hi_bg_z_n_reg%i",i+1), Form("hi_bg_z_n_reg%i",i+1),100, -200.,6500.);
    hi_bg_z_n_reg[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_n_reg[i]->SetTitle("Rate (MHz)");
        
    hi_bg_z_p_reg[i]= new TH1F(Form("hi_bg_z_p_reg%i",i+1), Form("hi_bg_z_p_reg%i",i+1),100, -200.,6500.);
    hi_bg_z_p_reg[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_p_reg[i]->SetTitle("Rate (MHz)");

    hi_bg_z_pi_reg[i]= new TH1F(Form("hi_bg_z_pi_reg%i",i+1), Form("hi_bg_z_pi_reg%i",i+1),100, -200.,6500.);
    hi_bg_z_pi_reg[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_pi_reg[i]->SetTitle("Rate (MHz)");
        
    hi_bg_z_o_reg[i]= new TH1F(Form("hi_bg_z_o_reg%i",i+1), Form("hi_bg_z_o_reg%i",i+1),100, -200.,6500.);
    hi_bg_z_o_reg[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_o_reg[i]->SetTitle("Rate (MHz)");
        
        
    }


    TH1F *hi_bg_z_e   = new TH1F("Vz of Bg electrons", "",400, -500.,3500.);
    hi_bg_z_e->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_e->GetYaxis()->SetTitle("Rate (MHz)");

    TH1F *hi_bg_z_g   = new TH1F("Vz of Bg Photons", "",400, -500.,3500.);
    hi_bg_z_g->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_g->GetYaxis()->SetTitle("Rate (MHz)");
    
    TH1F *hi_bg_z_o   = new TH1F("Vz of Bg Other", "",400, -500.,3500.);
    hi_bg_z_o->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_o->GetYaxis()->SetTitle("Rate (MHz)");

    TH1F *hi_bg_E  = new TH1F("momentum of Bg", "momentum of Bg",200,0.,10.);
   // hi_bg_E->SetTitle("momentum of Bg");
    hi_bg_E->GetXaxis()->SetTitle("p(MeV)");
    hi_bg_E->GetYaxis()->SetTitle("Rate (MHz)");


// in standard run mode (1 electron at a time)
// float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (124000 electrons in 250ns window)
    float levents=124000;
    int   nsum=124000./levents;
    int   nfull=0;
 
    cout << "number of events to integrate = " << nsum  << endl;

    int nint=0;
    int ngoodentries=0;
    double mass=0;
    double dc_weight;

    for(Long64_t jentry=0; jentry < nentries; jentry++) {
//      for(Long64_t jentry=0; jentry < 1000; jentry++) {

        
        // clear vectors
        pid_gen->clear();
        sector->clear();
        layer->clear();
        wire->clear();
        pid->clear();
        mpid->clear();
        E->clear();
        Edep->clear();
        vx->clear();
        vy->clear();
        vz->clear();
        mvx->clear();
        mvy->clear();
        mvz->clear();
        
        d->GetEntry(jentry);
        g->GetEntry(jentry);
        ngen=pid_gen->size();
        ndchit=sector->size();
        if(ngen>0) ngoodentries++;
        if(int(jentry/1000)*1000==jentry) cout << "Analyzed " << jentry << " events of " << nentries << endl;

	for(int i=0; i<ndchit; i++) {
	    int it=(*hitn)[i]-1;
	    if((*pid)[it]==2112 || (*pid)[it]==2212) {
                mass=938;
            }
            else{
                mass=0;
            }
            int dc_reg=int(((*layer)[i]-1)/12)+1;
            if(dc_reg==1) dc_weight=1;
            else          dc_weight=2;
            hi_dcocc_all->Fill((*wire)[i],(*layer)[i],dc_weight);
	    
	    if((*Edep)[it]>0.00005) {
                hi_dcocc_ecut->Fill((*wire)[i],(*layer)[i],dc_weight);
	        hi_dcocc_region[dc_reg-1]->Fill((*sector)[i],dc_weight);
	        if(dc_reg==1) {
                    hi_bg_origin->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]));
                    hi_bg_z->Fill((*vz)[it]);
                    if((*pid)[it]==11) {
                        hi_bg_z_e->Fill((*vz)[it]);
                    }
                    else if((*pid)[it]==22) {
                        hi_bg_z_g->Fill((*vz)[it]);
                    }
                    else {
                        hi_bg_z_o->Fill((*vz)[it]);
                    }
                    if((*E)[it]>mass) hi_bg_energy_tmp->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]),sqrt((*E)[it]*(*E)[it]-mass*mass));
                  
                    hi_bg_E->Fill(sqrt((*E)[it]*(*E)[it]-mass*mass));
                }
                
                hi_bg_r_vs_z_vs_ene_reg_temp[dc_reg-1]->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]),sqrt((*E)[it]*(*E)[it]-mass*mass));
                hi_bg_r_vs_z_reg[dc_reg-1]->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]));
                
                if((*pid)[it]==11) hi_bg_z_e_reg[dc_reg-1]->Fill((*vz)[it]);
                if((*pid)[it]==22) hi_bg_z_g_reg[dc_reg-1]->Fill((*vz)[it]);
                if((*pid)[it]==2112) hi_bg_z_n_reg[dc_reg-1]->Fill((*vz)[it]);
                if((*pid)[it]==2212) hi_bg_z_p_reg[dc_reg-1]->Fill((*vz)[it]);
                if((*pid)[it]==211 || (*pid)[it]==-211 ) hi_bg_z_pi_reg[dc_reg-1]->Fill((*vz)[it]);
                if((*pid)[it]!=211 && (*pid)[it]!=-211 && (*pid)[it]!=11 && (*pid)[it]!=22 && (*pid)[it]!=2212) hi_bg_z_o_reg[dc_reg-1]->Fill((*vz)[it]);
                
                if((*vz)[it]>1000*(dc_reg+1)) hi_bg_y_vs_x_reg[dc_reg-1]->Fill((*vx)[it],(*vy)[it]);
                hi_bg_z_reg[dc_reg-1]->Fill((*vz)[it]);
                
		}
        }
        
    }

    // calculating normalization factors based on number of good events
    float norm=124000/levents/ngoodentries;
    float lumi=ngoodentries*250/10*levents/124000; // (mbarn^-1 or 10^27 cm^-2)
    float time=250/norm;
    cout << norm/6 << " "<< ngoodentries << endl;

    cout << "normalization factor          = " << norm  << endl;
    cout << "Run time                      = " << time  << " ns" << endl;
    cout << nfull << " " << 1/norm << endl;

    // normalizing rate histogram to 10^35 luminosity
    hi_dcocc_all->Scale(norm/6.);
    hi_dcocc_ecut->Scale(norm/6.);

    
    hi_bg_z->Scale(1000./time);
    hi_bg_z_e->Scale(1000./time);
    hi_bg_z_g->Scale(1000./time);
    hi_bg_z_o->Scale(1000./time);
    hi_bg_E->Scale(1000./time);
    for(int i=0; i<3; i++) {
      hi_bg_z_reg[i]->Scale(1000./time);
      hi_bg_z_e_reg[i]->Scale(1000./time);
      hi_bg_z_g_reg[i]->Scale(1000./time);
      hi_bg_z_p_reg[i]->Scale(1000./time);
      hi_bg_z_pi_reg[i]->Scale(1000./time);
      hi_bg_z_n_reg[i]->Scale(1000./time);
      hi_bg_z_o_reg[i]->Scale(1000./time);
      hi_dcocc_region[i]->Sumw2();
      hi_dcocc_region[i]->Scale(100*norm/112/12);
        
        hi_bg_r_vs_z_vs_ene_reg[i]->Divide(hi_bg_r_vs_z_vs_ene_reg_temp[i],hi_bg_r_vs_z_reg[i]);
    }
    hi_bg_energy->Divide(hi_bg_energy_tmp,hi_bg_origin);

    TCanvas *c1=new TCanvas("c1","Occupancy",750,1000);
    c1->Divide(1,2);
    c1->cd(1);
    //   gPad->SetLogz();
    hi_dcocc_all->Draw("COLZ");
    c1->cd(2);
    //   gPad->SetLogz();
    hi_dcocc_ecut->Draw("COLZ");
    c1->Print("dc_occ.pdf(");
    
    
    TCanvas *m1=new TCanvas("m1","Background Origin Region 1",750,1000);
    m1->Divide(1,3);
    m1->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg[0]->Draw("COLZ");
    m1->cd(2);
    gPad->SetLogz();
    hi_bg_r_vs_z_vs_ene_reg[0]->SetMinimum(0);
    hi_bg_r_vs_z_vs_ene_reg[0]->SetMaximum(100);
    hi_bg_r_vs_z_vs_ene_reg[0]->Draw("COLZ");
    m1->cd(3);
     gPad->SetLogy();
    hi_bg_z_reg[0]->Draw("H");
    
    hi_bg_z_e_reg[0]->SetLineColor(2);
    hi_bg_z_e_reg[0]->Draw("SAMEH");
    
    hi_bg_z_g_reg[0]->SetLineColor(4);
    hi_bg_z_g_reg[0]->Draw("SAMEH");
    
    hi_bg_z_p_reg[0]->SetLineColor(91);
    hi_bg_z_p_reg[0]->Draw("SAMEH");

    hi_bg_z_n_reg[0]->SetLineColor(8);
    hi_bg_z_n_reg[0]->Draw("SAMEH");

    hi_bg_z_pi_reg[0]->SetLineColor(3);
    hi_bg_z_pi_reg[0]->Draw("SAMEH");
    
    hi_bg_z_o_reg[0]->SetLineColor(6);
    hi_bg_z_o_reg[0]->Draw("SAMEH");
    
     TLegend *leg1 = new TLegend(0.7,0.75,0.96,0.96);
      leg1->SetTextSize(.04);
      leg1->AddEntry(hi_bg_z_reg[0],"All","l");
      leg1->AddEntry(hi_bg_z_e_reg[0],"electrons","l");
      leg1->AddEntry(hi_bg_z_g_reg[0],"photons","l");
      leg1->AddEntry(hi_bg_z_p_reg[0],"proton","l");
      leg1->AddEntry(hi_bg_z_pi_reg[0],"pion","l");
      leg1->AddEntry(hi_bg_z_n_reg[0],"neutron","l");
      leg1->AddEntry(hi_bg_z_o_reg[0],"other","l");
      leg1->Draw();
    
    m1->Print("dc_occ.pdf");
    
    TCanvas *m2=new TCanvas("m2","Background Origin Region 2",750,1000);
    m2->Divide(1,3);
    m2->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg[1]->Draw("COLZ");
    m2->cd(2);
    gPad->SetLogz();
    hi_bg_r_vs_z_vs_ene_reg[1]->SetMinimum(0);
    hi_bg_r_vs_z_vs_ene_reg[1]->SetMaximum(100);
    hi_bg_r_vs_z_vs_ene_reg[1]->Draw("COLZ");
    m2->cd(3);
       gPad->SetLogy();
    hi_bg_z_reg[1]->Draw("H");
    
    hi_bg_z_e_reg[1]->SetLineColor(2);
    hi_bg_z_e_reg[1]->Draw("SAMEH");
    
    hi_bg_z_g_reg[1]->SetLineColor(4);
    hi_bg_z_g_reg[1]->Draw("SAMEH");
    
    hi_bg_z_p_reg[1]->SetLineColor(91);
    hi_bg_z_p_reg[1]->Draw("SAMEH");

    hi_bg_z_n_reg[1]->SetLineColor(8);
    hi_bg_z_n_reg[1]->Draw("SAMEH");

    hi_bg_z_pi_reg[1]->SetLineColor(3);
    hi_bg_z_pi_reg[1]->Draw("SAMEH");
    
    hi_bg_z_o_reg[1]->SetLineColor(6);
    hi_bg_z_o_reg[1]->Draw("SAMEH");
    
     TLegend *leg2 = new TLegend(0.7,0.75,0.96,0.96);
      leg2->SetTextSize(.04);
      leg2->AddEntry(hi_bg_z_reg[1],"All","l");
      leg2->AddEntry(hi_bg_z_e_reg[1],"electrons","l");
      leg2->AddEntry(hi_bg_z_g_reg[1],"photons","l");
      leg2->AddEntry(hi_bg_z_p_reg[1],"proton","l");
      leg2->AddEntry(hi_bg_z_pi_reg[1],"pion","l");
      leg2->AddEntry(hi_bg_z_n_reg[1],"neutron","l");
      leg2->AddEntry(hi_bg_z_o_reg[1],"other","l");
      leg2->Draw();
    
    
    m2->Print("dc_occ.pdf");
    
    TCanvas *m3=new TCanvas("m3","Background Origin Region 3",750,1000);
    m3->Divide(1,3);
    m3->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg[2]->Draw("COLZ");
    m3->cd(2);
    gPad->SetLogz();
    hi_bg_r_vs_z_vs_ene_reg[2]->SetMinimum(0);
    hi_bg_r_vs_z_vs_ene_reg[2]->SetMaximum(100);
    hi_bg_r_vs_z_vs_ene_reg[2]->Draw("COLZ");
    m3->cd(3);
    
       gPad->SetLogy();
    
    hi_bg_z_reg[2]->Draw("H");
    
    hi_bg_z_e_reg[2]->SetLineColor(2);
    hi_bg_z_e_reg[2]->Draw("SAMEH");
    
    hi_bg_z_g_reg[2]->SetLineColor(4);
    hi_bg_z_g_reg[2]->Draw("SAMEH");
    
    hi_bg_z_p_reg[2]->SetLineColor(91);
    hi_bg_z_p_reg[2]->Draw("SAMEH");

    hi_bg_z_n_reg[2]->SetLineColor(8);
    hi_bg_z_n_reg[2]->Draw("SAMEH");

    hi_bg_z_pi_reg[2]->SetLineColor(3);
    hi_bg_z_pi_reg[2]->Draw("SAMEH");
    
    hi_bg_z_o_reg[2]->SetLineColor(6);
    hi_bg_z_o_reg[2]->Draw("SAMEH");
    
     TLegend *leg3 = new TLegend(0.7,0.75,0.96,0.96);
      leg3->SetTextSize(.04);
      leg3->AddEntry(hi_bg_z_reg[2],"All","l");
      leg3->AddEntry(hi_bg_z_e_reg[2],"electrons","l");
      leg3->AddEntry(hi_bg_z_g_reg[2],"photons","l");
      leg3->AddEntry(hi_bg_z_p_reg[2],"proton","l");
      leg3->AddEntry(hi_bg_z_pi_reg[2],"pion","l");
      leg3->AddEntry(hi_bg_z_n_reg[2],"neutron","l");
      leg3->AddEntry(hi_bg_z_o_reg[2],"other","l");
      leg3->Draw();
    m3->Print("dc_occ.pdf");
    
    TCanvas *c3=new TCanvas("c3","Background Origin",750,1000);
    c3->Divide(1,2);
    c3->cd(1);
    gPad->SetLogz();
    hi_bg_origin->Draw("COLZ");
    c3->cd(2);
    gPad->SetLogy();
    hi_bg_z->SetLineColor(1);
    //    hi_bg_z->SetMaximum(15);
    hi_bg_z->Draw("H");
    hi_bg_z_e->SetLineColor(2);
    hi_bg_z_e->Draw("HSAME");
    hi_bg_z_g->SetLineColor(4);
    hi_bg_z_g->Draw("HSAME");
    hi_bg_z_o->SetLineColor(3);
    hi_bg_z_o->Draw("HSAME");
    
    
    TLegend *leg = new TLegend(0.7,0.75,0.96,0.96);

    // leg->SetFillStyle(0);
   //  leg ->SetBorderSize(0);
     leg->SetTextSize(.04);
     leg->AddEntry(hi_bg_z,"Vz of Bg","l");
     leg->AddEntry(hi_bg_z_e,"Vz of Bg electrons","l");
     leg->AddEntry(hi_bg_z_g,"Vz of Bg photons","l");
     leg->AddEntry(hi_bg_z_o,"Vz of Bg other","l");
     leg->Draw();
    c3->Print("dc_occ.pdf");
    
    TCanvas *c31=new TCanvas("c31","Background Origin",750,500);
    gPad->SetLogz();
    hi_bg_origin->Draw("COLZ");
    c3->Print("dc_occ_origin.pdf");

    TCanvas *c4=new TCanvas("c4","Background Energy",750,1000);
    c4->Divide(1,2);
    c4->cd(1);
    gPad->SetLogz();
    hi_bg_energy->Draw("COLZ");
    c4->cd(2);
    hi_bg_E->Draw("");
    c4->Print("dc_occ.pdf");
    
    TCanvas *c5=new TCanvas("c5","Region 1 Background Origin",750,1000);
    c5->Divide(1,3);
    c5->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg[0]->Draw("COLZ");
    c5->cd(2);
    gPad->SetLogz();
    hi_bg_y_vs_x_reg[0]->Draw("COLZ");
    c5->cd(3);
    hi_bg_z_reg[0]->Draw("");
    c5->Print("dc_occ.pdf");
    

    TCanvas *c6=new TCanvas("c6","Region 2 Background Origin",750,1000);
    c6->Divide(1,3);
    c6->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg[1]->Draw("COLZ");
    c6->cd(2);
    gPad->SetLogz();
    hi_bg_y_vs_x_reg[1]->Draw("COLZ");
    c6->cd(3);
    hi_bg_z_reg[1]->Draw("");
    c6->Print("dc_occ.pdf");

    TCanvas *c7=new TCanvas("c7","Region 3 Background Origin",750,1000);
    c7->Divide(1,3);
    c7->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg[2]->Draw("COLZ");
    c7->cd(2);
    gPad->SetLogz();
    hi_bg_y_vs_x_reg[2]->Draw("COLZ");
    c7->cd(3);
    hi_bg_z_reg[2]->Draw("");
    c7->Print("dc_occ.pdf");
 
/*
    TCanvas *c8=new TCanvas("c8","Occupancy",750,500);
    hi_dcocc_ecut->Draw("COLZ");
    c8->Print("dc_occ_map.pdf");
    c8->Print("dc_occ.pdf");
*/
    
    TCanvas *c9=new TCanvas("c9","Sector Occupancy",500,500);
    double region_ave_occ[3]={0};
    double region_max_occ[3]={0};
    for(int iy=0; iy<hi_dcocc_ecut->GetYaxis()->GetNbins(); iy++) {
        int dc_reg=int(iy/12);
        for(int ix=0; ix<hi_dcocc_ecut->GetXaxis()->GetNbins(); ix++) {
            region_ave_occ[dc_reg]+=hi_dcocc_ecut->GetBinContent(ix+1,iy+1);
            if(hi_dcocc_ecut->GetBinContent(ix+1,iy+1)>region_max_occ[dc_reg]) region_max_occ[dc_reg]=hi_dcocc_ecut->GetBinContent(ix+1,iy+1);
        }
    }
    TF1 *mypol0[3];
    for(int i=0; i<3; i++) {
        hi_dcocc_region[i]->SetMarkerStyle(20+i);
        hi_dcocc_region[i]->SetMinimum(0.);
        hi_dcocc_region[i]->SetMaximum(6);
        mypol0[i]= new TF1(Form("mypol%i",i),"pol0",-0.5,6.5);
        mypol0[i]->SetLineWidth(1);
        mypol0[i]->SetLineColor(1);
    }
    hi_dcocc_region[0]->SetMarkerColor(kBlue+2);
    hi_dcocc_region[0]->SetLineColor(kBlue+2);
    hi_dcocc_region[0]->Draw("PE");
    mypol0[0]->SetLineColor(kBlue+2);
    hi_dcocc_region[0]->Fit("mypol0","SAME");
    hi_dcocc_region[1]->SetMarkerColor(kRed+2);
    hi_dcocc_region[1]->SetLineColor(kRed+2);
    hi_dcocc_region[1]->Draw("SAMEPE");
    mypol0[1]->SetLineColor(kRed+2);
    hi_dcocc_region[1]->Fit("mypol1","SAME");
    hi_dcocc_region[2]->SetMarkerColor(kGreen+2);
    hi_dcocc_region[2]->SetLineColor(kGreen+2);
    hi_dcocc_region[2]->Draw("SAMEPE");
    mypol0[2]->SetLineColor(kGreen+2);
    hi_dcocc_region[2]->Fit("mypol2","SAME");
    TLegend *l_occ=new TLegend(0.5,0.75,0.96,0.96);
    for(int i=0; i<3; i++) {
      l_occ->AddEntry(hi_dcocc_region[i],Form("region %i occupancy: %.2f %%",i+1,mypol0[i]->GetParameter(0)),"p");
    }
    l_occ->SetTextFont(72);
    l_occ->Draw();
    c9->Print("dc_region_occ.pdf");
    c9->Print("dc_occ.pdf)");
    FILE *fp = fopen("dc_occ.txt","w");
    for(int i=0; i<3;i++) {
        region_ave_occ[i]=region_ave_occ[i]/12.0/hi_dcocc_ecut->GetXaxis()->GetNbins();
	fprintf(fp,"%d   %5.3f  %5.3f  %5.3f  \n",i+1,100*region_ave_occ[i],100*region_max_occ[i],mypol0[i]->GetParameter(0));
    }
    fclose(fp);    

    // saving histograms to file
	TIter next(gDirectory->GetList());
	TObject* obj;
	TH1 *my_hi[1000];
	int ihi=0;
	while(obj= ( TObject* ) next()){
		if(obj->InheritsFrom(TH1::Class())){
			my_hi[ihi] = (TH1*)obj;
			ihi++;
		}
	}
	TFile myfile("dc.root", "recreate");
	for(int i=0; i<ihi; i++) my_hi[i]->Write();
	myfile.Close();  
	
    gui.Run(1);

}   


