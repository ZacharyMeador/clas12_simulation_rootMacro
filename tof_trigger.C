#include <iostream>
#include <cmath>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLine.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TChain.h>
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
   gStyle->SetPadLeftMargin(0.1);
   gStyle->SetPadRightMargin(0.12);
   gStyle->SetPadTopMargin(0.1);
   gStyle->SetPadBottomMargin(0.13);
   gStyle->SetPadColor(10);


   gStyle->SetTitleFont(72,"X");
   gStyle->SetTitleFont(72,"Y");
   gStyle->SetTitleOffset(0.7,"X");
   gStyle->SetTitleOffset(0.5,"Y");
   gStyle->SetTitleSize(0.05,"X");
   gStyle->SetTitleSize(0.05,"Y");

   gStyle->SetLabelFont(72,"X");
   gStyle->SetLabelFont(72,"Y");
   gStyle->SetLabelFont(72,"Z");
   //   gStyle->SetLabelSize(0.04,"X");
   //   gStyle->SetLabelSize(0.04,"Y");
   //   gStyle->SetLabelSize(0.04,"Z");
   gStyle->SetPalette(1);
   gStyle->SetOptFit(111);
   gStyle->SetOptStat("nemriou");
   gStyle->SetOptStat("");

   int   evn;
   int   ngen;
   double mass=0;

   const int NTOFMAX = 100000;
// FTOF_Panel_1A
   int ntofhit;
   vector<double> *tof_sector=new vector<double>;
   vector<double> *tof_layer=new vector<double>;
   vector<double> *tof_paddle=new vector<double>;
   vector<double> *tof_ADCL=new vector<double>;
   vector<double> *tof_ADCR=new vector<double>;
   vector<double> *tof_TDCL=new vector<double>;
   vector<double> *tof_TDCR=new vector<double>;
   vector<double> *tof_Edep=new vector<double>;
   vector<double> *tof_E=new vector<double>;
   vector<double> *tof_x=new vector<double>;
   vector<double> *tof_y=new vector<double>;
   vector<double> *tof_z=new vector<double>;
   vector<double> *tof_lx=new vector<double>;
   vector<double> *tof_ly=new vector<double>;
   vector<double> *tof_lz=new vector<double>;
   vector<double> *tof_t=new vector<double>;
   vector<double> *tof_pid=new vector<double>;
   vector<double> *tof_mpid=new vector<double>;
   vector<double> *tof_vx=new vector<double>;
   vector<double> *tof_vy=new vector<double>;
   vector<double> *tof_vz=new vector<double>;
   vector<double> *tof_mvx=new vector<double>;
   vector<double> *tof_mvy=new vector<double>;
   vector<double> *tof_mvz=new vector<double>;



   int   hitS[3][6][60][64];
   float edepS[3][6][60][64];
   float timeS[3][6][60][64];
   int   ntof[64];

   cout << "Creating Tree chains" << endl;
   TChain *g= new TChain("generated");
   g->Add("out*.root");
   TChain *tof= new TChain("ftof");
   tof->Add("out*.root");


// GENERATED
//   g->SetBranchAddress("evn" ,&evn);
//   g->SetBranchAddress("ngen",&ngen);

// FTOF_Panel_1A
   tof->SetBranchAddress("sector" ,&tof_sector);
   tof->SetBranchAddress("layer"  ,&tof_layer);
   tof->SetBranchAddress("paddle" ,&tof_paddle);
   tof->SetBranchAddress("ADCL"   ,&tof_ADCL);
   tof->SetBranchAddress("ADCR"   ,&tof_ADCR);
   tof->SetBranchAddress("TDCL"   ,&tof_TDCL);
   tof->SetBranchAddress("TDCR"   ,&tof_TDCR);
   tof->SetBranchAddress("trackE" ,&tof_E);
   tof->SetBranchAddress("totEdep",&tof_Edep);
   tof->SetBranchAddress("avg_t"  ,&tof_t);
   tof->SetBranchAddress("pid"    ,&tof_pid);
   tof->SetBranchAddress("avg_x"  ,&tof_x);
   tof->SetBranchAddress("avg_y"  ,&tof_y);
   tof->SetBranchAddress("avg_z"  ,&tof_z);
   tof->SetBranchAddress("avg_lx" ,&tof_lx);
   tof->SetBranchAddress("avg_ly" ,&tof_ly);
   tof->SetBranchAddress("avg_lz" ,&tof_lz); 
   tof->SetBranchAddress("vx"     ,&tof_vx);
   tof->SetBranchAddress("vy"     ,&tof_vy);
   tof->SetBranchAddress("vz"     ,&tof_vz);
//   tof->SetBranchAddress("mpid"   ,tof_mpid);
//   tof->SetBranchAddress("mvx"    ,tof_mvx);
//   tof->SetBranchAddress("mvy"    ,tof_mvy);
//   tof->SetBranchAddress("mvz"    ,tof_mvz);


   Long64_t nentries = tof->GetEntries();
   cout << "N. entries:" << nentries << " " << tof->GetEntries() << endl;



// Create histos
   TH2F *hi_tof_occ = new TH2F("hi_tof_occ", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ->GetYaxis()->SetTitle("Sector");
   hi_tof_occ->GetZaxis()->SetTitle("MHz");

   TH2F *hi_tof_occ_norm = new TH2F("hi_tof_occ_norm", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ_norm->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ_norm->GetYaxis()->SetTitle("Sector");
   hi_tof_occ_norm->GetZaxis()->SetTitle("MHz/cm");

   TH2F *hi_tof_occ_5 = new TH2F("hi_tof_occ_5", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ_5->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ_5->GetYaxis()->SetTitle("Sector");
   hi_tof_occ_5->GetZaxis()->SetTitle("MHz");
    
   TH2F *hi_tof_occ_norm_5 = new TH2F("hi_tof_occ_norm_5", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ_norm_5->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ_norm_5->GetYaxis()->SetTitle("Sector");
   hi_tof_occ_norm_5->GetZaxis()->SetTitle("MHz/cm");
    

   TH2F *hi_tof_occ_5charged = new TH2F("hi_tof_occ_5charged", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ_5charged->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ_5charged->GetYaxis()->SetTitle("Sector");

   TH1F *hi_tof_pid = new TH1F("hi_tof_pid", "",3000, -600.,2400.);
   hi_tof_pid->GetXaxis()->SetTitle("PID");
   hi_tof_pid->GetYaxis()->SetTitle("Rate (MHz)");

   TH1F *hi_tof_pid_5 = new TH1F("hi_tof_pid_5", "",3000, -600.,2400.);
   hi_tof_pid_5->GetXaxis()->SetTitle("PID");
   hi_tof_pid_5->GetYaxis()->SetTitle("Rate (MHz)");

   TH1F *hi_tof1a_pid = new TH1F("hi_tof1a_pid", "",200, -600.,2400.);
   hi_tof1a_pid->GetXaxis()->SetTitle("PID");
   hi_tof1a_pid->GetYaxis()->SetTitle("Rate (MHz)");

   TH1F *hi_tof1b_pid = new TH1F("hi_tof1b_pid", "",200, -600.,2400.);
   hi_tof1b_pid->GetXaxis()->SetTitle("PID");
   hi_tof1b_pid->GetYaxis()->SetTitle("Rate (MHz)");

   TH1F *hi_tof2_pid = new TH1F("hi_tof2_pid", "",200, -600.,2400.);
   hi_tof2_pid->GetXaxis()->SetTitle("PID");
   hi_tof2_pid->GetYaxis()->SetTitle("Rate (MHz)");

   TH1F *hi_tof1a_pid_5 = new TH1F("hi_tof1a_pid_5", "",200, -600.,2400.);
   hi_tof1a_pid_5->GetXaxis()->SetTitle("PID");
   hi_tof1a_pid_5->GetYaxis()->SetTitle("Rate (MHz)");

   TH1F *hi_tof1b_pid_5 = new TH1F("hi_tof1b_pid_5", "",200, -600.,2400.);
   hi_tof1b_pid_5->GetXaxis()->SetTitle("PID");
   hi_tof1b_pid_5->GetYaxis()->SetTitle("Rate (MHz)");

   TH1F *hi_tof2_pid_5 = new TH1F("hi_tof2_pid_5", "",200, -600.,2400.);
   hi_tof2_pid_5->GetXaxis()->SetTitle("PID");
   hi_tof2_pid_5->GetYaxis()->SetTitle("Rate (MHz)");

   TH1F *hi_tof_vz_all = new TH1F("hi_tof_vz_all", "",200, 0.,800.);
   hi_tof_vz_all->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof_vz_all->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof_vz_e = new TH1F("hi_tof_vz_e", "",200, 0.,800.);
   hi_tof_vz_e->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof_vz_e->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof_vz_g = new TH1F("hi_tof_vz_g", "",200, 0.,800.);
   hi_tof_vz_g->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof_vz_g->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof_vz_h = new TH1F("hi_tof_vz_h", "",200, 0.,800.);
   hi_tof_vz_h->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof_vz_h->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof_vz_n = new TH1F("hi_tof_vz_n", "",200, 0.,800.);
   hi_tof_vz_n->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof_vz_n->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH2F *hi_tof_origin_all = new TH2F("hi_tof_origin_all", "",200, 0.,800.,200, 0.,100.);
   hi_tof_origin_all->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof_origin_all->GetYaxis()->SetTitle("R_{vertex} (cm)");
    
     TH2F *hi_tof_origin_all_temp = new TH2F("hi_tof_origin_all_temp", "",200, 0.,800.,200, 0.,100.);

  TH2F *hi_tof_origin_all_ene = new TH2F("hi_tof_origin_all_ene", "",200, 0.,800.,200, 0.,100.);
  hi_tof_origin_all_ene->GetXaxis()->SetTitle("Z_{vertex} (cm)");
  hi_tof_origin_all_ene->GetYaxis()->SetTitle("R_{vertex} (cm)");
  hi_tof_origin_all_ene->GetZaxis()->SetTitle("E (MeV)");
  
  TH2F *hi_tof_origin_all_ene_temp = new TH2F("hi_tof_origin_all_ene_temp", "",200, 0.,800.,200, 0.,100.);

   TH1F *hi_tof1a_vz = new TH1F("hi_tof1a_vz", "",200, 0.,800.);
   hi_tof1a_vz->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof1a_vz->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof1b_vz = new TH1F("hi_tof1b_vz", "",200, 0.,800.);
   hi_tof1b_vz->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof1b_vz->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof2_vz = new TH1F("hi_tof2_vz", "",200, 0.,800.);
   hi_tof2_vz->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof2_vz->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof1a_vz_5 = new TH1F("hi_tof1a_vz_5", "",200, 0.,800.);
   hi_tof1a_vz_5->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof1a_vz_5->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof1b_vz_5 = new TH1F("hi_tof1b_vz_5", "",200, 0.,800.);
   hi_tof1b_vz_5->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof1b_vz_5->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof2_vz_5 = new TH1F("hi_tof2_vz_5", "",200, 0.,800.);
   hi_tof2_vz_5->GetXaxis()->SetTitle("Z_{vertex} (cm)");
   hi_tof2_vz_5->GetYaxis()->SetTitle("Rate (MHz/4cm)");

   TH1F *hi_tof1a_edep = new TH1F("hi_tof1a_edep", "",200, 0.,10.);
   hi_tof1a_edep->GetXaxis()->SetTitle("E(MeV)");
   hi_tof1a_edep->GetYaxis()->SetTitle("Rate (MHz/50keV)");

   TH1F *hi_tof1b_edep = new TH1F("hi_tof1b_edep", "",200, 0.,10.);
   hi_tof1b_edep->GetXaxis()->SetTitle("E(MeV)");
   hi_tof1b_edep->GetYaxis()->SetTitle("Rate (MHz/50keV)");

   TH1F *hi_tof2_edep = new TH1F("hi_tof2_edep", "",200, 0.,10.);
   hi_tof2_edep->GetXaxis()->SetTitle("E(MeV)");
   hi_tof2_edep->GetYaxis()->SetTitle("Rate (MHz/50keV)");

   TH1F *hi_tof1a_edep_5 = new TH1F("hi_tof1a_edep_5", "",200, 0.,10.);
   hi_tof1a_edep_5->GetXaxis()->SetTitle("E(MeV)");
   hi_tof1a_edep_5->GetYaxis()->SetTitle("Rate (MHz/50keV)");
   hi_tof1a_edep_5->SetLineColor(2);

   TH1F *hi_tof1b_edep_5 = new TH1F("hi_tof1b_edep_5", "",200, 0.,10.);
   hi_tof1b_edep_5->GetXaxis()->SetTitle("E(MeV)");
   hi_tof1b_edep_5->GetYaxis()->SetTitle("Rate (MHz/50keV)");
   hi_tof1b_edep_5->SetLineColor(2);

   TH1F *hi_tof2_edep_5 = new TH1F("hi_tof2_edep_5", "",200, 0.,10.);
   hi_tof2_edep_5->GetXaxis()->SetTitle("E(MeV)");
   hi_tof2_edep_5->GetYaxis()->SetTitle("Rate (MHz/50keV)");
   hi_tof2_edep_5->SetLineColor(2);

   TH2F *hi_tof_occ_edep = new TH2F("hi_tof_occ_edep (MeV/us)", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ_edep->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ_edep->GetYaxis()->SetTitle("Sector");

   TH2F *hi_tof_occ_edep_norm = new TH2F("hi_tof_occ_edep_norm (MeV/us)", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ_edep_norm->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ_edep_norm->GetYaxis()->SetTitle("Sector");
   
   TH2F *hi_tof_occ_curr = new TH2F("hi_tof_occ_curr (uA)", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ_curr->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ_curr->GetYaxis()->SetTitle("Sector");
   
   TH2F *hi_tof_occ_adc = new TH2F("hi_tof_occ_adc (uA)", "",120, 0.,120.,6,1.,7.);
   hi_tof_occ_adc->GetXaxis()->SetTitle("Paddle");
   hi_tof_occ_adc->GetYaxis()->SetTitle("Sector");
   
   TH1F *hi_ntof = new TH1F("hi_ntof", "",10, 0.,10.);
   hi_ntof->GetXaxis()->SetTitle("Nhits");
   hi_ntof->GetYaxis()->SetTitle("Counts");





// in standard run mode (1 electron at a time)
//   float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (50*2300 or 118600 electrons in 250ns window)
   float levents=124000;
   int   nsum=124000./levents;
   float norm=124000/levents/(nentries);
   float lumi=nentries*250/10*levents/124000; // (mbarn^-1 or 10^27 cm^-2)
   float time=250/norm/1000;

   cout << "number of events to integrate = " << nsum  << endl;
   cout << "normalization factor          = " << norm  << endl;
   cout << "Run time                      = " << time  << " us" << endl;

   // initializing matrix
   for(int l=0; l<64; l++) ntof[l]=0;
   for(int i=0; i<3; i++) {
     for(int j=0; j<6; j++) {
         for(int k=0; k<60; k++) {
             for(int l=0; l<64; l++) {
                 hitS[i][j][k][l]=0;
                 edepS[i][j][k][l]=0;
                 timeS[i][j][k][l]=0;
             }
         }
      }
   }

   int   nfull=0;
   int   ntrigger=0;
   int   trigger_window=64; // coincidence time window used for TOF; has to be a multiple of 4 ns
   float Ethr=1;            // energy thresholds in MeV
   int nint=0;
   int ngoodentries=0;
   for (Long64_t jentry=0; jentry < nentries; jentry++) {
  
       tof_sector->clear();
       tof_paddle->clear();
       tof_ADCL->clear();
       tof_ADCR->clear();
       tof_TDCL->clear();
       tof_TDCR->clear();
       tof_E->clear();
       tof_Edep->clear();
       tof_t->clear();
       tof_pid->clear();
       tof_x->clear();
       tof_y->clear();
       tof_z->clear();
       tof_lx->clear();
       tof_ly->clear();
       tof_lz->clear(); 
       tof_vx->clear();
       tof_vy->clear();
       tof_vz->clear();
     
       int nb = g->GetEntry(jentry); 
       tof->GetEntry(jentry); 
       ngoodentries++;
       if(int(jentry/1000)*1000==jentry) cout << "Analyzed " << jentry << " events of " << nentries << endl;
	
       ntofhit=tof_pid->size();

       for(int i=0; i<ntofhit; i++) {
           
           if((*tof_pid)[i]==2112 || (*tof_pid)[i]==2212) {
               mass=938;
           }
           else{
               mass=0;
           }
           
           
	 int offset=0;
	 if      ((*tof_layer)[i]==1) offset=70;
	 else if ((*tof_layer)[i]==3) offset=100;

	 // Q L/R = ADC L/R *( adode/dynode) * fADC_conv * dT * 50 Ohm
	 double weightADC =  (*tof_ADCL)[i]*3.*0.242E-3*4E-3/50;
	 // numbers for conversions from energy in MeV to Charge in uC are taken from Dan Carman's note
	 double weightCurrent = (*tof_Edep)[i]*373.*1.6E-7/10.;	 
	 if((*tof_layer)[i]==2) weightCurrent = (*tof_Edep)[i]*1158.*1.6E-7/12.;
	 double width = 15;
	 if     ((*tof_layer)[i]==2) width = 6;
	 else if((*tof_layer)[i]==3) width = 22;

	 bool charged = !(sqrt((*tof_vx)[i]*(*tof_vx)[i]/100.+(*tof_vy)[i]*(*tof_vy)[i]/100.)<100. && (*tof_vz)[i]/10.<150. 
		      && ((*tof_pid)[i]==22 || (*tof_pid)[i]==2112));
	 hi_tof_occ->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.);
	 hi_tof_occ_norm->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.,1./width);
	 hi_tof_occ_edep->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.,(*tof_Edep)[i]);
	 hi_tof_occ_edep_norm->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.,(*tof_Edep)[i]/width);
	 hi_tof_occ_curr->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.,weightCurrent); 
	 hi_tof_occ_adc->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.,weightADC);
	 hi_tof_pid->Fill((*tof_pid)[i]*1.);
	 hi_tof_vz_all->Fill((*tof_vz)[i]/10.);
	 if(abs((*tof_pid)[i])==11) hi_tof_vz_e->Fill((*tof_vz)[i]/10.);
	 else if((*tof_pid)[i]==22) hi_tof_vz_g->Fill((*tof_vz)[i]/10.);
	 else                       hi_tof_vz_h->Fill((*tof_vz)[i]/10.);
	 if((*tof_pid)[i]==2112)    hi_tof_vz_n->Fill((*tof_vz)[i]/10.);
	 hi_tof_origin_all->Fill((*tof_vz)[i]/10.,sqrt((*tof_vx)[i]*(*tof_vx)[i]/100.+(*tof_vy)[i]*(*tof_vy)[i]/100.));
     hi_tof_origin_all_temp->Fill((*tof_vz)[i]/10.,sqrt((*tof_vx)[i]*(*tof_vx)[i]/100.+(*tof_vy)[i]*(*tof_vy)[i]/100.));
    if((*tof_E)[i]>mass) hi_tof_origin_all_ene_temp->Fill((*tof_vz)[i]/10.,sqrt((*tof_vx)[i]*(*tof_vx)[i]/100.+(*tof_vy)[i]*(*tof_vy)[i]/100.), sqrt((*tof_E)[i]*(*tof_E)[i]-mass*mass) );
	 if((*tof_Edep)[i]>Ethr) {
	   hi_tof_occ_5->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.);
       hi_tof_occ_norm_5->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.,1./width);
	   if(charged) hi_tof_occ_5charged->Fill((*tof_paddle)[i]*1.+offset,(*tof_sector)[i]*1.);
	   hi_tof_pid_5->Fill((*tof_pid)[i]*1.);
	 }
	 if((*tof_layer)[i]==1) {
	   hi_tof1a_pid->Fill((*tof_pid)[i]*1.);
	   hi_tof1a_vz->Fill((*tof_vz)[i]/10.);
	   hi_tof1a_edep->Fill((*tof_Edep)[i]);
	   if((*tof_Edep)[i]>Ethr) {
	     hi_tof1a_pid_5->Fill((*tof_pid)[i]*1.);
	     hi_tof1a_vz_5->Fill((*tof_vz)[i]/10.);
	     hi_tof1a_edep_5->Fill((*tof_Edep)[i]);
	   }
	 }
	 else if((*tof_layer)[i]==2) {
	   hi_tof1b_pid->Fill((*tof_pid)[i]*1.);
	   hi_tof1b_vz->Fill((*tof_vz)[i]/10.);
           hi_tof1b_edep->Fill((*tof_Edep)[i]);
	   if((*tof_Edep)[i]>Ethr) {
	     hi_tof1b_pid_5->Fill((*tof_pid)[i]*1.);
	     hi_tof1b_vz_5->Fill((*tof_vz)[i]/10.);
	     hi_tof1b_edep_5->Fill((*tof_Edep)[i]);
	   }
	 }
	 else if((*tof_layer)[i]==3) {
	   hi_tof2_pid->Fill((*tof_pid)[i]*1.);
	   hi_tof2_vz->Fill((*tof_vz)[i]/10.);
           hi_tof2_edep->Fill((*tof_Edep)[i]);
	   if((*tof_Edep)[i]>Ethr) {
	     hi_tof2_pid_5->Fill((*tof_pid)[i]*1.);
	     hi_tof2_vz_5->Fill((*tof_vz)[i]/10.);
	     hi_tof2_edep_5->Fill((*tof_Edep)[i]);
	   }
	 }
       }

//	cout << nint << " " << nsum << endl;
       if(nint<nsum){
          for(int i=0; i<ntofhit; i++){
	     int ix=(*tof_layer)[i]-1;
	     int iy=(*tof_sector)[i]-1;
	     int ik=(*tof_paddle)[i]-1;
	     int il=(int) (*tof_t)[i]/4;
	     if(il<64 && (*tof_t)[i]>0 && (*tof_t)[i]<256 && !(sqrt((*tof_vx)[i]*(*tof_vx)[i]/100.+(*tof_vx)[i]*(*tof_vx)[i]/100.)<100. && (*tof_vz)[i]/10.<150. && ((*tof_pid)[i]==22 || (*tof_pid)[i]==2112))) {
//	     if(il<64 && (*tof_t)[i]>0 && (*tof_t)[i]<256 && (sqrt((*tof_vx)[i]*(*tof_vx)[i]/100.+(*tof_vx)[i]*(*tof_vx)[i]/100.)<100. && (*tof_vz)[i]/10.<150. && (abs((*tof_pid)[i])==11 || abs((*tof_pid)[i])==211 || (*tof_pid)[i]==2212)))  {
//	     if(il<64 && (*tof_t)[i]>0 && (*tof_t)[i]<256 && (sqrt((*tof_vx)[i]*(*tof_vx)[i]/100.+(*tof_vx)[i]*(*tof_vx)[i]/100.)<100. && (*tof_vz)[i]/10.<150. && (abs((*tof_pid)[i])==11)))  {
//	     if(il<64 && (*tof_t)[i]>0 && (*tof_t)[i]<256) {
	        hitS[ix][iy][ik][il]++;
	        edepS[ix][iy][ik][il]=edepS[ix][iy][ik][il]+(*tof_Edep)[i];
	        timeS[ix][iy][ik][il]=timeS[ix][iy][ik][il]+(*tof_t)[i];
	     }
          }
	  nint++;
       }
       if(nint==nsum) {
	 for(int i=0; i<3; i++) {
	   for(int j=0; j<6; j++) {
	     for(int k=0; k<60; k++) {
	       for(int l=0; l<64; l++) {
		 if(hitS[i][j][k][l]>0) {
		   timeS[i][j][k][l]/=hitS[i][j][k][l];
		   if(timeS[i][j][k][l]<l*4 || timeS[i][j][k][l]>(l+1)*4) 
		     cout << "warning: error in time clustering  t=" << timeS[i][j][k][l] << " i,y,k,l=" << i << " " << j << " " << k << " " << l << " " << endl;
		 }
	       }
	     }
	   }
	 }
	 int it=0;
	//        while(it<64-(trigger_window/4)) {
	 while(it<1) {
	   for(int i=1; i<3; i++) {
	     for(int j=0; j<6; j++) {
	       for(int k=0; k<60; k++) {
		 //                        for(int l=it; l<it+(trigger_window/4); l++) {
		 for(int l=8; l<8+(trigger_window/4); l++) {
		   if(edepS[i][j][k][l]>5.) ntof[it]++;
		 }
	       }
	     }
	   }
	   hi_ntof->Fill(ntof[it]);
	   if(ntof[it]>=3) {
	     ntrigger++;
	     it+=(trigger_window/4);
	   }
	   else {
	     it++;
	   }
	 }
	 for(int i=0; i<3; i++) {
	   for(int j=0; j<6; j++) {
	     for(int k=0; k<60; k++) {
	       for(int l=0; l<64; l++) {
		 hitS[i][j][k][l]=0;
		 edepS[i][j][k][l]=0;
		 timeS[i][j][k][l]=0;
	       }
	     }
	   }
	 }
	 for(int l=0; l<64; l++) ntof[l]=0;
	 nint=0;
	 nfull++;
       }
   }
   
// in standard run mode (1 electron at a time)
//   float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (50*2300 or 118600 electrons in 250ns window)
   norm=124000/levents/ngoodentries;
   lumi=ngoodentries*250/10*levents/124000; // (mbarn^-1 or 10^27 cm^-2)
   time=250/norm/1000;

   cout << "number of events to integrate = " << nsum  << endl;
   cout << "normalization factor          = " << norm  << endl;
   cout << "Run time                      = " << time  << " us" << endl;
   
   cout << ngoodentries << " events analyzed" << endl;
   cout << ntrigger << " trigger found" << endl;
   cout << "Trigger rate: " << ntrigger*1./(ngoodentries*trigger_window/1000.) << " MHz" << endl;

   // normalizing rate histogram to 10^35 luminosity
   hi_tof_occ->Scale(1/time);
   hi_tof_occ_norm->Scale(1/time);
   hi_tof_occ_edep->Scale(1/time);
   hi_tof_occ_edep_norm->Scale(1/time);
   hi_tof_occ_curr->Scale(1E6/time);
   hi_tof_occ_adc->Scale(1E6/time);
   hi_tof_occ_5->Scale(1/time);
   hi_tof_occ_norm_5->Scale(1/time);
   hi_tof_occ_5charged->Scale(1/time);
   hi_tof_pid->Scale(1/time);
   hi_tof_pid_5->Scale(1/time);
   hi_tof1a_pid->Scale(1/time);
   hi_tof1b_pid->Scale(1/time);
   hi_tof2_pid->Scale(1/time);
   hi_tof1a_pid_5->Scale(1/time);
   hi_tof1b_pid_5->Scale(1/time);
   hi_tof2_pid_5->Scale(1/time);
   hi_tof_vz_all->Scale(1/time);
   hi_tof_vz_e->Scale(1/time);
   hi_tof_vz_g->Scale(1/time);
   hi_tof_vz_h->Scale(1/time);
   hi_tof_vz_n->Scale(1/time);
   hi_tof_origin_all->Scale(1/time);
   hi_tof1a_vz->Scale(1/time);
   hi_tof1b_vz->Scale(1/time);
   hi_tof2_vz->Scale(1/time);
   hi_tof1a_vz_5->Scale(1/time);
   hi_tof1b_vz_5->Scale(1/time);
   hi_tof2_vz_5->Scale(1/time);
   hi_tof1a_edep->Scale(1/time);
   hi_tof1b_edep->Scale(1/time);
   hi_tof2_edep->Scale(1/time);
   hi_tof1a_edep_5->Scale(1/time);
   hi_tof1b_edep_5->Scale(1/time);
   hi_tof2_edep_5->Scale(1/time);
    
    hi_tof_origin_all_ene->Divide(hi_tof_origin_all_ene_temp, hi_tof_origin_all_temp);
   

   TCanvas *c_occ=new TCanvas("c_occ","Occupancy",750,1000);
   c_occ->Divide(1,2);
   c_occ->cd(1);
   hi_tof_occ->Draw("COLZ");
   c_occ->cd(2);
   hi_tof_occ_norm->Draw("COLZ");
   c_occ->Print("tof_occupancy.pdf");

   TCanvas *c_occ_cut=new TCanvas("c_occ_cut","Occupancy_Cuts",750,1000);
   c_occ_cut->Divide(1,2);
   c_occ_cut->cd(1);
   hi_tof_occ_5->Draw("COLZ");
   c_occ_cut->cd(2);
  // hi_tof_occ_5charged->Draw("COLZ");
    hi_tof_occ_norm_5->Draw("COLZ");
   c_occ_cut->Print("tof_occupancy_cut.pdf");

   TCanvas *c_pid0=new TCanvas("c_pid0","PID0",750,1000);
   c_pid0->Divide(1,2);
   c_pid0->cd(1);
   gPad->SetLogy();
   hi_tof_pid->Draw("H");
   c_pid0->cd(2);
   gPad->SetLogy();
   hi_tof_pid_5->Draw("H");
   c_pid0->Print("tof_pid0.pdf");


   TCanvas *c_pid=new TCanvas("c_pid","PID",750,1000);
   c_pid->Divide(2,3);
   c_pid->cd(1);
   gPad->SetLogy();
   hi_tof1a_pid->Draw("H");
   c_pid->cd(2);
   gPad->SetLogy();
   hi_tof1a_pid_5->Draw("H");
   c_pid->cd(3);
   gPad->SetLogy();
   hi_tof1b_pid->Draw("H");
   c_pid->cd(4);
   gPad->SetLogy();
   hi_tof1b_pid_5->Draw("H");
   c_pid->cd(5);
   gPad->SetLogy();
   hi_tof2_pid->Draw("H");
   c_pid->cd(6);
   gPad->SetLogy();
   hi_tof2_pid_5->Draw("H");
   c_pid->Print("tof_pid.pdf");

   TCanvas *c_origin=new TCanvas("c_origin","Origin",750,1000);
   c_origin->Divide(1,2);
   c_origin->cd(1);
   gPad->SetLogz();
   hi_tof_origin_all->Draw("COLZ");
   c_origin->cd(2);
   gPad->SetLogy();
   hi_tof_vz_all->SetMinimum(0.001);
   hi_tof_vz_all->Draw("H");
   hi_tof_vz_e->SetLineColor(2);
   hi_tof_vz_e->Draw("SAME");
   hi_tof_vz_g->SetLineColor(4);
   hi_tof_vz_g->Draw("SAME");
   hi_tof_vz_h->SetLineColor(kGreen);
   hi_tof_vz_h->Draw("SAME");
   hi_tof_vz_n->SetLineColor(kGreen+2);
   hi_tof_vz_n->Draw("SAME");
    
    TLegend *leg1 = new TLegend(0.7,0.75,0.96,0.96);
     leg1->SetTextSize(.04);
     leg1->AddEntry(hi_tof_vz_all,"All","l");
     leg1->AddEntry(hi_tof_vz_e,"electrons","l");
     leg1->AddEntry(hi_tof_vz_g,"photons","l");
     leg1->AddEntry(hi_tof_vz_n,"neutron","l");
     leg1->AddEntry(hi_tof_vz_h,"hadron","l");
     leg1->Draw();
    
   c_origin->Print("tof_origin.pdf");

    TCanvas *c_origin_e=new TCanvas("c_origin_e","Origin_e",750,1000);
    c_origin_e->Divide(1,2);
    c_origin_e->cd(1);
    gPad->SetLogz();
    hi_tof_origin_all->SetMinimum(1E-3);
    hi_tof_origin_all->Draw("COLZ");
    c_origin_e->cd(2);
    gPad->SetLogz();
     hi_tof_origin_all_ene->Draw("COLZ");
     c_origin_e->Print("tof_origin_ene.pdf");
    
   TCanvas *c_vz=new TCanvas("c_vz","VZ",750,1000);
   c_vz->Divide(2,3);
   c_vz->cd(1);
   hi_tof1a_vz->Draw();
   c_vz->cd(2);
   hi_tof1a_vz_5->Draw();
   c_vz->cd(3);
   hi_tof1b_vz->Draw();
   c_vz->cd(4);
   hi_tof1b_vz_5->Draw();
   c_vz->cd(5);
   hi_tof2_vz->Draw();
   c_vz->cd(6);
   hi_tof2_vz_5->Draw();
   c_vz->Print("tof_vertex.pdf");

   TCanvas *c_edep=new TCanvas("c_edep","EDep",750,1000);
   c_edep->Divide(1,3);
   c_edep->cd(1);
   gPad->SetLogy();
   hi_tof1a_edep->Draw();
   hi_tof1a_edep_5->Draw("SAME");
   c_edep->cd(2);
   gPad->SetLogy();
   hi_tof1b_edep->Draw();
   hi_tof1b_edep_5->Draw("SAME");
   c_edep->cd(3);
   gPad->SetLogy();
   hi_tof2_edep->Draw();
   hi_tof2_edep_5->Draw("SAME");
   c_edep->Print("tof_edep.pdf");

   TCanvas *c_occ_edep=new TCanvas("c_occ_edep","Occ_edep",750,1000);
   c_occ_edep->Divide(1,2);
   c_occ_edep->cd(1);
   hi_tof_occ_edep->Draw("COLZ");
   c_occ_edep->cd(2);
   hi_tof_occ_edep_norm->Draw("COLZ");
   c_occ_edep->Print("tof_occ_edep.pdf");
   
   TCanvas *c_occ_curr=new TCanvas("c_occ_curr","Occ_curr",750,1000);
   c_occ_curr->Divide(1,2);
   c_occ_curr->cd(1);
   hi_tof_occ_curr->Draw("COLZ");
   c_occ_curr->cd(2);
   hi_tof_occ_adc->Draw("COLZ");
   c_occ_curr->Print("tof_currents.pdf");
   
   TCanvas *c3=new TCanvas("c3","HITS",250,250);
   hi_ntof->Draw();
   c3->Print("tof_hits.pdf)");

   TCanvas *c_occ_kpp=new TCanvas("c_occ_kpp","Occ_kpp",750,1000);
   c_occ_kpp->Divide(1,2);
   c_occ_kpp->cd(1);
   TH2F *hi_tof_occ_kHz=(TH2F*)hi_tof_occ_5->Clone();
   hi_tof_occ_kHz->Scale(1000.);
   hi_tof_occ_kHz->Draw("COLZ");
   c_occ_kpp->cd(2);
   hi_tof_occ_adc->Draw("COLZ");
   c_occ_kpp->Print("tof_kpp.pdf");
   
   gui.Run(1);

}   


