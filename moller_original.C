
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
#include <TLine.h>
#include "TFile.h"
#include "TApplication.h"
#include <TROOT.h>

using namespace std; 
TApplication gui("GUI",0,NULL);

#define NPLANES 40
#define NGOOD 16
#define NMIN  15
#define NMAX  27
#define NX 4
#define NY 3
#define NLUMI 15625

//float e_moller(float,float);
//float moller_xsec(float,float);
float e_moller(float ebeam, float th)
{
	float me=0.511/1000;
	float re=2.817;          // x 10^-15 m
	float degrad=3.14159/180.;
	// mandelstan s
	float s=2*me*(ebeam+me);
	// gamma factor
	float gcm=(ebeam+me)/sqrt(s);	
	// symmetric pair angle (corresponds to th_cm=90 deg)
	float th_pair=atan(2*me/(ebeam+me))/degrad;
	// cm angle in radians
	float th_cm=atan(gcm*tan(th*degrad))*2; 
	
	float e_moller=(ebeam+me)/2+(ebeam-me)*cos(th_cm)/2;
	return e_moller;
}


float moller_xsec(float ebeam, float th)
{
	float me=0.511/1000;
	float re=2.817;          // x 10^-15 m
	float degrad=3.14159/180.;
	// mandelstan s
	float s=2*me*(ebeam+me);
	// gamma factor
	float gcm=(ebeam+me)/sqrt(s);
	// cm angle in radians
	float th_cm=atan(gcm*tan(th*degrad))*2; // cm angle in radians
	// dsigma/dth_cm
	float xsec_cm=(re*re/100)*(me*me/s)*pow((3+pow(cos(th_cm),2)),2)/pow(sin(th_cm),3)*2*3.14159; //  in barn/rad
	float moller_xsec=gcm*(1+cos(th_cm))/pow(cos(th*degrad),2)*xsec_cm*1000;
	return moller_xsec;
}


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
   gStyle->SetPalette(1);
   gStyle->SetOptFit(111);
   gStyle->SetOptStat("nemri");
//   gStyle->SetOptStat("");

  
   TChain *g= new TChain("generated");
   g->Add("out*.root");

   TChain *f = new TChain("flux");
   f->Add("out*.root");



   // header  & generated
   vector<int>   *evn = new vector<int>;
   vector<double> *px = new vector<double>;
   vector<double> *py = new vector<double>;
   vector<double> *pz = new vector<double>;
   double         p,theta,phi;

   //flux
   vector<double>   *flux_ID = new vector<double>;
   vector<double>  *flux_pid = new vector<double>;
   vector<double> *flux_mpid = new vector<double>;
   vector<double> *flux_edep = new vector<double>;
   vector<double>    *flux_E = new vector<double>;
   vector<double>    *flux_t = new vector<double>;
   vector<double>    *flux_x = new vector<double>;
   vector<double>    *flux_y = new vector<double>;
   vector<double>    *flux_z = new vector<double>;
   vector<double>   *flux_px = new vector<double>;
   vector<double>   *flux_py = new vector<double>;
   vector<double>   *flux_pz = new vector<double>;
   vector<double>   *flux_vx = new vector<double>;
   vector<double>   *flux_vy = new vector<double>;
   vector<double>   *flux_vz = new vector<double>;
   vector<double>  *flux_mvx = new vector<double>;
   vector<double>  *flux_mvy = new vector<double>;
   vector<double>  *flux_mvz = new vector<double>;
  
   // generated
   g->SetBranchAddress("px",&px);
   g->SetBranchAddress("py",&py);
   g->SetBranchAddress("pz",&pz);
 
   // flux
   f->SetBranchAddress("id",&flux_ID);
   f->SetBranchAddress("trackE",&flux_E);
   f->SetBranchAddress("totEdep",&flux_edep);
   f->SetBranchAddress("avg_t",&flux_t);
   f->SetBranchAddress("pid",&flux_pid);
   f->SetBranchAddress("mpid",&flux_mpid);
   f->SetBranchAddress("mvx",&flux_mvx);
   f->SetBranchAddress("mvy",&flux_mvy);
   f->SetBranchAddress("mvz",&flux_mvz);
   f->SetBranchAddress("px",&flux_px);
   f->SetBranchAddress("py",&flux_py);
   f->SetBranchAddress("pz",&flux_pz);
   f->SetBranchAddress("vx",&flux_vx);
   f->SetBranchAddress("vy",&flux_vy);
   f->SetBranchAddress("vz",&flux_vz);
   f->SetBranchAddress("avg_x",&flux_x);
   f->SetBranchAddress("avg_y",&flux_y);
   f->SetBranchAddress("avg_z",&flux_z);

   // Variables
   double rate_25_50[NPLANES]={0};
   double rate_50_350[NPLANES]={0};
   double rate_25_50_emin[NPLANES]={0};
   double rate_50_350_emin[NPLANES]={0};
   double spot_m_maxpos[NPLANES]={0};
   double spot_m_maxval[NPLANES]={0};
   double spot_m_size[NPLANES]={0};
   double spot_m_angle[NPLANES]={0};
   double rclear[NPLANES]={0};
   double zplane[NPLANES]={0};

   for(int i=0; i<NPLANES; i++) {
     zplane[i]=10+2*i;
     rclear[i]=zplane[i]*tan(5/57.296)-1.129;
   }

   
// Create histos
// Electrons
   TH2F *hi_elec_xy   = new TH2F("x vs. y for electrons", "",100, -50.,50., 100, -50.,50.);
   hi_elec_xy->GetXaxis()->SetTitle("x (cm)");
   hi_elec_xy->GetYaxis()->SetTitle("y (cm)");
   TH1F *hi_elec_e    = new TH1F("Energy of electrons", "", 100, 0., 11.);
   hi_elec_e->GetXaxis()->SetTitle("E (GeV)");
   hi_elec_e->GetYaxis()->SetTitle("Counts");
   TH1F *hi_elec_th   = new TH1F("Theta for electrons", "", 100, 0., 90.);
   hi_elec_th->GetXaxis()->SetTitle("#theta (deg)");
   hi_elec_th->GetYaxis()->SetTitle("Counts");
   TH1F *hi_elec_th5   = new TH1F("Theta for electrons (<5 deg)", "", 100, 0., 5.);
   hi_elec_th5->GetXaxis()->SetTitle("#theta (deg)");
   hi_elec_th5->GetYaxis()->SetTitle("Counts");
   TH2F *hi_elec_eth  = new TH2F("Energy vs. Theta for electrons", "",100, 0., 11., 100, 0., 11.);
   hi_elec_eth->GetXaxis()->SetTitle("#theta (deg)");
   hi_elec_eth->GetYaxis()->SetTitle("Energy (GeV)");

// Gammas
   TH2F *hi_first_xy   = new TH2F("x vs. y for the first particle", "",100, -50.,50., 100, -50.,50.);
   hi_first_xy->GetXaxis()->SetTitle("x (cm)");
   hi_first_xy->GetYaxis()->SetTitle("y (cm)");
   TH1F *hi_first_e    = new TH1F("Energy of the first particle", "", 100, 0., 11.);
   hi_first_e->GetXaxis()->SetTitle("E (GeV)");
   hi_first_e->GetYaxis()->SetTitle("Counts");
   TH1F *hi_first_th   = new TH1F("Theta for the first particle", "", 100, 0., 90.);
   hi_first_th->GetXaxis()->SetTitle("#theta (deg)");
   hi_first_th->GetYaxis()->SetTitle("Counts");

// Primaries
   TH2F *hi_prim_xy   = new TH2F("x vs. y for primaries", "",100, -50.,50., 100, -50.,50.);
   hi_prim_xy->GetXaxis()->SetTitle("x (cm)");
   hi_prim_xy->GetYaxis()->SetTitle("y (cm)");
   TH1F *hi_prim_e    = new TH1F("Energy of primaries", "", 100, 0., 11.);
   hi_prim_e->GetXaxis()->SetTitle("E (GeV)");
   hi_prim_e->GetYaxis()->SetTitle("Counts");
   TH1F *hi_prim_th   = new TH1F("Theta for primaries", "", 100, 0., 90.);
   hi_prim_th->GetXaxis()->SetTitle("#theta (deg)");
   hi_prim_th->GetYaxis()->SetTitle("Counts");

// Secondaries
   TH2F *hi_seco_xy   = new TH2F("x vs. y for secondaries", "",100, -50.,50., 100, -50.,50.);
   hi_seco_xy->GetXaxis()->SetTitle("x (cm)");
   hi_seco_xy->GetYaxis()->SetTitle("y (cm)");
   TH1F *hi_seco_e    = new TH1F("Energy of secondaries", "", 100, 0., 11.);
   hi_seco_e->GetXaxis()->SetTitle("E (GeV)");
   hi_seco_e->GetYaxis()->SetTitle("Counts");
   TH1F *hi_seco_th   = new TH1F("Theta for secondaries", "", 100, 0., 90.);
   hi_seco_th->GetXaxis()->SetTitle("#theta (deg)");
   hi_seco_th->GetYaxis()->SetTitle("Counts");
   TH1F *hi_seco_th5  = new TH1F("Theta for secondaries (<5 deg)", "", 100, 0., 5.);
   hi_seco_th5->GetXaxis()->SetTitle("#theta (deg)");
   hi_seco_th5->GetYaxis()->SetTitle("Counts");
   TH2F *hi_seco_eth  = new TH2F("Energy vs. Theta for secondaries", "",100, 0., 11., 100, 0., 5.);
   hi_seco_eth->GetXaxis()->SetTitle("#theta (deg)");
   hi_seco_eth->GetYaxis()->SetTitle("Energy (GeV)");

// Rate in FT range
   TH1F *hi_rate_th5  = new TH1F("Electron Rate (<5 deg)", "", 100, 0., 5.);
   hi_rate_th5->GetXaxis()->SetTitle("#theta (deg)");
   hi_rate_th5->GetYaxis()->SetTitle("Counts/s");
   TH1F *hi_rate_ft  = new TH1F("Electron Rate", "", 100, 0., 5.);
   hi_rate_ft->GetXaxis()->SetTitle("#theta (deg)");
   hi_rate_ft->GetYaxis()->SetTitle("Counts/s");
   
// Electron Spots
   TH1F *hi_spot[NPLANES];
   TH2F *hi_spot_2d[NPLANES];
   TH1F *hi_spot_e[NPLANES];
   TH1F *hi_spot_ene[NPLANES];
   TH2F *hi_spot_e_2d[NPLANES];
   TH1F *hi_spot_e_ene[NPLANES];
   TH1F *hi_spot_m[NPLANES];
   TH2F *hi_spot_m_2d[NPLANES];
   TH1F *hi_spot_m_ene[NPLANES];
   TH1F *hi_spot_m_gt200[NPLANES];
   TH1F *hi_spot_m_gt200_clone[NPLANES];
   TH1F *hi_spot_m_gt200_ene[NPLANES];
   TH1F *hi_spot_g[NPLANES];
   TH2F *hi_spot_g_2d[NPLANES];
   TH1F *hi_spot_g_ene[NPLANES];
   TH1F *hi_spot_o[NPLANES];
   TH2F *hi_spot_o_2d[NPLANES];
   TH1F *hi_spot_o_ene[NPLANES];
   TH1F *hi_spotz[NPLANES];

   char htitle[50]; // create a char of ample size
   for(int i=0; i<NPLANES; i++){
     //     zplane[i]=10.+2.*i;
     //     if(i==31) zplane[31]=179.4;
     sprintf (htitle,"Background Spot at %f cm",zplane[i]); 
     hi_spot[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot[i]->GetYaxis()->SetTitle("Counts/s");

     sprintf (htitle,"Background E vs. R at %f cm",zplane[i]); 
     hi_spot_2d[i] = new TH2F(htitle,"",100,0.,10.,100,0.,1000.);
     hi_spot_2d[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_2d[i]->GetYaxis()->SetTitle("E (MeV)");

     sprintf (htitle,"Background Energy at %f cm",zplane[i]); 
     hi_spot_ene[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_ene[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_ene[i]->GetYaxis()->SetTitle("Edep(MeV)/ns");

     sprintf (htitle,"Electron Spot at %f cm",zplane[i]); 
     hi_spot_e[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_e[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_e[i]->GetYaxis()->SetTitle("Counts/s");

     sprintf (htitle,"Electron E vs. R at %f cm",zplane[i]); 
     hi_spot_e_2d[i] = new TH2F(htitle,"",100,0.,10.,100,0.,1000.);
     hi_spot_e_2d[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_e_2d[i]->GetYaxis()->SetTitle("E (MeV)");

     sprintf (htitle,"Electron Energy at %f cm",zplane[i]); 
     hi_spot_e_ene[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_e_ene[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_e_ene[i]->GetYaxis()->SetTitle("Edep(MeV)/ns");

     sprintf (htitle,"Moller Spot at %f cm",zplane[i]); 
     hi_spot_m[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_m[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_m[i]->GetYaxis()->SetTitle("Counts/s");

     sprintf (htitle,"Moller E vs. R at %f cm",zplane[i]); 
     hi_spot_m_2d[i] = new TH2F(htitle,"",100,0.,10.,100,0.,1000.);
     hi_spot_m_2d[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_m_2d[i]->GetYaxis()->SetTitle("E (MeV)");

     sprintf (htitle,"Moller Energy at %f cm",zplane[i]); 
     hi_spot_m_ene[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_m_ene[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_m_ene[i]->GetYaxis()->SetTitle("Edep(MeV)/ns");

     sprintf (htitle,"Moller E> 200 MeV Spot at %f cm",zplane[i]); 
     hi_spot_m_gt200[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_m_gt200[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_m_gt200[i]->GetYaxis()->SetTitle("Counts/s");

     sprintf (htitle,"Moller E> 200 MeV Energy at %f cm",zplane[i]); 
     hi_spot_m_gt200_ene[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_m_gt200_ene[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_m_gt200_ene[i]->GetYaxis()->SetTitle("Edep(MeV)/ns");

    sprintf (htitle,"Gamma Spot at %f cm",zplane[i]); 
     hi_spot_g[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_g[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_g[i]->GetYaxis()->SetTitle("Counts/s");

     sprintf (htitle,"Gamma E vs. R at %f cm",zplane[i]); 
     hi_spot_g_2d[i] = new TH2F(htitle,"",100,0.,10.,100,0.,1000.);
     hi_spot_g_2d[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_g_2d[i]->GetYaxis()->SetTitle("E (MeV)");

     sprintf (htitle,"Gamma Energy at %f cm",zplane[i]); 
     hi_spot_g_ene[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_g_ene[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_g_ene[i]->GetYaxis()->SetTitle("Edep(MeV)/ns");

     sprintf (htitle,"Other Spot at %f cm",zplane[i]); 
     hi_spot_o[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_o[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_o[i]->GetYaxis()->SetTitle("Counts/s");

     sprintf (htitle,"Other E vs. R at %f cm",zplane[i]); 
     hi_spot_o_2d[i] = new TH2F(htitle,"",100,0.,10.,100,0.,1000.);
     hi_spot_o_2d[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_o_2d[i]->GetYaxis()->SetTitle("E (MeV)");

     sprintf (htitle,"Other Energy at %f cm",zplane[i]); 
     hi_spot_o_ene[i] = new TH1F(htitle,"",100,0.,10.);
     hi_spot_o_ene[i]->GetXaxis()->SetTitle("R(cm)");
     hi_spot_o_ene[i]->GetYaxis()->SetTitle("Edep(MeV)/ns");

     sprintf (htitle,"Background z at %f cm",zplane[i]); 
     hi_spotz[i] = new TH1F(htitle,"",100,0.,200.);
     hi_spotz[i]->GetXaxis()->SetTitle("z(cm)");
     hi_spotz[i]->GetYaxis()->SetTitle("Counts");
   }
   
   
// Neutrons
   TH1F *hi_n_th = new TH1F("Theta for neutrons", "", 100, 0., 90.);
   hi_n_th->GetXaxis()->SetTitle("#theta (deg)");
   hi_n_th->GetYaxis()->SetTitle("Counts");
   TH1F *hi_n_e = new TH1F("Energy for neutrons", "", 100, 0., 1.);
   hi_n_e->GetXaxis()->SetTitle("E (GeV)");
   hi_n_e->GetYaxis()->SetTitle("Counts");
   TH2F *hi_n_eth = new TH2F("Energy vs. Theta for neutrons", "", 100, 0., 90., 100, 0., 1.);
   hi_n_eth->GetXaxis()->SetTitle("#theta (deg)");
   hi_n_eth->GetYaxis()->SetTitle("E (GeV)");


   
   Long64_t nentries = g->GetEntries();
	
   cout << nentries << endl;
// in standard run mode (1 electron at a time)
// float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (118600 electrons in 250ns window)
   float levents=1;
   float norm=NLUMI/levents/nentries;
   float lumi=nentries*250/10*levents/124000; // (mbarn^-1 or 10^27 cm^-2)
   float time=250/norm;
   cout << "Run Luminosity = " << lumi << " mbarn^-1" << endl;
   cout << "Run time       = " << time << " ns      " << endl;

   

   
   for (Long64_t jentry=0; jentry < nentries; jentry++) {
  
     px->clear();
     py->clear();
     pz->clear();
     flux_ID->clear();
     flux_pid->clear();
     flux_mpid->clear();
     flux_E->clear();
     flux_edep->clear();
     flux_t->clear();
     flux_x->clear();
     flux_y->clear();
     flux_z->clear();
     flux_px->clear();
     flux_py->clear();
     flux_pz->clear();
     flux_vx->clear();
     flux_vy->clear();
     flux_vz->clear();
     flux_mvx->clear();
     flux_mvy->clear();
     flux_mvz->clear();
     if(int(jentry/124000)*124000==jentry) cout << "Analyzed " << jentry << " events of " << nentries << endl;
     //		g->GetEntry(jentry);
     f->GetEntry(jentry);

     // Flux
     // analyzing hits
     int nfhit=flux_pid->size();
     for(int i=0; i<nfhit; i++) {
         float theta=57.3*atan(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/(*flux_z)[i]);
	 float radius=sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i]);  
         if((*flux_ID)[i]==1) {
	    if((*flux_pid)[i]==11) {
	       hi_elec_xy->Fill((*flux_x)[i]/10.,(*flux_y)[i]/10.);
	       hi_elec_e->Fill((*flux_E)[i]/1000);
	       hi_elec_th->Fill(theta);
	       hi_elec_th5->Fill(theta);
	       hi_elec_eth->Fill(theta,(*flux_E)[i]/1000);
	       if((*flux_E)[i]>200. && (*flux_E)[i]<4000. && theta>2. && theta < 5.) hi_rate_th5->Fill(theta);
	       if((*flux_E)[i]>500. && (*flux_E)[i]<4000. && theta>2. && theta < 5.) hi_rate_ft->Fill(theta);
	    }		
	    if(i==0 && (*flux_ID)[0]==1 && (*flux_pid)[0]==11) {
	       hi_first_xy->Fill((*flux_x)[i]/10.,(*flux_y)[i]/10.);
	       hi_first_e->Fill((*flux_E)[i]/1000);
	       hi_first_th->Fill(theta);
	    }		
	    if((*flux_pid)[i]==11 && (*flux_mpid)[i]==0) {
	       hi_prim_xy->Fill((*flux_x)[i]/10.,(*flux_y)[i]/10.);
	       hi_prim_e->Fill((*flux_E)[i]/1000);
	       hi_prim_th->Fill(theta);
	    }	
	    if((*flux_pid)[i]==11 && (*flux_mpid)[i]==11) {
	       hi_seco_xy->Fill((*flux_x)[i]/10.,(*flux_y)[i]/10.);
	       hi_seco_e->Fill((*flux_E)[i]/1000);
	       hi_seco_th->Fill(theta);
	       hi_seco_th5->Fill(theta);			
	       hi_seco_eth->Fill(theta,(*flux_E)[i]/1000);
	    }
	    if((*flux_pid)[i]==2112) {
	       hi_n_th->Fill(theta);
	       hi_n_e->Fill(sqrt((*flux_E)[i]*(*flux_E)[i]-940.*940.)/1000.);
	       hi_n_eth->Fill(theta,sqrt((*flux_E)[i]*(*flux_E)[i]-940.*940.)/1000.);
	    }
	 }
	 if((*flux_ID)[i]<NPLANES+1) {
	   //	   cout << (*flux_z)[i] << " " << radius << " " << zplane[(int)(*flux_ID)[i]-1] << " " << theta << endl;
	    if(radius>zplane[(int)(*flux_ID)[i]-1]*10*tan(2.5/57.3)-2. && theta<5) {
	       rate_25_50[(int) (*flux_ID)[i]-1]++;
	       if((*flux_E)[i]>0.00005) rate_25_50_emin[(int) (*flux_ID)[i]-1]++;
	    }
	    if(theta>5.0 && theta<35.) {
	       rate_50_350[(int) (*flux_ID)[i]-1]++;
	       if((*flux_E)[i]>0.00005) rate_50_350_emin[(int) (*flux_ID)[i]-1]++;
	    }
	    hi_spot[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.);
	    hi_spot_2d[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
	    hi_spot_ene[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
	    hi_spotz[(int) (*flux_ID)[i]-1]->Fill((*flux_z)[i]/10.);
	    if((*flux_pid)[i]==11 && (*flux_mpid)[i]==11) {
		hi_spot_m[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.);
		hi_spot_m_2d[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
		hi_spot_m_ene[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
		if((*flux_E)[i]>200) {
			hi_spot_m_gt200[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.);
			hi_spot_m_gt200_ene[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
		}
	    }
	    if((*flux_pid)[i]==11) {
	       hi_spot_e[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.);
	       hi_spot_e_2d[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
	       hi_spot_e_ene[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
	    }
	    else if((*flux_pid)[i]==22) {
	       hi_spot_g[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.);
	       hi_spot_g_2d[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
	       hi_spot_g_ene[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
	    }
	    else {
	       hi_spot_o[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.);
	       hi_spot_o_2d[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
	       hi_spot_o_ene[(int) (*flux_ID)[i]-1]->Fill(sqrt((*flux_x)[i]*(*flux_x)[i]+(*flux_y)[i]*(*flux_y)[i])/10.,(*flux_E)[i]);
	    }
	 }	
     }
   }
   
   // normalizing rate histogram to 10^35 luminosity
   hi_rate_th5->Scale(100000000/lumi);
   for(int i=0; i<NPLANES; i++) {
	hi_spot[i]->Scale(1.E9/time);
	hi_spot_e[i]->Scale(1.E9/time);
	hi_spot_m[i]->Scale(1.E9/time);
	hi_spot_m_gt200[i]->Scale(1.E9/time);
	hi_spot_g[i]->Scale(1.E9/time);
	hi_spot_o[i]->Scale(1.E9/time);
	hi_spot_ene[i]->Scale(1./time);
	hi_spot_e_ene[i]->Scale(1./time);
	hi_spot_m_ene[i]->Scale(1./time);
	hi_spot_m_gt200_ene[i]->Scale(1./time);
	hi_spot_g_ene[i]->Scale(1./time);
	hi_spot_o_ene[i]->Scale(1./time);
	rate_25_50[i]/=time;
	rate_50_350[i]/=time;
	rate_25_50_emin[i]/=time;
	rate_50_350_emin[i]/=time;
  }
   
   const int nxsec=1000;
   float thv[nxsec],xv[nxsec],ev[nxsec],rv[nxsec];
   for(int i=0; i<nxsec; i++) {
	   thv[i]=(i+0.5)*(90./nxsec);
	   xv[i]=moller_xsec(11.,thv[i]);
	   ev[i]=e_moller(11.,thv[i]);
	   rv[i]=xv[i]*(3.14156*5./180./100)*100000000;
   }
   TGraph *gr_eth  = new TGraph(nxsec,thv,ev);
   TGraph *gr_xsec = new TGraph(nxsec,thv,xv);
   TGraph *gr_rate = new TGraph(nxsec,thv,rv);
   
   TH1F *hi_elec_xsec = (TH1F*)hi_elec_th->Clone();
   hi_elec_xsec->Scale(1/lumi/(3.14156/2/100));
   hi_elec_xsec->GetXaxis()->SetTitle("#theta (deg)");
   hi_elec_xsec->GetYaxis()->SetTitle("d#sigma/d#theta (mbarn/rad)");

   TH1F *hi_elec_xsec5 = (TH1F*)hi_elec_th5->Clone();
   hi_elec_xsec5->Scale(1/lumi/(3.14156*5/180/100));
   hi_elec_xsec5->GetXaxis()->SetTitle("#theta (deg)");
   hi_elec_xsec5->GetYaxis()->SetTitle("d#sigma/d#theta (mbarn/rad)");

   TH1F *hi_seco_xsec = (TH1F*)hi_seco_th->Clone();
   hi_seco_xsec->Scale(1/lumi/(3.14156/2/100));
   hi_seco_xsec->GetXaxis()->SetTitle("#theta (deg)");
   hi_seco_xsec->GetYaxis()->SetTitle("d#sigma/d#theta (mbarn/rad)");

   TH1F *hi_seco_xsec5 = (TH1F*)hi_seco_th5->Clone();
   hi_seco_xsec5->Scale(1/lumi/(3.14156*5/180/100));
   hi_seco_xsec5->GetXaxis()->SetTitle("#theta (deg)");
   hi_seco_xsec5->GetYaxis()->SetTitle("d#sigma/d#theta (mbarn/rad)");

   TCanvas *c1=new TCanvas("c1","Plot1",750,1000);
   c1->Divide(2,3);
   c1->cd(1);
   gPad->SetLogy();
   hi_elec_e->Draw("");
   c1->cd(2);
   gPad->SetLogy();
   hi_first_e->Draw("");
   c1->cd(3);
   gPad->SetLogy();
   hi_elec_th->Draw("");
   c1->cd(4);
   gPad->SetLogy();
   hi_first_th->Draw("");
   c1->cd(5);
   hi_elec_xy->Draw("COLZ");
   c1->cd(6);
   hi_first_xy->Draw("COLZ");
   c1->Print("moller.pdf(");


   TCanvas *c2=new TCanvas("c2","Plot2",750,1000);
   c2->Divide(2,3);
   c2->cd(1);
   gPad->SetLogy();
   hi_prim_e->Draw("");
   c2->cd(2);
   gPad->SetLogy();
   hi_seco_e->Draw("");
   c2->cd(3);
   gPad->SetLogy();
   hi_prim_th->Draw("");
   c2->cd(4);
   gPad->SetLogy();
   hi_seco_th->Draw("");
   c2->cd(5);
   hi_prim_xy->Draw("COLZ");
   c2->cd(6);
   hi_seco_xy->Draw("COLZ");
   c2->Print("moller.pdf");


   TCanvas *c3=new TCanvas("c3","Xsec",750,1000);
   c3->Divide(2,3);
   c3->cd(1);
   hi_elec_eth->Draw("COLZ");
   gr_eth->SetLineColor(2);
   gr_eth->Draw("L");
   c3->cd(2);
   hi_seco_eth->Draw("COLZ");
   gr_eth->SetLineColor(2);
   gr_eth->Draw("L");
   c3->cd(3);
   gPad->SetLogy();
   hi_elec_xsec->Draw("");
   c3->cd(4);
   gPad->SetLogy();
   hi_seco_xsec->Draw("");
   gr_xsec->SetLineColor(2);
   gr_xsec->Draw("L");
   c3->cd(5);
   gPad->SetLogy();
   hi_elec_xsec5->Draw("");
   gr_xsec->SetLineColor(2);
   gr_xsec->Draw("L");
   c3->cd(6);
   gPad->SetLogy();
   hi_seco_xsec5->Draw("");
   gr_xsec->SetLineColor(2);
   gr_xsec->Draw("L");
   c3->Print("moller.pdf");
   
   TCanvas *c4=new TCanvas("c4","Rate",750,1000);
   c4->Divide(1,2);
   c4->cd(1);
   hi_rate_th5->Draw("");
   gr_rate->SetLineColor(2);
   gr_rate->Draw("L");
   c4->cd(2);
   hi_rate_ft->Draw("");
   c4->Print("moller.pdf");
   
   TCanvas *c5=new TCanvas("c5","Spot",1000,1500);
   c5->Divide(NX,NY);
   for(int i=0; i<NPLANES; i++){
        int spot_m_maxbin = hi_spot_m_gt200[i]->GetMaximumBin();
	spot_m_maxpos[i] = hi_spot_m_gt200[i]->GetBinCenter(spot_m_maxbin);
	spot_m_maxval[i] = hi_spot_m_gt200[i]->GetBinCenter(spot_m_maxbin);
	TF1 *myg = new TF1("myg","gaus",spot_m_maxpos[i]-0.5,spot_m_maxpos[i]+2);
	myg->SetLineWidth(1);
	myg->SetLineColor(1);
	hi_spot_m_gt200[i]->Fit("myg","R0");
	spot_m_size[i]=myg->GetParameter(1)+2*myg->GetParameter(2);
	spot_m_angle[i]=57.3*atan(spot_m_size[i]/zplane[i]);
	if(i>=NMIN && i<NMAX) {
	   c5->cd(i-NMIN+1);
	   gPad->SetLogy();
	   hi_spot[i]->SetMinimum(1000.);
	   hi_spot[i]->Draw("");
	   hi_spot_e[i]->SetLineColor(2);
	   hi_spot_e[i]->Draw("SAME");
	   hi_spot_m[i]->SetLineColor(kOrange+7);
	   hi_spot_m[i]->Draw("SAME");
	   hi_spot_m_gt200[i]->SetLineColor(kOrange-2);
	   hi_spot_m_gt200_clone[i]=(TH1F*)hi_spot_m_gt200[i]->Clone();
	   myg->Draw("SAME");
	   hi_spot_g[i]->SetLineColor(4);
	   hi_spot_g[i]->Draw("SAME");
	   hi_spot_o[i]->SetLineColor(3);
	   hi_spot_o[i]->Draw("SAME");
	   TLine *lplane = new TLine(rclear[i],0.001,rclear[i],100000000.); lplane->Draw();
	}
   }
   c5->Print("moller.pdf");
   
   TCanvas *c6=new TCanvas("c6","Background Edep",1000,1500);
   c6->Divide(NX,NY);
   for(int i=NMIN; i<NMAX; i++){
	c6->cd(i-NMIN+1);
	gPad->SetLogy();
	hi_spot_ene[i]->Draw("");
	hi_spot_e_ene[i]->SetLineColor(2);
	hi_spot_e_ene[i]->Draw("SAME");
	hi_spot_m_ene[i]->SetLineColor(kOrange+7);
	hi_spot_m_ene[i]->Draw("SAME");
	hi_spot_m_gt200_ene[i]->SetLineColor(kOrange-2);
	hi_spot_m_gt200_ene[i]->Draw("SAME");
	hi_spot_g_ene[i]->SetLineColor(4);
	hi_spot_g_ene[i]->Draw("SAME");
	hi_spot_o_ene[i]->SetLineColor(3);
	hi_spot_o_ene[i]->Draw("SAME");
	TLine *lplane = new TLine(rclear[i],0.001,rclear[i],100000000.); lplane->Draw();
   }
   c6->Print("moller.pdf");
   
   TCanvas *c61=new TCanvas("c61","Background 2d",1000,1500);
   gPad->SetLogz();
   c61->Divide(NX,NY);
   for(int i=NMIN; i<NMAX; i++){
	c61->cd(i-NMIN+1);
	gPad->SetLogz();
	hi_spot_2d[i]->Draw("COLZ");
	TLine *lplane = new TLine(rclear[i],0.001,rclear[i],100000000.); lplane->Draw();
   }
   c61->Print("moller.pdf");
   
   TCanvas *c62=new TCanvas("c62","Electron 2d",1000,1500);
   gPad->SetLogz();
   c62->Divide(NX,NY);
   for(int i=NMIN; i<NMAX; i++){
	c62->cd(i-NMIN+1);
	hi_spot_e_2d[i]->Draw("COLZ");
	TLine *lplane = new TLine(rclear[i],0.001,rclear[i],100000000.); lplane->Draw();
   }
   c62->Print("moller.pdf");
   
   TCanvas *c63=new TCanvas("c63","Moller 2d",1000,1500);
   gPad->SetLogz();
   c63->Divide(NX,NY);
   for(int i=NMIN; i<NMAX; i++){
	c63->cd(i-NMIN+1);
	hi_spot_m_2d[i]->Draw("COLZ");
	TLine *lplane = new TLine(rclear[i],0.001,rclear[i],100000000.); lplane->Draw();
   }
   c63->Print("moller.pdf");
   
   TCanvas *c64=new TCanvas("c64","Gamma 2d",1000,1500);
   gPad->SetLogz();
   c64->Divide(NX,NY);
   for(int i=NMIN; i<NMAX; i++){
	c64->cd(i-NMIN+1);
	hi_spot_g_2d[i]->Draw("COLZ");
	TLine *lplane = new TLine(rclear[i],0.001,rclear[i],100000000.); lplane->Draw();
  }
   c64->Print("moller.pdf");
   
   TCanvas *c65=new TCanvas("c65","Other 2d",1000,1500);
   c65->Divide(NX,NY);
   for(int i=NMIN; i<NMAX; i++){
	c65->cd(i-NMIN+1);
	hi_spot_o_2d[i]->Draw("COLZ");
	TLine *lplane = new TLine(rclear[i],0.001,rclear[i],100000000.); lplane->Draw();
   }
   c65->Print("moller.pdf");
   
   TCanvas *c66=new TCanvas("c66","Good PLane",750,1000);
   c66->Divide(2,3);
   c66->cd(1);
   gPad->SetLogy();
   hi_spot[NGOOD]->Draw("");
   hi_spot_e[NGOOD]->SetLineColor(2);
   hi_spot_e[NGOOD]->Draw("SAME");
   hi_spot_m[NGOOD]->SetLineColor(kOrange+7);
   hi_spot_m[NGOOD]->Draw("SAME");
   hi_spot_m_gt200_clone[NGOOD]->SetLineColor(kOrange-2);
   hi_spot_m_gt200_clone[NGOOD]->Draw("SAME");
   hi_spot_g[NGOOD]->SetLineColor(4);
   hi_spot_g[NGOOD]->Draw("SAME");
   hi_spot_o[NGOOD]->SetLineColor(3);
   hi_spot_o[NGOOD]->Draw("SAME");
   TLine *lplane = new TLine(rclear[NGOOD],0.001,rclear[NGOOD],100000000.); lplane->Draw();
   c66->cd(2);
   gPad->SetLogy();
   hi_spot_ene[NGOOD]->Draw("");
   hi_spot_e_ene[NGOOD]->SetLineColor(2);
   hi_spot_e_ene[NGOOD]->Draw("SAME");
   hi_spot_m_ene[NGOOD]->SetLineColor(kOrange+7);
   hi_spot_m_ene[NGOOD]->Draw("SAME");
   hi_spot_m_gt200_ene[NGOOD]->SetLineColor(kOrange-2);
   hi_spot_m_gt200_ene[NGOOD]->Draw("SAME");
   hi_spot_g_ene[NGOOD]->SetLineColor(4);
   hi_spot_g_ene[NGOOD]->Draw("SAME");
   hi_spot_o_ene[NGOOD]->SetLineColor(3);
   hi_spot_o_ene[NGOOD]->Draw("SAME");
   lplane->Draw();
   c66->cd(3);
   gPad->SetLogz();
   hi_spot_2d[NGOOD]->Draw("COLZ");
   lplane->Draw();
   c66->cd(4);
   gPad->SetLogz();
   hi_spot_e_2d[NGOOD]->Draw("COLZ");
   lplane->Draw();
   c66->cd(5);
   gPad->SetLogz();
   hi_spot_g_2d[NGOOD]->Draw("COLZ");
   lplane->Draw();
   c66->cd(6);
   gPad->SetLogz();
   hi_spot_o_2d[NGOOD]->Draw("COLZ");
   lplane->Draw();
   c66->Print("moller.pdf");

   TCanvas *c67=new TCanvas("c67","Good PLane 2",1000,500);
   c67->Divide(2,1);
   c67->cd(1);
   gPad->SetLogy();
   hi_spot[NGOOD]->Draw("");
   hi_spot_e[NGOOD]->SetLineColor(2);
   hi_spot_e[NGOOD]->Draw("SAME");
   hi_spot_m[NGOOD]->SetLineColor(kOrange+7);
   hi_spot_m[NGOOD]->Draw("SAME");
   hi_spot_m_gt200_clone[NGOOD]->SetLineColor(kOrange-2);
   hi_spot_m_gt200_clone[NGOOD]->Draw("SAME");
   hi_spot_g[NGOOD]->SetLineColor(4);
   hi_spot_g[NGOOD]->Draw("SAME");
   hi_spot_o[NGOOD]->SetLineColor(3);
   hi_spot_o[NGOOD]->Draw("SAME");
   c67->cd(2);
   gPad->SetLogz();
   hi_spot_2d[NGOOD]->Draw("COLZ");
   c67->Print("moller.pdf");
     
  TCanvas *c7=new TCanvas("c7","SpotZ",1000,1000);
   c7->Divide(NX,NY);
   for(int i=NMIN; i<NMAX; i++){
	c7->cd(i-NMIN+1);
	gPad->SetLogy();
	hi_spotz[i]->Draw("");
   }
   c7->Print("moller.pdf");

   TCanvas *c8=new TCanvas("c8","Neutrons",500,1000);
   c8->Divide(1,3);
   c8->cd(1);
   hi_n_e->Draw("");
   c8->cd(2);
   hi_n_th->Draw("");
   c8->cd(3);
   hi_n_eth->Draw("COLZ");
   c8->Print("moller.pdf");

   TCanvas *c9=new TCanvas("c9","Rates",700,1000);
   c9->Divide(1,2);
   c9->cd(1);
   TGraph *gr_25_50 = new TGraph(32,zplane,rate_25_50);
   gr_25_50->SetLineColor(4);
   gr_25_50->SetLineWidth(1);
   gr_25_50->SetMarkerColor(4);
   gr_25_50->SetMarkerStyle(21);
   gr_25_50->GetXaxis()->SetTitle("Distance (mm)");
   gr_25_50->GetYaxis()->SetTitle("Rate (Hz)");
   gr_25_50->SetMinimum(0.);
   TGraph *gr_50_350 = new TGraph(32,zplane,rate_50_350);
   gr_50_350->SetLineColor(2);
   gr_50_350->SetLineWidth(1);
   gr_50_350->SetMarkerColor(2);
   gr_50_350->SetMarkerStyle(20);
   gr_50_350->GetXaxis()->SetTitle("Distance (mm)");
   gr_50_350->GetYaxis()->SetTitle("Rate (Hz)");
   gr_50_350->SetMinimum(0.);
   if(gr_50_350->GetHistogram()->GetMaximum()>gr_25_50->GetHistogram()->GetMaximum()) gr_25_50->SetMaximum(gr_50_350->GetHistogram()->GetMaximum()*1.1);
   gr_25_50->Draw("ACP");
   gr_50_350->Draw("SP");
   c9->cd(2);
   TGraph *gr_25_50_emin = new TGraph(32,zplane,rate_25_50_emin);
   gr_25_50_emin->SetLineColor(4);
   gr_25_50_emin->SetLineWidth(1);
   gr_25_50_emin->SetMarkerColor(4);
   gr_25_50_emin->SetMarkerStyle(21);
   gr_25_50_emin->GetXaxis()->SetTitle("Distance (mm)");
   gr_25_50_emin->GetYaxis()->SetTitle("Rate (GHz)");
   gr_25_50_emin->SetMinimum(0.);
   gr_25_50_emin->Draw("ACP");
   TGraph *gr_50_350_emin = new TGraph(32,zplane,rate_50_350_emin);
   gr_50_350_emin->SetLineColor(2);
   gr_50_350_emin->SetLineWidth(1);
   gr_50_350_emin->SetMarkerColor(2);
   gr_50_350_emin->SetMarkerStyle(20);
   gr_50_350_emin->GetXaxis()->SetTitle("Distance (mm)");
   gr_50_350_emin->GetYaxis()->SetTitle("Rate (GHz)");
   gr_50_350_emin->SetMinimum(0.);
   if(gr_50_350_emin->GetHistogram()->GetMaximum()>gr_25_50_emin->GetHistogram()->GetMaximum()) gr_25_50_emin->SetMaximum(gr_50_350_emin->GetHistogram()->GetMaximum()*1.1);
   gr_25_50_emin->Draw("ACP");
   gr_50_350_emin->Draw("SP");
   c9->Print("moller.pdf");

   TCanvas *c10=new TCanvas("c10","Rates",700,1000);
   c10->Divide(1,2);
   c10->cd(1);
   TGraph *gr_spot_size = new TGraph(NPLANES,zplane,spot_m_size);
   gr_spot_size->SetTitle("Spot size");
   gr_spot_size->SetLineColor(4);
   gr_spot_size->SetLineWidth(1);
   gr_spot_size->SetMarkerColor(4);
   gr_spot_size->SetMarkerStyle(21);
   gr_spot_size->GetXaxis()->SetTitle("Distance (cm)");
   gr_spot_size->GetYaxis()->SetTitle("R (cm)");
   gr_spot_size->SetMinimum(0.);
   gr_spot_size->Draw("ACP");
   TF1 *clearance = new TF1("clearance","(x+3)*tan(5/57.3)-1.129",zplane[NGOOD],zplane[NPLANES-1]); clearance->Draw("SAME");
   TLine *lcone = new TLine(71.64,3.81,71.64,5.40); lcone->Draw();
   TLine *ltip  = new TLine(57.00,3.90,57.00,4.12); ltip->Draw();
   TLine *llin1 = new TLine(57.00,3.90,71.64,3.81); llin1->Draw();
   TLine *llin2 = new TLine(57.00,4.12,71.64,5.40); llin2->Draw();
   TLine *llin3 = new TLine(71.64,3.81,95.00,3.81); llin3->Draw();
   TLine *llin4 = new TLine(78.74,3.38,95.00,3.38); llin4->Draw();
   gr_spot_size->GetYaxis()->SetRangeUser(0,5.4);
   c10->cd(2);
   TGraph *gr_spot_angle = new TGraph(NPLANES,zplane,spot_m_angle);
   gr_spot_angle->SetTitle("Spot angle view");
   gr_spot_angle->SetLineColor(4);
   gr_spot_angle->SetLineWidth(1);
   gr_spot_angle->SetMarkerColor(4);
   gr_spot_angle->SetMarkerStyle(21);
   gr_spot_angle->GetXaxis()->SetTitle("Distance (cm)");
   gr_spot_angle->GetYaxis()->SetTitle("Angle (deg)");
   gr_spot_angle->SetMinimum(0.);
   gr_spot_angle->Draw("ACP");
   c10->Print("moller.pdf)");

// saving histograms to file
   TIter next1(gDirectory->GetList());
   TH1 *my_hi[1000];
   int ihi=0;
   TObject* obj;
   while((obj= (TObject*)next1())){
     if(obj->InheritsFrom(TH1::Class())){
       my_hi[ihi] = (TH1*)obj;
       ihi++;
     }
   }
   TFile myfile("moller.root", "recreate");
   for(int i=0; i<ihi; i++) my_hi[i]->Write();
   myfile.Close();

   FILE *fp = fopen("moller.txt","w");
   fprintf(fp,"%d\n",NPLANES);
   for(int i=0; i<NPLANES;i++) {
     fprintf(fp,"%5.3f  %5.3f  %5.3f  %5.3f  %5.3f \n",zplane[i],rate_25_50[i]*1000,rate_50_350[i]*1000,spot_m_size[i],spot_m_angle[i]);
   }
   fclose(fp);
   gui.Run(1);
}   

