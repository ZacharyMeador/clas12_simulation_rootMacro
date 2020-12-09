#include <cstdlib>
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
#include <TMath.h>
#include "TFile.h"
#include "TApplication.h"
#include "TLatex.h"
#include <TROOT.h>

void plot_ic_55();

using namespace std; 
TApplication gui("GUI",0,NULL);

int main()
{
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasColor(10);

   gStyle->SetPadBorderMode(0);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadRightMargin(0.16);
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
  




   // header  & generated
   vector<int>   *evn = new vector<int>;
   vector<double> *px = new vector<double>;
   vector<double> *py = new vector<double>;
   vector<double> *pz = new vector<double>;
   
   
   //ft-cal
   vector<double>  *cal_pid = new vector<double>;
   vector<double> *cal_Edep = new vector<double>;
   vector<double>    *cal_E = new vector<double>;
   vector<double>    *cal_t = new vector<double>;
   vector<double>    *cal_x = new vector<double>;
   vector<double>    *cal_y = new vector<double>;
   vector<double>    *cal_z = new vector<double>;
   vector<double>  *cal_adc = new vector<double>;
   vector<double>  *cal_tdc = new vector<double>;
   vector<double>   *cal_id = new vector<double>;
	
   //ft-hodo
   vector<double>   *hodo_id = new vector<double>;
   vector<double>  *hodo_pid = new vector<double>;
   vector<double> *hodo_Edep = new vector<double>;
   vector<double>    *hodo_E = new vector<double>;
   vector<double>    *hodo_t = new vector<double>;
   vector<double>    *hodo_x = new vector<double>;
   vector<double>    *hodo_y = new vector<double>;
   vector<double>    *hodo_z = new vector<double>;
	

   
   float Edep_tot,Edep_tot_int,Edep_tot_time,Emax;
   int   nhit_int, nhit_time, nthr_time;
   float e_pd[32*32];
   float t_pd[32*32];
   float adc[32*32];
   float tdc[32*32];
   int   hitS[32*32][64];
   float edepS[32*32][64];
   float timeS[32*32][64];
   float theta_hit,theta_ic;
   float threshold;
   int   trigger,fill_single,n;

   fill_single=1;

   TChain *h = new TChain("header");
   h->Add("out*.root");
   TChain *g = new TChain("generated");
   g->Add("ou*.root");
   TChain *ftcal = new TChain("ft_cal");
   ftcal->Add("out*.root");
   TChain *fthodo = new TChain("ft_hodo");
   fthodo->Add("out*.root"); 
  
   // header bank
   h->SetBranchAddress("evn",&evn);

   // generated
   g->SetBranchAddress("px",&px);
   g->SetBranchAddress("py",&py);
   g->SetBranchAddress("pz",&pz);

   // ft-cal
   ftcal->SetBranchAddress("trackE",&cal_E);
   ftcal->SetBranchAddress("totEdep",&cal_Edep);
   ftcal->SetBranchAddress("avg_t",&cal_t);
   ftcal->SetBranchAddress("pid",&cal_pid);
   ftcal->SetBranchAddress("avg_x",&cal_x);
   ftcal->SetBranchAddress("avg_y",&cal_y);
   ftcal->SetBranchAddress("avg_z",&cal_z);
   ftcal->SetBranchAddress("adc",&cal_adc);
   ftcal->SetBranchAddress("tdc",&cal_tdc);
   ftcal->SetBranchAddress("component",&cal_id);

   // ft-hodo
   fthodo->SetBranchAddress("component",&hodo_id);
   fthodo->SetBranchAddress("trackE",&hodo_E);
   fthodo->SetBranchAddress("totEdep",&hodo_Edep);
   fthodo->SetBranchAddress("avg_t",&hodo_t);
   fthodo->SetBranchAddress("pid",&hodo_pid);
   fthodo->SetBranchAddress("avg_x",&hodo_x);
   fthodo->SetBranchAddress("avg_y",&hodo_y);
   fthodo->SetBranchAddress("avg_z",&hodo_z);
 


   
// Create histos

   TH1F *hi_etot  = new TH1F("Etot", "", 100, 0., 2.5);
   hi_etot->GetXaxis()->SetTitle("E_{tot} (GeV)");
   hi_etot->GetYaxis()->SetTitle("Counts");
   TH1F *hi_iedep = new TH1F("FT Edep","",200,0.,1000.);
   hi_iedep->GetXaxis()->SetTitle("E_{hit} (MeV)");
   hi_iedep->GetYaxis()->SetTitle("Counts");
   TH1F *hi_time = new TH1F("FT Time","",100,0.,250.);
   hi_time->GetXaxis()->SetTitle("Time (ns)");
   hi_time->GetYaxis()->SetTitle("Counts");
   TH2F *hi_e_zic   = new TH2F("Edep vs Z_ic", "", 100,  180., 210., 100, 0., 50.0);
   hi_e_zic->GetXaxis()->SetTitle("z_{hit} (cm)");
   hi_e_zic->GetYaxis()->SetTitle("E_{hit} (MeV)");
   TH1F *hi_zic     = new TH1F("Z_ic", "", 100,  180., 210.);
   hi_zic->GetXaxis()->SetTitle("z_{hit} (cm)");
   hi_zic->GetYaxis()->SetTitle("Counts");
   TH1F *hi_ith   = new TH1F("FT Theta","",100,0.,10.);
   hi_ith->GetXaxis()->SetTitle("#Theta_{hit} (deg)");
   hi_ith->GetYaxis()->SetTitle("Counts");
   TH2F *hi_ixy   = new TH2F("FT Occupancy","",25,0.,25.,25,0.,25.);
   hi_ixy->GetXaxis()->SetTitle("X");
   hi_ixy->GetYaxis()->SetTitle("Y");
   TH2F *hi_iratexy   = new TH2F("FT Rate","",25,0.,25.,25,0.,25.);
   hi_iratexy->GetXaxis()->SetTitle("X");
   hi_iratexy->GetYaxis()->SetTitle("Y");
   TH2F *hi_rxy   = new TH2F("FT Dose","",25,0.,25.,25,0.,25.);
   hi_rxy->GetXaxis()->SetTitle("X");
   hi_rxy->GetYaxis()->SetTitle("Y");
   TH1F *hi_icrystal = new TH1F("FT Occupancy 1D","",600,0.,600.);
   hi_icrystal->GetXaxis()->SetTitle("X");
   hi_icrystal->GetYaxis()->SetTitle("E(MeV/ns)");
   TH1F *hi_itime = new TH1F("FT Edep vs, Time","",100,0.,250.);
   hi_itime->GetXaxis()->SetTitle("Time(ns)");
   hi_itime->GetYaxis()->SetTitle("E(MeV/ns)");

   TH1F *hi_etot_int  = new TH1F("Etot Integrated", "", 100, 0., 2.5);
   hi_etot_int->GetXaxis()->SetTitle("E_{tot} (GeV)");
   hi_etot_int->GetYaxis()->SetTitle("Counts");
   TH1F *hi_nhit_int= new TH1F("Nhit Integrated","",300,0.,3000.);
   hi_nhit_int->GetXaxis()->SetTitle("Nr. Hits");
   hi_nhit_int->GetYaxis()->SetTitle("Counts");
   TH1F *hi_iedep_int = new TH1F("FT Edep Integrated","",200,0.,1000.);
   hi_iedep_int->GetXaxis()->SetTitle("E_{hit} (MeV)");
   hi_iedep_int->GetYaxis()->SetTitle("Counts");
   TH2F *hi_nxy_int   = new TH2F("FT Occupancy Integrated","",25,0.,25.,25,0.,25.);
   hi_nxy_int->GetXaxis()->SetTitle("X");
   hi_nxy_int->GetYaxis()->SetTitle("Y");
   TH2F *hi_exy_int   = new TH2F("FT Energy Occupancy Integrated","",25,0.,25.,25,0.,25.);
   hi_exy_int->GetXaxis()->SetTitle("X");
   hi_exy_int->GetYaxis()->SetTitle("Y");
   TH1F *hi_etime_int = new TH1F("Edep vs, Time Integrated","",64,0.,256.);
   hi_etime_int->GetXaxis()->SetTitle("Time(ns)");
   hi_etime_int->GetYaxis()->SetTitle("E(MeV)");
   TH1F *hi_ntime_int = new TH1F("FT Hits vs, Time Integrated","",64,0.,256.);
   hi_ntime_int->GetXaxis()->SetTitle("Time(ns)");
   hi_ntime_int->GetYaxis()->SetTitle("N. hit");

   TH1F *hi_etot_time  = new TH1F("Etot in Time", "", 100, 0., 0.1);
   hi_etot_time ->GetXaxis()->SetTitle("E_{tot} (GeV)");
   hi_etot_time ->GetYaxis()->SetTitle("Counts");
   TH1F *hi_nhit_time= new TH1F("Nhits in Time","",50,0.,50.);
   hi_nhit_time ->GetXaxis()->SetTitle("Nr. Hits");
   hi_nhit_time ->GetYaxis()->SetTitle("Counts");
   TH1F *hi_nthr_time= new TH1F("Nhit above threshold","",10,0.,10.);
   hi_nthr_time ->GetXaxis()->SetTitle("Nr. Hits");
   hi_nthr_time ->GetYaxis()->SetTitle("Counts");
   TH1F *hi_iedep_time = new TH1F("FT Edep Time","",100,0.,100.);
   hi_iedep_time->GetXaxis()->SetTitle("E_{hit} (MeV)");
   hi_iedep_time->GetYaxis()->SetTitle("Counts");
   TH2F *hi_nxy_time     = new TH2F("FT Occupancy Time","",25,0.,25.,25,0.,25.);
   hi_nxy_time ->GetXaxis()->SetTitle("X");
   hi_nxy_time ->GetYaxis()->SetTitle("Y");
   TH2F *hi_exy_time  = new TH2F("FT Energy Occupancy Time","",25,0.,25.,25,0.,25.);
   hi_exy_time ->GetXaxis()->SetTitle("X");
   hi_exy_time ->GetYaxis()->SetTitle("Y");

   TH1F *hi_adc = new TH1F("FT ADC","",200,0.,5000.);
   hi_adc->GetXaxis()->SetTitle("ADC (Ch.)");
   hi_adc->GetYaxis()->SetTitle("Counts");
   TH1F *hi_tdc = new TH1F("FT Channel","",200,0.,8000.);
   hi_tdc->GetXaxis()->SetTitle("TDC (Ch.)");
   hi_tdc->GetYaxis()->SetTitle("Counts");
   TH2F *hi_adc_edep = new TH2F("FT ADC vs. Edep","",100,0.,1000.,100,0.,5000.);
   hi_adc_edep->GetXaxis()->SetTitle("E_{hit} (MeV)");
   hi_adc_edep->GetYaxis()->SetTitle("ADC (Ch.)");
   TH2F *hi_tdc_time = new TH2F("FT TDC vs. Time","",100,0.,250.,100,0.,8000.);
   hi_tdc_time->GetXaxis()->SetTitle("Time (ns)");
   hi_tdc_time->GetYaxis()->SetTitle("TDC (Ch.)");

   //hodoscope
   TH2F *hi_h_ixy   = new TH2F("hi_h_ixy","",40,-20., 20.,40,-20.,20.);
   hi_h_ixy->GetXaxis()->SetTitle("X");
   hi_h_ixy->GetYaxis()->SetTitle("Y");
   TH2F *hi_h_rxy   = new TH2F("hi_h_rxy","",40,-20., 20.,40,-20.,20.);
   hi_h_rxy->GetXaxis()->SetTitle("X");
   hi_h_rxy->GetYaxis()->SetTitle("Y");


   TH2F *hi_se[12];
   float p_se[12];
   float th_se[12];
   float p_se0 = 0;
   float th_se0 = 0;
   char ptext[40];
   char thtext[40];
   TLatex l_p_se;
   TLatex l_th_se;
   l_p_se.SetTextSize(0.05);
   l_th_se.SetTextSize(0.05);
   l_p_se.SetTextFont(72);
   l_th_se.SetTextFont(72);	    
   for(int i=0; i<12; i++){
     sprintf(ptext,"FTC Occupancy %d",i);
     hi_se[i] = new TH2F(ptext,"",25,0.,25.,25,0.,25.);
     hi_se[i]->GetXaxis()->SetTitle("X");
     hi_se[i]->GetYaxis()->SetTitle("Y");
   }
   TH2F *hi_se0_nocut = new TH2F("FTC Occupancy","",25,0.,25.,25,0.,25.);
   hi_se0_nocut->GetXaxis()->SetTitle("X");
   hi_se0_nocut->GetYaxis()->SetTitle("Y");
   TH2F *hi_se0_time = new TH2F("FTC Occupancy time","",25,0.,25.,25,0.,25.);
   hi_se0_time->GetXaxis()->SetTitle("X");
   hi_se0_time->GetYaxis()->SetTitle("Y");
   
   int nint=0;
   int nfull=0;
   int ngoodentries=0;
  
   Long64_t nentries = g->GetEntries();	
   cout << nentries << endl;
   threshold=10;
   int   nelec=124000;
   float twin=250;
   cout << " Number of events: " << nentries << endl;
   cout << " Number of electrons per event: " << nelec << endl;
   cout << " Time window: " << twin << " ns" << endl;
   float normL=(124000./250.)/(nelec/twin);
   int   nsum=(int) normL;
   cout << " Number of events to integrate = " << nsum  << endl;
   cout << " Luminosity Scale Factor: " << normL << endl;

   for(Long64_t jentry=0; jentry < nentries; jentry++) {
  
       if(int(jentry/1000)*1000==jentry) cout << "Analyzed " << jentry << " events of " << nentries << endl;

        // clear vectors
        evn->clear();
	px->clear();
	py->clear();
	pz->clear();
	cal_pid->clear();
	cal_E->clear();
	cal_Edep->clear();
	cal_t->clear();
	cal_x->clear();
	cal_y->clear();
	cal_z->clear();
	cal_adc->clear();
	cal_tdc->clear();
	cal_id->clear();
	hodo_id->clear();
	hodo_pid->clear();
	hodo_E->clear();
	hodo_Edep->clear();
	hodo_t->clear();
	hodo_x->clear();
	hodo_y->clear();
	hodo_z->clear();

	g->GetEntry(jentry);
	ftcal->GetEntry(jentry); 
       	fthodo->GetEntry(jentry); 

     trigger=0;
     n=0;

//     cout << " " << endl;
//     cout << "new event" << endl;
     //     cout << p[0] << " " <<  theta[0] << " " << vz[0] << endl;
     // IC
     // analyzing hits
     theta_ic=0;
     Emax=0;
     Edep_tot=0;
     n=0;
     for(int i=0; i<32*32; i++) {
       e_pd[i]=0;
       t_pd[i]=600;
       adc[i]=0;
       tdc[i]=0;
     }
     
     /*     p[0]=p[0]/1000;
     if( p[0]-2. <0.1 && p[0]-2.>-0.1 && theta[0]> 2.7 && theta[0]<2.9 && fill_single==0 ) {
       cout << p[0] << " " <<  theta[0] << " " << (p[0]-2.) << endl;
       fill_single++;
       p_se0=p[0];
       th_se0=theta[0];
       }*/
     
     //     cout << " new event" << endl;
     int inhit=cal_adc->size();
     for(int i=0; i<inhit; i++){
       // complete hit digitization summing hits in the same crystal
       int idy = (*cal_id)[i]/22;
       int idx = (*cal_id)[i]-idy*22+1;
       idy=idy+1;
       
       // looking for maximum energy deposition
       if((*cal_Edep)[i]>Emax)  Emax=(*cal_Edep)[i];
       
       // get integrated quantities
       // cout << Edep_tot << " " <<  (*cal_Edep)[i] << endl;
       Edep_tot=Edep_tot+(*cal_Edep)[i];
       // hit related info
       theta_hit=atan(sqrt((*cal_x)[i]*(*cal_x)[i]+(*cal_y)[i]*(*cal_y)[i])/(*cal_z)[i])*57.3;
       hi_iedep->Fill((*cal_Edep)[i],1/twin);
       hi_time->Fill((*cal_t)[i],1/twin);
       hi_e_zic->Fill((*cal_z)[i]/10,(*cal_Edep)[i],(*cal_Edep)[i]/twin);
       hi_zic->Fill((*cal_z)[i]/10,(*cal_Edep)[i]/twin);
       hi_ith->Fill(theta_hit,(*cal_Edep)[i]/twin);
       hi_ixy->Fill(idx*1.,idy*1.,(*cal_Edep)[i]/twin);
       if((*cal_Edep)[i]>threshold) hi_iratexy->Fill(idx*1.,idy*1.,1/twin);
       hi_rxy->Fill(idx*1.,idy*1.,(*cal_Edep)[i]*154./twin);
       hi_icrystal->Fill(idx*1.+(idy-1)*23.,(*cal_Edep)[i]/twin);
       hi_itime->Fill((*cal_t)[i],(*cal_Edep)[i]/twin);
       hi_adc->Fill((*cal_adc)[i],1/twin);
       hi_tdc->Fill((*cal_tdc)[i],1/twin);
       hi_adc_edep->Fill((*cal_Edep)[i],(*cal_adc)[i],1/twin);
       hi_tdc_time->Fill((*cal_t)[i],(*cal_tdc)[i],1/twin);
     }
     hi_etot->Fill(Edep_tot/1000.);

     int hnhit=hodo_pid->size();
     //     cout << hnhit << endl;
     for(int i=0; i<hnhit; i++){
       if((*hodo_id)[i]<5000) {
	 hi_h_ixy->Fill((*hodo_x)[i]/10.,(*hodo_y)[i]/10.,(*hodo_Edep)[i]/twin); 
	 hi_h_rxy->Fill((*hodo_x)[i]/10.,(*hodo_y)[i]/10.,(*hodo_Edep)[i]*57600./2/twin);
       }
     }

     // summing up events to get the right luminosity
     if(nint<nsum){
       for(int i=0; i<inhit; i++){
	 int idy = (*cal_id)[i]/22;
	 int idx = (*cal_id)[i]-idy*22+1;
	 idy=idy+1;
	 int icrystal=32*idy+idx;
         int itime=(*cal_t)[i]/4.;
         if(itime<64 && (*cal_Edep)[i]>0) {
	   hitS[icrystal][itime] =hitS[icrystal][itime]+1;
	   edepS[icrystal][itime]=edepS[icrystal][itime]+(*cal_Edep)[i];
	   timeS[icrystal][itime]=timeS[icrystal][itime]+(*cal_t)[i];
	 }
       }
       nint++;
     }
     if(nint==nsum) {
       Edep_tot_int=0;
       nhit_int=0;
       Edep_tot_time=0;
       nhit_time=0;
       nthr_time=0;
       for(int i=0; i<32*32; i++) {
	 for(int j=0; j<64; j++) {
	   if(hitS[i][j]>0) {
             nhit_int++;
	     Edep_tot_int=Edep_tot_int+edepS[i][j];
	     int iiy=int(i/32.);
	     int iix=i-32*iiy;
             timeS[i][j]=timeS[i][j]/hitS[i][j];
             hi_iedep_int->Fill(edepS[i][j],1./twin);
             hi_etime_int->Fill(timeS[i][j],edepS[i][j]);
             hi_ntime_int->Fill(timeS[i][j]);
             hi_nxy_int->Fill(iix*1.,iiy*1.,1./twin);
             hi_exy_int->Fill(iix*1.,iiy*1.,edepS[i][j]/twin);
	     if(fill_single>0) hi_se0_nocut->Fill(iix*1.,iiy*1.,edepS[i][j]);
//	     if(j==32 || j==33 || j==34 || j==35 || j==36) {
	     if(j>15 && j <=44) {
	       nhit_time++;
	       if(edepS[i][j]>threshold) nthr_time++;
	       Edep_tot_time=Edep_tot_time+edepS[i][j];
	       hi_iedep_time->Fill(edepS[i][j]);
	       hi_nxy_time->Fill(iix*1.,iiy*1.);
	       hi_exy_time->Fill(iix*1.,iiy*1.,edepS[i][j]);
	       if(nfull<12) {
		 hi_se[nfull]->Fill(iix*1.,iiy*1.,edepS[i][j]);
		 //		 cout << nhit_time << " " << iix << " " << iiy << " " << edepS[i][j] << endl;
	       }
	       if(fill_single>0) hi_se0_time->Fill(iix*1.,iiy*1.,edepS[i][j]);
	     }
 	   }
	   hitS[i][j]=0;
	   edepS[i][j]=0;
	   timeS[i][j]=0;
	 }
       }
       fill_single=0;
       hi_etot_int->Fill(Edep_tot_int/1000.);
       hi_nhit_int->Fill(nhit_int*1.);
       hi_etot_time->Fill(Edep_tot_time/1000.);
       hi_nhit_time->Fill(nhit_time*1.);
       hi_nthr_time->Fill(nthr_time*1.);
       if(nfull<12) {
	 p_se[nfull]=Edep_tot_time/1000.;
	 th_se[nfull]=nhit_time;
       }
       nint=0;
       nfull++;
     }   
   }
   ngoodentries=nentries;
   float normS=normL/ngoodentries;
   float time=twin/normS;
   cout << " Run time: " << time  << " ns" << endl;
   cout << " Number of good events: " << ngoodentries << endl;
   cout << " Histogram Normalization Factor: " << normS << endl;
   cout << " Run time: " << time  << " ns" << endl;
   cout << " Number of full events: " << nfull << endl;
   cout << 1/normS << " " << nfull << endl;
   
   
   if(nfull>0) {
     hi_iedep->Scale(normS);
     hi_time->Scale(normS);
     hi_e_zic->Scale(normS);
     hi_zic->Scale(normS);
     hi_ith->Scale(normS);
     hi_ixy->Scale(normS);
     hi_iratexy->Scale(normS*1E6);
     hi_rxy->Scale(normS);
     hi_icrystal->Scale(normS);
     hi_adc->Scale(normS);
     hi_tdc->Scale(normS);
     hi_adc_edep->Scale(normS);
     hi_tdc_time->Scale(normS);

     hi_iedep_int->Scale(1./nfull);
     hi_etime_int->Scale(1./nfull);
     hi_ntime_int->Scale(1./nfull);
     hi_exy_int->Scale(1./nfull);
     hi_nxy_int->Scale(1./nfull);

     hi_iedep_time->Scale(1./nfull);
     hi_nxy_time->Scale(1./nfull);
     hi_exy_time->Scale(1./nfull);

     hi_h_ixy->Scale(normS);
     hi_h_rxy->Scale(normS);
   }
   

   
   TCanvas *c1=new TCanvas("c1","FT General",750,1000);
   c1->Divide(2,3);
   c1->cd(1);
   gPad->SetLogy();
   hi_etot->Draw("");
   c1->cd(2);
   gPad->SetLogy();
   hi_iedep->Draw("");  
   c1->cd(3);
   gPad->SetLogz();
   hi_e_zic->Draw("COLZ");
   c1->cd(4);
   hi_zic->Draw("");
   c1->cd(5);
   hi_ixy->Draw("COLZ");
   plot_ic_55();
   c1->cd(6);
   hi_ith->Draw("");
   c1->Print("ft_rate.ps(");
  
   TCanvas *c2=new TCanvas("c2","FT Integral",750,1000);
   c2->Divide(2,3);
   c2->cd(1);
   gPad->SetLogy();
   hi_etot_int->Draw("");
   c2->cd(2);
   hi_nhit_int->Draw("");
   c2->cd(3);
   gPad->SetLogy();
   hi_iedep_int->Draw("");
   c2->cd(4);
   hi_etime_int->Draw("");
   c2->cd(5);
   hi_nxy_int->Draw("COLZ");
   plot_ic_55();
   c2->cd(6);
   hi_exy_int->Draw("COLZ");
   plot_ic_55();
   c2->Print("ft_rate.ps");

   TCanvas *c3=new TCanvas("c3","FT Time",750,1000);
   c3->Divide(2,3);
   c3->cd(1);
   gPad->SetLogy();
   hi_etot_time->Draw("");
   c3->cd(2);
   gPad->SetLogy();
   hi_iedep_time->Draw("");
   c3->cd(3);
   hi_nhit_time->Draw("");
   c3->cd(4);
   hi_nthr_time->Draw("");
   c3->cd(5);
   gPad->SetLogz();
   hi_nxy_time->Draw("COLZ");
   plot_ic_55();
   c3->cd(6);
   gPad->SetLogz();
   hi_exy_time->Draw("COLZ");
   plot_ic_55();
   c3->Print("ft_rate.ps");

   TCanvas *c4=new TCanvas("c4","FT Digitization",750,750);
   c4->Divide(2,2);
   c4->cd(1);
   gPad->SetLogy();
   hi_adc->Draw("");
   c4->cd(2);
   gPad->SetLogy();
   hi_tdc->Draw("");
   c4->cd(3);
   hi_adc_edep->Draw("COLZ");
   c4->cd(4);
   hi_tdc_time->Draw("COLZ");
   c4->Print("ft_rate.ps");

  TCanvas *c5=new TCanvas("c5","Single Events",750,1000);
   c5->Divide(3,4);
   for(int i=0; i<12; i++){
     c5->cd(i+1);
     //     gPad->SetLogz();
     hi_se[i]->Draw("COLZ");
     sprintf(ptext,"E_{tot} = %1.2f GeV",p_se[i]);
     sprintf(thtext,"N_{crystals} = %d",(int) th_se[i]);
     l_p_se.DrawLatex(3.,28.,ptext);
     l_th_se.DrawLatex(3.,25.,thtext);
     plot_ic_55();
   }
   c5->Print("ft_rate.ps)");

   TCanvas *c7=new TCanvas("c7","gemc_moller",750,750);
   c7->Divide(2,2);
   c7->cd(1);
   gPad->SetLogy();
   hi_etot->Draw("");
   c7->cd(2);
   gPad->SetLogy();
   hi_iedep->Draw("");
   c7->cd(3);
   gPad->SetLogz();
   hi_ixy->Draw("COLZ");
   plot_ic_55();
   c7->cd(4);
   gPad->SetLogy();
   hi_icrystal->Draw("");
   c7->Print("gemc_moller.pdf");
   
   TCanvas *c8=new TCanvas("c8","coincidence",750,1000);
   c8->Divide(2,3);
   c8->cd(1);
   hi_etime_int->Draw("");
   c8->cd(2);
   hi_ntime_int->Draw("");
   c8->cd(3);
   gPad->SetLogz();
   hi_exy_int->Draw("COLZ");
   plot_ic_55();
   c8->cd(4);
   gPad->SetLogz();
   hi_exy_time->Draw("COLZ");
   plot_ic_55();
   c8->cd(5);
   gPad->SetLogz();
   hi_se0_nocut->Draw("COLZ");
   sprintf(ptext,"p = %1.2f GeV",p_se0);
   sprintf(thtext,"#Theta = %1.2f deg",th_se0);
   l_p_se.DrawLatex(3.,28.,ptext);
   l_th_se.DrawLatex(3.,25.,thtext);
   plot_ic_55();
   c8->cd(6);
   gPad->SetLogz();
   hi_se0_time->Draw("COLZ");
   sprintf(ptext,"p = %1.2f GeV",p_se0);
   sprintf(thtext,"#Theta = %1.2f deg",th_se0);
   l_p_se.DrawLatex(3.,28.,ptext);
   l_th_se.DrawLatex(3.,25.,thtext);
   plot_ic_55();
   c8->Print("gemc_coinc.pdf");

   gStyle->SetOptStat("");
   TCanvas *c0=new TCanvas("c0","Radiation Dose",500,500);
   hi_rxy->SetMaximum(5.);
   hi_rxy->Draw("COLZ");
   plot_ic_55();
   double radv[332]={0};
   double radx[332]={0};
   double rady[332]={0};
   double radxx[332]={0};
   double radyy[332]={0};
   int irad[332]={0};
   int radi=-1;
   for(int i=0; i<hi_rxy->GetNbinsX(); i++) {
       for(int j=0; j<hi_rxy->GetNbinsY(); j++) {
   	   int    rad_bin = hi_rxy->GetBin(i,j);
   	   double rad_val = hi_rxy->GetBinContent(rad_bin);
	   if(rad_val>0) {
	      radi++;
	      radv[radi] = rad_val;
	      radxx[radi]=i;
	      radyy[radi]=j;
	      if(i>12) radx[radi] = i-12;
	      else     radx[radi] = i-13;
	      if(j>12) rady[radi] = j-12;
	      else     rady[radi] = j-13;
	   }
	}
   }
   TMath::Sort(332,radv,irad);
   for(int i=0; i<332; i++) {
       cout << radxx[irad[i]] << "\t" << radyy[irad[i]] << "\t" << radx[irad[i]] << "\t" << rady[irad[i]] << "\t" << radv[irad[i]] << endl;
   }
   c0->Print("ft_rad.pdf");

   TCanvas *c9=new TCanvas("c9","Background Rate and Energy",1000,500);
   c9->Divide(2,1);
   c9->cd(1);
   gPad->SetLogz();
   hi_iratexy->Draw("COLZ");
   plot_ic_55();
   c9->cd(2);
   gPad->SetLogz();
   hi_ixy->Draw("COLZ");
   plot_ic_55();
   c9->Print("ft_rate.pdf");

   TCanvas *c10=new TCanvas("c10","Background Rate and Energy in Time",1000,500);
   c10->Divide(2,1);
   c10->cd(1);
   gPad->SetLogy();
   hi_nthr_time->Draw("");
   c10->cd(2);
   gPad->SetLogz();
   hi_exy_time->Draw("COLZ");
   plot_ic_55();
   c10->Print("ft_time.pdf");

   gStyle->SetOptStat("");
   TCanvas *c11=new TCanvas("c11","Radiation Dose - Hodoscope",500,750);
   c11->Divide(1,2);
   c11->cd(1);
   hi_h_ixy->Draw("COLZ");
   c11->cd(2);
   hi_h_rxy->Draw("COLZ");
   c11->Print("ft_rad_hodo.pdf");
  

   gui.Run(1);

}




void plot_ic_55(){
     TLine *icl= new TLine();
     icl= new TLine( 1., 9., 1.,15.); icl->Draw();
     icl= new TLine( 1.,15., 2.,15.); icl->Draw();
     icl= new TLine( 2.,15., 2.,18.); icl->Draw();
     for(int i=0; i<4; i++) { icl= new TLine(i*1.+2.,18.+i,i*1.+3.,18.+i); icl->Draw();}
     for(int i=0; i<4; i++) { icl= new TLine(i*1.+3.,18.+i,i*1.+3.,19.+i); icl->Draw();}
     icl= new TLine( 6.,22., 9.,22.); icl->Draw();
     icl= new TLine( 9.,22., 9.,23.); icl->Draw();
     icl= new TLine( 9.,23.,15.,23.); icl->Draw();
     icl= new TLine(15.,23.,15.,22.); icl->Draw();
     icl= new TLine(15.,22.,18.,22.); icl->Draw();
     for(int i=0; i<4; i++) { icl= new TLine(i*1.+18.,21.-i,i*1.+19.,21.-i); icl->Draw();}
     for(int i=0; i<4; i++) { icl= new TLine(i*1.+18.,22.-i,i*1.+18.,21.-i); icl->Draw();}
     icl= new TLine(22.,18.,22.,15.); icl->Draw();
     icl= new TLine(22.,15.,23.,15.); icl->Draw();
     icl= new TLine(23.,15.,23., 9.); icl->Draw();
     icl= new TLine(23., 9.,22., 9.); icl->Draw();
     icl= new TLine(22., 9.,22., 6.); icl->Draw();
     for(int i=0; i<4; i++) { icl= new TLine(i*1.+18.,3.+i,i*1.+19.,3.+i); icl->Draw();}
     for(int i=0; i<4; i++) { icl= new TLine(i*1.+18.,2.+i,i*1.+18.,3.+i); icl->Draw();}
     icl= new TLine(18., 2.,15., 2.); icl->Draw();
     icl= new TLine(15., 2.,15., 1.); icl->Draw();
     icl= new TLine( 9., 1.,15., 1.); icl->Draw();
     icl= new TLine( 9., 1., 9., 2.); icl->Draw();
     icl= new TLine( 9., 2., 6., 2.); icl->Draw();
     for(int i=0; i<4; i++) { icl= new TLine(i*1.+2.,6.-i,i*1.+3.,6.-i); icl->Draw();}
     for(int i=0; i<4; i++) { icl= new TLine(i*1.+3.,6.-i,i*1.+3.,5.-i); icl->Draw();}
     icl= new TLine( 2., 6., 2., 9.); icl->Draw();
     icl= new TLine( 2., 9., 1., 9.); icl->Draw();

     icl= new TLine( 8.,14., 9.,14.); icl->Draw();
     icl= new TLine( 9.,14., 9.,15.); icl->Draw();
     icl= new TLine( 9.,15.,10.,15.); icl->Draw();
     icl= new TLine(10.,15.,10.,16.); icl->Draw();
     icl= new TLine(10.,16.,14.,16.); icl->Draw();
     icl= new TLine(14.,16.,14.,15.); icl->Draw();
     icl= new TLine(14.,15.,15.,15.); icl->Draw();
     icl= new TLine(15.,15.,15.,14.); icl->Draw();
     icl= new TLine(15.,14.,16.,14.); icl->Draw();
     icl= new TLine(16.,14.,16.,10.); icl->Draw();
     icl= new TLine(16.,10.,15.,10.); icl->Draw();
     icl= new TLine(15.,10.,15., 9.); icl->Draw();
     icl= new TLine(15., 9.,14., 9.); icl->Draw();
     icl= new TLine(14., 9.,14., 9.); icl->Draw();
     icl= new TLine(14., 9.,14., 8.); icl->Draw();
     icl= new TLine(14., 8.,10., 8.); icl->Draw();
     icl= new TLine(10., 8.,10., 9.); icl->Draw();
     icl= new TLine(10., 9., 9., 9.); icl->Draw();
     icl= new TLine( 9., 9., 9.,10.); icl->Draw();
     icl= new TLine( 9.,10., 8.,10.); icl->Draw();
     icl= new TLine( 8.,10., 8.,14.); icl->Draw();
}

