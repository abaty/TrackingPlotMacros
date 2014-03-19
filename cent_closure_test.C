#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"  
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLine.h"


void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
   if (canv==0) {
      Error("makeMultiPanelCanvas","Got null canvas.");
      return;
   }
   canv->Clear();
   
   TPad* pad[columns][rows];

   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth = 
   (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
   (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
   (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
            Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}
 
void cent_closure_test(){
TH1D::SetDefaultSumw2();
TFile * f= new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/track_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0_akVs3Calo_testforfake_big.root");
TTree * nt_track = (TTree*)f->Get("nt_track");
TTree * nt_particle = (TTree*)f->Get("nt_particle");


TLegend *leg = new TLegend(0.6,0.75,0.95,0.95);
leg->SetBorderSize(0);
leg->SetFillStyle(0);

double bin_pt_min=0.5;
double bin_pt_max=300;
double nbins = 50;
double low   =0;
double high  = 200;


/*
const int ny=50;
double x[ny+1];
double inix=log(bin_pt_min)/log(10);
double delta=(log(bin_pt_max)-log(bin_pt_min))/(50*log(10));
for(int ix=0; ix<ny+1;ix++){
 x[ix]=pow(10,inix); 
 inix+=delta;
}*/

/*
const int ny=50;
double x[ny+1];
double inix=log(bin_pt_min)/log(10);
double delta=(log(bin_pt_max)-log(bin_pt_min))/(50*log(10));
int maxbin=ny;
for(int ix=0; ix<ny+1;ix++){
 x[ix]=pow(10,inix);
 if(x[ix]>100){
        x[ix]=bin_pt_max;
        maxbin=ix;
        break;
 }
 inix+=delta;
}*/

//eff correction
TH1D * h_gen = new  TH1D("h_gen",";cent;N_{evt}",nbins,low,high);
TH1D * h_gen_select = new  TH1D("h_gen_select",";cent;N_{evt}",nbins,low,high);
TH1D * h_gen_matched_select_corr = new TH1D("h_gen_matched_select_corr",";cent;N_{evt}",nbins,low,high);

nt_particle->Draw("cent>>h_gen","pt>0.5");
nt_particle->Draw("cent>>h_gen_select","(trackselect && pt>0.5)");
nt_particle->Draw("cent>>h_gen_matched_select_corr","(1/eff)*(trackselect && pt>0.5)");

h_gen->SetMarkerColor(1);
h_gen->SetMarkerStyle(25);
h_gen->SetLineWidth(1);
h_gen_select->SetMarkerColor(1);
h_gen_matched_select_corr->SetMarkerColor(kRed);
h_gen_matched_select_corr->SetLineColor(kRed);

TLegend *leg2 = new TLegend(0.3,0.15,0.75,0.55);
leg2->SetBorderSize(0);
leg2->SetFillStyle(0);

leg2->AddEntry(h_gen,"gen","p");
leg2->AddEntry(h_gen_select,"gen select matched uncorr","p");
leg2->AddEntry(h_gen_matched_select_corr,"gen select matched corr","p");

TCanvas * c2 = new TCanvas("c2","",600,800);
makeMultiPanelCanvas(c2,1,2,0.0,0.0,0.15,0.15,0.02);
c2->cd(1);
//c2->cd(1)->SetLogx();
c2->cd(1)->SetLogy();
h_gen->Draw();
//h_gen->Draw("same hist");
h_gen_select->Draw("same");
h_gen_matched_select_corr->Draw("same");
leg2->Draw("same");
c2->cd(2);
//c2->cd(2)->SetLogx();

TH1D * hgen_corr_rat = (TH1D*)h_gen_matched_select_corr->Clone("hgen_corr_rat");
hgen_corr_rat->Divide(h_gen);
//hgen_corr_rat->GetYaxis()->SetTitle("N_{gen}/N_{gen,sel,match}");

TH1D * hgen_uncorr_rat = (TH1D*)h_gen_select->Clone("hgen_uncorr_rat");
hgen_uncorr_rat->Divide(h_gen);
hgen_uncorr_rat->SetMaximum(1.1);
hgen_uncorr_rat->SetMinimum(0.2);

hgen_uncorr_rat->GetYaxis()->SetTitle("N_{gen,sel,match}/N_{gen}");
//hgen_corr_rat->SetMaximum(1.1);
//hgen_corr_rat->SetMinimum(0.7);
hgen_uncorr_rat->Draw();
hgen_corr_rat->Draw("same");

TLegend *leg2a = new TLegend(0.7,0.2,0.95,0.4);
leg2a->SetBorderSize(0);
leg2a->SetFillStyle(0);

leg2a->AddEntry(hgen_corr_rat,"corr","p");
leg2a->AddEntry(hgen_uncorr_rat,"uncorr","p");
//leg2a->Draw("same");

c2->SaveAs("compare_gen_select_corr_cent.png");
c2->SaveAs("compare_gen_select_corr_cent.pdf");

//fake correction
TH1D * h_reco = new  TH1D("h_reco",";cent;N_{evt}",nbins,low,high);
TH1D * h_reco_matched = new  TH1D("h_reco_matched",";cent;N_{evt}",nbins,low,high);
TH1D * h_reco_fakecorr = new TH1D("h_reco_fakecorr",";cent;N_{evt}",nbins,low,high);
nt_track->Draw("cent>>h_reco","(trackselect && pt>0.5)"); 
nt_track->Draw("cent>>h_reco_matched","(trackselect && !trkfake && pt>0.5)"); 
nt_track->Draw("cent>>h_reco_fakecorr","(1-fake)*(trackselect && pt>0.5)"); 

//h_reco->SetMarkerSize(1);
//h_reco->SetMarkerStyle(25);
h_reco_matched->SetMarkerColor(1);
h_reco_matched->SetMarkerStyle(25);
h_reco_fakecorr->SetMarkerColor(kRed);
h_reco_fakecorr->SetLineColor(kRed);

TLegend *leg3 = new TLegend(0.3,0.15,0.75,0.55);
leg3->SetBorderSize(0);
leg3->SetFillStyle(0);

leg3->AddEntry(h_reco,"reco","p");
leg3->AddEntry(h_reco_matched,"reco matched","p");
leg3->AddEntry(h_reco_fakecorr,"reco fake corr","p");

TCanvas * c3 = new TCanvas("c3","",600,800);
makeMultiPanelCanvas(c3,1,2,0.0,0.0,0.15,0.15,0.02);
c3->cd(1);
//c3->cd(1)->SetLogx();
c3->cd(1)->SetLogy();
h_reco->Draw();
//h_reco->Draw("same hist");
h_reco_matched->Draw("same");
h_reco_fakecorr->Draw("same");
leg3->Draw("same");
c3->cd(2);
//c3->cd(2)->SetLogx();
TH1D * hreco_fakecorr_rat = (TH1D*)h_reco_fakecorr->Clone("hreco_fakecorr_rat");
hreco_fakecorr_rat->Divide(h_reco_matched);
//hreco_fakecorr_rat->GetYaxis()->SetTitle("N_{reco,matched}/N_{reco,fake}");

TH1D * hreco_fakeuncorr_rat = (TH1D*)h_reco->Clone("hreco_fakeuncorr_rat");
hreco_fakeuncorr_rat->Divide(h_reco_matched);
hreco_fakeuncorr_rat->SetMaximum(1.4);
hreco_fakeuncorr_rat->SetMinimum(0.9);
hreco_fakeuncorr_rat->GetYaxis()->SetTitle("N_{reco}/N_{reco,matched}");
//hreco_fakecorr_rat->SetMaximum(1.1);
//hreco_fakecorr_rat->SetMinimum(0.85);
hreco_fakeuncorr_rat->Draw();
hreco_fakecorr_rat->Draw("same");

TLegend *leg3a = new TLegend(0.5,0.7,0.95,0.9);
leg3a->SetBorderSize(0);
leg3a->SetFillStyle(0);

leg3a->AddEntry(hreco_fakecorr_rat,"corr","p");
leg3a->AddEntry(hreco_fakeuncorr_rat,"uncorr","p");
//leg3a->Draw("same");

c3->SaveAs("compare_reco_fake_corr_cent.png");
c3->SaveAs("compare_reco_fake_corr_cent.pdf");

//full correction
TH1D * h_reco_fakecorr_effcorr = new TH1D("h_reco_fakecorr_effcorr",";cent;N_{evt}",nbins,low,high);
TH1D * h_reco_fakecorr_effcorr_unfold = new TH1D("h_reco_fakecorr_effcorr_unfold",";cent;N_{evt}",nbins,low,high);
nt_track->Draw("cent>>h_reco_fakecorr_effcorr","((1-fake)/eff)*(trackselect && pt>0.5)"); 

h_reco_fakecorr_effcorr->SetMarkerColor(kRed);
h_reco_fakecorr_effcorr->SetLineColor(kRed);
TLegend *leg4 = new TLegend(0.3,0.15,0.75,0.55);
leg4->SetBorderSize(0);
leg4->SetFillStyle(0);

leg4->AddEntry(h_gen,"gen","p");
leg4->AddEntry(h_reco_fakecorr_effcorr,"reco corr","p");
leg4->AddEntry(h_reco,"reco","p");


TCanvas *c4 = new TCanvas("c4","",600,800);
makeMultiPanelCanvas(c4,1,2,0.0,0.0,0.15,0.15,0.02);

c4->cd(1);
//c4->cd(1)->SetLogx();
c4->cd(1)->SetLogy();
h_gen->Draw();
h_reco->Draw("same");
h_reco_fakecorr_effcorr->Draw("same");
leg4->Draw("same");
c4->cd(2);
//c4->cd(2)->SetLogx(); 
TH1D * h_genreco_fullcorr=(TH1D*)h_reco_fakecorr_effcorr->Clone("h_genreco_fullcorr");
h_genreco_fullcorr->Divide(h_gen);
//h_genreco_fullcorr->GetYaxis()->SetTitle("N_{reco,select,corr}/N_{gen}");
//h_genreco_fullcorr->SetMaximum(1.2);
//h_genreco_fullcorr->SetMinimum(0.8);

TH1D * h_genreco=(TH1D*)h_reco->Clone("h_genreco");
h_genreco->Divide(h_gen);
h_genreco->GetYaxis()->SetTitle("N_{reco,select,corr}/N_{gen}");
h_genreco->SetMaximum(1.1);
h_genreco->SetMinimum(0.2);
h_genreco->Draw();
h_genreco_fullcorr->Draw("same");

TLegend *leg4a = new TLegend(0.7,0.2,0.95,0.4);
leg4a->SetBorderSize(0);
leg4a->SetFillStyle(0);

leg4a->AddEntry(h_genreco_fullcorr,"corr","p");
leg4a->AddEntry(h_genreco,"uncorr","p");
//leg4a->Draw("same");

c4->SaveAs("compare_select_fullcorr_cent.png");
c4->SaveAs("compare_select_fullcorr_cent.pdf");

//this part is for crosscheck to see if the matched particle spectrum agrees with spectrum of not-fake tracks
/*
TH1D *h_matched_pt_particle= new TH1D("h_matched_pt_particle","",nbins,low,high);
TH1D *h_matched_pt_trk= new TH1D("h_matched_pt_trk","",nbins,low,high);
nt_particle->Draw("matchedpt>>h_matched_pt_particle","");
nt_track->Draw("pt>>h_matched_pt_trk","(trkfake!=1)");
h_matched_pt_trk->SetMarkerColor(kRed);
h_matched_pt_trk->SetLineColor(kRed);
h_matched_pt_particle->SetMarkerStyle(25);

TCanvas *c5 = new TCanvas("c5","",600,800);
makeMultiPanelCanvas(c5,1,2,0.0,0.0,0.15,0.15,0.02);

c5->cd(1);
c5->cd(1)->SetLogx();
c5->cd(1)->SetLogy();
h_matched_pt_trk->Draw();
h_matched_pt_particle->Draw("same");
c5->cd(2); 
c5->cd(2)->SetLogx();
TH1D * h_reco_gen_trk=(TH1D*)h_matched_pt_trk->Clone("h_reco_gen_trk");
h_reco_gen_trk->Divide(h_matched_pt_particle);
h_reco_gen_trk->GetYaxis()->SetTitle("N_{reco,select}/N_{gen,matched,select}");
h_reco_gen_trk->SetMaximum(1.1);
h_reco_gen_trk->SetMinimum(0.85);
h_reco_gen_trk->Draw();
c5->SaveAs("compare_matchedpt_recopt.png");
*/
}
