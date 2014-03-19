#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom1.h"
#include "TH1F.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"
#include "../ntupler/trackTree.C"

void makeNtuple(){
 TH1D::SetDefaultSumw2();
 
 //input file
 //TString directory="/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/";
 // TString infname="Dijet100_HydjetDrum_v27_mergedV1";
 TString algo="akVs3Calo";
 TString directory="/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet21_STARTHI53_LV1/merged3/";
 TString infname="HiForest_Pythia_Hydjet_Jet80_Track8_Jet21_STARTHI53_LV1_merged_forest_0"; 
 //TString directory=" /mnt/hadoop/cms/store/user/dgulhan/HIMC/MB/Track8_Jet21_STARTHI53_LV1/merged/";
 //TString infname="HiForest_HYDJET_Track8_Jet21_STARTHI53_LV1_merged_forest_0";  


 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()));
 
 //pt bins for track efficiency correction
 int npt_fake=14; 
 double ptmin_fake[]={0.4,0.4,0.4,0.4,0.4, 1, 1, 1,  1,  1, 3, 3,  3,  8};
 double ptmax_fake[]={  1,  1,  1,  1,  1, 3, 3, 3,  3,  3, 8, 8,  8,300};
 
 int cent_min_fake[]={  0, 20, 40, 60,100, 0,20,40, 60,100, 0,20, 40,  0};
 int cent_max_fake[]={ 20, 40, 60,100,200,20,40,60,100,200,20,40,200,200};
 
 int npt_eff=14; 
 double ptmin_eff[]={0.4,0.4,0.4,0.4,0.4, 1, 1, 1,  1,  1, 3, 3,  3,  8};
 double ptmax_eff[]={  1,  1,  1,  1,  1, 3, 3, 3,  3,  3, 8, 8,  8,300};
 
 int cent_min[]={  0, 20, 40, 60,100, 0,20,40, 60,100, 0,20, 40,  0};
 int cent_max[]={ 20, 40, 60,100,200,20,40,60,100,200,20,40,200,200};
  cout<<0<<endl;

 //getting histograms for track efficiency correction 
 TFile *f_eff[npt_eff];
 TProfile *p_eff_cent[npt_eff]; 
 TProfile2D *p_eff_accept[npt_eff]; 
 TProfile *p_eff_pt[npt_eff]; 
 TProfile *p_eff_rmin[npt_eff]; 
 for(int ipt=0; ipt<npt_eff;ipt++){
   if(ipt<13)f_eff[ipt]= new TFile(Form("../final_hists_vsCalo/eff_pt%d_%d_cent%d_%d_step_cent4accept4pt4rmin3.root",(int)ptmin_eff[ipt],(int)ptmax_eff[ipt],(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   else f_eff[ipt]= new TFile(Form("../final_hists_vsCalo/eff_pt%d_%d_cent%d_%d_step_cent3accept3pt3rmin3.root",(int)ptmin_eff[ipt],(int)ptmax_eff[ipt],(int)(0.5*cent_min[ipt]),(int)(0.5*cent_max[ipt])));
   p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }

 TFile *f_fake[npt_fake];
 TProfile *p_fake_cent[npt_fake]; 
 TProfile2D *p_fake_accept[npt_fake]; 
 TProfile *p_fake_pt[npt_fake]; 
 TProfile *p_fake_rmin[npt_fake]; 
 for(int ipt=0; ipt<npt_fake;ipt++){
   if(ipt==0)f_fake[ipt]= new TFile(Form("../final_hists_vsCalo/fake/fake_pt%d_%d_cent%d_%d_%s_dogenjet0.root",(int)ptmin_fake[ipt],(int)ptmax_fake[ipt],(int)(0.5*cent_min_fake[ipt]),(int)(0.5*cent_max_fake[ipt]),algo.Data()));
   else if(ipt<4)f_fake[ipt]= new TFile(Form("../final_hists_vsCalo/fake/fake_pt%d_%d_cent%d_%d_%s_dogenjet0.root",(int)ptmin_fake[ipt],(int)ptmax_fake[ipt],(int)(0.5*cent_min_fake[ipt]),(int)(0.5*cent_max_fake[ipt]),algo.Data()));
   else if(ipt<13)f_fake[ipt]= new TFile(Form("../final_hists_vsCalo/fake/fake_pt%d_%d_cent%d_%d_%s_dogenjet0.root",(int)ptmin_fake[ipt],(int)ptmax_fake[ipt],(int)(0.5*cent_min_fake[ipt]),(int)(0.5*cent_max_fake[ipt]),algo.Data()));
   else f_fake[ipt]= new TFile(Form("../final_hists_vsCalo/fake/fake_pt%d_%d_cent%d_%d_%s_dogenjet0.root",(int)ptmin_fake[ipt],(int)ptmax_fake[ipt],(int)(0.5*cent_min_fake[ipt]),(int)(0.5*cent_max_fake[ipt]),algo.Data()));
   p_fake_cent[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_cent");
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }
 p_fake_pt[9]->GetBinContent(1);
 //output file and tree
 TFile *outf= new TFile(Form("/export/d00/scratch/abaty/trackingEff/closure_ntuples/track_ntuple_%s_%s_testforfake_small.root",infname.Data(),algo.Data()),"recreate");
 std::string particleVars="pt:matchedpt:eta:phi:rmin:trackselect:cent:eff";

 TNtuple *nt_particle = new TNtuple("nt_particle","",particleVars.data());
 std::string trackVars="pt:eta:phi:rmin:trackselect:trackstatus:cent:eff:trkfake:fake";

 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());

 //loop over events
 int nentries = ftrk->GetEntriesFast();
 // for(int jentry=0;jentry<nentries;jentry++){
 for(int jentry=0;jentry<10000;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fjet->GetEntry(jentry);

  float cent=fhi->hiBin;
  //loop over tracks
  for(int itrk=0;itrk<ftrk->nParticle;itrk++){

   float trackselect=(ftrk->mtrkQual[itrk] && fabs(ftrk->mtrkDxy1[itrk]/ftrk->mtrkDxyError1[itrk])<3.0 && fabs(ftrk->mtrkDz1[itrk]/ftrk->mtrkDzError1[itrk])<3 && (ftrk->mtrkPtError[itrk]/ftrk->mtrkPt[itrk])<0.1);
   float eta=ftrk->pEta[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker
   float pt=ftrk->pPt[itrk];
   // if(pt<0.5) continue;
   // if(pt<0.5 || pt>20) continue; //acceptance of the tracker
   float mpt=ftrk->mtrkPt[itrk];
   float phi=ftrk->pPhi[itrk];
   //float vz = ftrk->vz[itrk];
   //float zVtx = ftrk->zVtx[itrk];
   float rmin=99;
 

   //for(int ijet=0;ijet<fjet->ngen;ijet++){
    //if(fabs(fjet->geneta[ijet])>2 || fjet->genpt[ijet]<30) continue;
    //float r_reco=sqrt(pow(eta-fjet->geneta[ijet],2)+pow(acos(cos(phi-fjet->genphi[ijet])),2));
    //if(r_reco<rmin)rmin=r_reco;
   //}
   for(int ijet=0;ijet<fjet->nref;ijet++){
     if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
     float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
     if(r_reco<rmin)rmin=r_reco;
    }

   
   //get efficiency correction for the track
   float eff_accept=1;
   float eff_pt=1;
   float eff_cent=1;
   float eff_rmin=1;

   for(int ipt=0;ipt<npt_eff;ipt++){
    if(pt>=ptmin_eff[ipt] && pt<ptmax_eff[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));//efficiency for rmin>3 is 1. 
     }     
   } 

   float eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   
   //fill in the output tree

   float entry[]={pt,mpt,eta,phi,rmin,trackselect,cent,eff};

   nt_particle->Fill(entry);

  }
  
  for(int itrk=0;itrk<ftrk->nTrk;itrk++){

   float trackselect=(ftrk->highPurity[itrk] && fabs(ftrk->trkDxy1[itrk]/ftrk->trkDxyError1[itrk])<3.0 && fabs(ftrk->trkDz1[itrk]/ftrk->trkDzError1[itrk])<3 && (ftrk->trkPtError[itrk]/ftrk->trkPt[itrk])<0.1);
   float eta=ftrk->trkEta[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker   
   float pt=ftrk->trkPt[itrk];

 //if(pt<0.5) continue; //acceptance of the tracker
 //acceptance of the tracker
   float phi=ftrk->trkPhi[itrk];
   float trkfake=ftrk->trkFake[itrk];
   float trackstatus=ftrk->trkStatus[itrk];
   //float vz = ftrk->vz[itrk];
   //float zVtx = ftrk->zVtx[itrk];
   float rmin=99;

   //find rmin; 
     for(int ijet=0;ijet<fjet->nref;ijet++){
     if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
     float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
     if(r_reco<rmin)rmin=r_reco;
    }
   
   //for(int ijet=0;ijet<fjet->ngen;ijet++){
    //if(fabs(fjet->geneta[ijet])>2 || fjet->genpt[ijet]<30) continue;
    //float r_reco=sqrt(pow(eta-fjet->geneta[ijet],2)+pow(acos(cos(phi-fjet->genphi[ijet])),2));
    //if(r_reco<rmin)rmin=r_reco;
   //}

   //get efficiency correction for the track
   float eff_accept=1;
   float eff_pt=1;
   float eff_cent=1;
   float eff_rmin=1;
   
   float fake_pt,fake_cent,fake_accept,fake_rmin;
   fake_pt=fake_cent=fake_accept=fake_rmin=0;

   for(int ipt=0;ipt<npt_eff;ipt++){
    if(pt>=ptmin_eff[ipt] && pt<ptmax_eff[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));//efficiency for rmin>3 is 1. 
     }     
   } 
   
   for(int ipt=0;ipt<npt_fake;ipt++){
    if(pt>=ptmin_fake[ipt] && pt<ptmax_fake[ipt] && cent>=cent_min_fake[ipt] && cent<cent_max_fake[ipt]){
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_cent=p_fake_cent[ipt]->GetBinContent(p_fake_cent[ipt]->FindBin(cent));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
     }     
   }

   float eff=1;
   eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   float fake=0;
   if(pt<100)fake=fake_accept+fake_cent+fake_pt+fake_rmin;
   if(eff==0){
    //cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<" cent="<<cent<<endl;
	  if(pt>100)eff=0.8;
	  else eff=1;
   }

   //fill in the output tree
   float entry[]={pt,eta,phi,rmin,trackselect,trackstatus,cent,eff,trkfake,fake};
   nt_track->Fill(entry);

  }
   
 }
 
  //nt_track->Write();
 // nt_particle->Write();
    outf->Write();
    outf->Close();
}
