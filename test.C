#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TFile.h>
#include <TChain.h>

using namespace std;

map<TString, TH1*> maphist;

/// Making Histogram Function //////////////////////////////////////////////// 
void MakeHistograms(TString hname, int nbins, float xmin, float xmax){

  maphist[hname] = new TH1F(hname.Data(), hname.Data(), nbins, xmin, xmax);

}
//////////////////////////////////////////////////////////////////////////////



/// Getting Histogram Function /////////////////////////////////////////////// 
TH1 * GetHist(TString hname){

  TH1 *h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit-> second;

  return h;

}
//////////////////////////////////////////////////////////////////////////////



/// Filling Histogram Function ///////////////////////////////////////////////
void FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins){

  if(GetHist(histname)) GetHist(histname) -> Fill(value, w);
  else{
    cout << "Making histogram..." <<endl;
    MakeHistograms(histname, nbins, xmin, xmax);
    if(GetHist(histname)) GetHist(histname) -> Fill(value, w);
  }

}
//////////////////////////////////////////////////////////////////////////////



/// Calculating delta Phi Function ///////////////////////////////////////////    
float delta_phi( float phi1, float phi2){
  float sub_phi = phi2 - phi1;
  float del_phi;
  if( cos(sub_phi / 2) >= 0){
    return sub_phi;
  }
  else if( sub_phi < 0){
    del_phi = sub_phi + 2 * TMath::Pi();
    return del_phi;
  }
  else{
    del_phi = sub_phi - 2 * TMath::Pi();
    return del_phi;
  }
}
//////////////////////////////////////////////////////////////////////////////




void test::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L test.C
//      Root > test t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //float EgEt, EgEta, EgPhi, EgN;
   //float GenElPT, GenElPhi;
   
   TH2F *EtEt = new TH2F("genEt & EgEt", "GenEt & EgEt", 80, 0, 80, 80, 0, 80);
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   //vectors to use to draw histograms
   TVector3 L1_hit;
   TVector3 L2_hit;
   TVector3 L3_hit;
   TVector3 L4_hit;
   TVector3 D1_hit;
   TVector3 D2_hit;
   TVector3 D3_hit;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry < nentries ; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
     //Vectors to pushback
     std::vector<TVector3> L1_hits;
     std::vector<TVector3> L2_hits;
     std::vector<TVector3> L3_hits;
     std::vector<TVector3> L4_hits;
     std::vector<TVector3> D1_hits;
     std::vector<TVector3> D2_hits;
     std::vector<TVector3> D3_hits;
     std::vector<int> hitted_layers;
     
     
     float EgN = egCrysN;
     float EgEt, EgEta, EgPhi, EgGx, EgGy, EgGz;
     float GenElPT, GenElPhi, GenElEta;;
     int junk;
     //hit number for each pixel layer
     int L1_N = 0;
     int L2_N = 0;
     int L3_N = 0;
     int L4_N = 0;
     int D1_N = 0;
     int D2_N = 0;
     int D3_N = 0;
     
     GenElPT = genPartPt->at(0);
     GenElPhi = genPartPhi->at(0);
     GenElEta = genPartEta->at(0);

     //if(GenElPT > 140) continue;
     
     vector<float> EGET;
     
     float mini = 999.;
     for(int k = 0; k < EgN; k++){
       EgEt = egCrysEt -> at(k);
       if(1){
	 EGET.push_back(EgEt);
	 sort(EGET.begin(), EGET.end(), greater<float>());
	 if(mini >= fabs(GenElPT - EGET[k])){
	   mini = fabs(GenElPT - EGET[k]);
	 }
       }
     }
     
     if(egCrysN != 1) continue;
     
     for(int q = 0; q < egCrysN; q++){
       
       EgEt = egCrysEt -> at(q);
       EgEta = egCrysEta -> at(q);
       EgPhi = egCrysPhi -> at(q);

       EgGx = egCrysGx -> at(q);
       EgGy = egCrysGy -> at(q);
       EgGz = egCrysGz -> at(q);
       
       
       //EgEt cut
       if(EgEt < 10) continue;
       /*
       if(EgEt>= 9 && EgEt<11 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.4)  /(-0.02-0.02)  *(GenElPhi-EgPhi +0.02  )  -0.2 ) continue; }
       else if(EgEt>=11 && EgEt<15 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.4)  /(-0.02-0.02)  *(GenElPhi-EgPhi +0.02  )  -0.2 ) continue; }
       else if(EgEt>=15 && EgEt<20 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.17) / -0.011       *(GenElPhi-EgPhi        )  -0.2 ) continue; }
       else if(EgEt>=20 && EgEt<30 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.17) /(0.018-0.026) *(GenElPhi-EgPhi -0.018 )  -0.2 ) continue; }
       else if(EgEt>=30 && EgEt<40 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.06) /(0.015-0.036) *(GenElPhi-EgPhi -0.015 )  -0.2 ) continue; }
       else if(EgEt>=40 && EgEt<50 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.03) /(0.020-0.042) *(GenElPhi-EgPhi -0.020 )  -0.2 ) continue; }
       else if(EgEt>=50 ) { if( ((GenElPT-EGET[q])/GenElPT) < -0.14 ) continue; }       
       */

       float ratio = EgEt / GenElPT ;
       //if( fabs(ratio - 1) > 0.5 ) continue;
       
       EtEt -> Fill( GenElPT, EgEt );

       float match_R = sqrt( (EgEta - GenElEta) * (EgEta - GenElEta) + (EgPhi - GenElPhi) * (EgPhi -GenElPhi) );

       if(match_R > 0.3) continue;
 
       //to use only R4 region
       //if( fabs(EgEta) < 1.9 || fabs(EgEta) > 2.5 ) continue;
       
       
       for(int bhit = 0; bhit < bHitN; bhit ++){
	 
	 float b_cl_r = sqrt( pow(bHitGx->at(bhit), 2) + pow(bHitGy->at(bhit), 2) + pow(bHitGz->at(bhit),2) );
	 
	 //if(bHitGz->at(bhit) > 28) continue; //z-length cut
	 
	 if(b_cl_r > 2 && b_cl_r < 5 && fabs(bHitGz->at(bhit)) < 28){ 
	   L1_hits.push_back( TVector3(bHitGx->at(bhit), bHitGy->at(bhit), bHitGz->at(bhit)));
	   L1_N ++;
	 }
	 else if(b_cl_r > 6 && b_cl_r < 8 && fabs(bHitGz->at(bhit)) < 28){
	   L2_hits.push_back( TVector3(bHitGx->at(bhit), bHitGy->at(bhit), bHitGz->at(bhit)));
           L2_N ++;
	 }
	 else if(b_cl_r> 10 && b_cl_r < 12 && fabs(bHitGz->at(bhit)) < 28){
	   L3_hits.push_back( TVector3(bHitGx->at(bhit), bHitGy->at(bhit), bHitGz->at(bhit)));
           L3_N ++;
	 }
	 else if(b_cl_r> 14 && b_cl_r < 18 && fabs(bHitGz->at(bhit)) < 28){
	   L4_hits.push_back( TVector3(bHitGx->at(bhit), bHitGy->at(bhit), bHitGz->at(bhit)));
           L4_N ++;
	 }

	 
       }//bhit for
       
       for(int fhit = 0; fhit < fHitN; fhit ++){
	 
	 float f_z = fHitGz->at(fhit);
	 
	 if( fabs(f_z) > 28 && fabs(f_z) < 36){
	   D1_hits.push_back( TVector3(fHitGx->at(fhit), fHitGy->at(fhit), fHitGz->at(fhit)));
	   D1_N ++;
	 }
	 else if( fabs(f_z) > 36 && fabs(f_z) < 44){
	   D2_hits.push_back( TVector3(fHitGx->at(fhit), fHitGy->at(fhit), fHitGz->at(fhit)));
	   D2_N ++;
	 }
	 else if( fabs(f_z) > 44 && fabs(f_z) < 54){
	   D3_hits.push_back( TVector3(fHitGx->at(fhit), fHitGy->at(fhit), fHitGz->at(fhit)));
	   D3_N ++;
	 }
	 else junk = 1;
	 
       }//fhit for
       
       
       for(int range = 10; range < 141; range ++){ 
	 if(EgEt >= range && EgEt < range+1 ){
	   char range_in_char[30];//TString range_in_string;   
	   int range_sub = range;
	   sprintf( range_in_char, "%d", range_sub);
	   TString range_in_string = range_in_char;
	   if(range <10){
	     range_in_string.Insert(0, "00");
	   }
	   if(range >9 && range < 100){
	     range_in_string.Insert(0, "0");
	   }
	   TString histname_phi = "_layer_delta_phi_";
	   TString histname_eta = "_layer_delta_eta_";
	   TString histname_R = "_layer_delta_R_";
	   
	   histname_phi.Append(range_in_string);
	   histname_eta.Append(range_in_string);
	   histname_R.Append(range_in_string);

	   TString L1_string = "L1_";
	   TString L2_string = "L2_";
	   TString L3_string = "L3_";
	   TString L4_string = "L4_";
	   TString D1_string = "D1_";
	   TString D2_string = "D2_";
	   TString D3_string = "D3_";

	   
	   //Region 1, L1234
	   if( fabs(EgEta) < 1.3 ){
	     /*
	     L1_string.Insert(0, "R1_");
	     L2_string.Insert(0, "R1_");
	     L3_string.Insert(0, "R1_");
	     L4_string.Insert(0, "R1_");
	     */
	     if(L1_N > 0){
	       for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
		 L1_hit.SetXYZ( L1_hits[hitnum_L1].X(), L1_hits[hitnum_L1].Y(), L1_hits[hitnum_L1].Z() );
		 FillHist(L1_string + histname_phi, delta_phi( EgPhi, L1_hit.Phi() ), 1., -0.3, 0.3, 2400);

		 TVector3 current_vector;
		 current_vector.SetXYZ( EgGx - L1_hits[hitnum_L1].X(), EgGy - L1_hits[hitnum_L1].Y(), EgGz - L1_hits[hitnum_L1].Z());
		 FillHist(L1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(L2_N > 0){
	       for(int hitnum_L2 = 0; hitnum_L2 <L2_N; hitnum_L2++){
		 L2_hit.SetXYZ( L2_hits[hitnum_L2].X(), L2_hits[hitnum_L2].Y(), L2_hits[hitnum_L2].Z() );
		 FillHist(L2_string + histname_phi, delta_phi( EgPhi, L2_hit.Phi() ), 1., -0.3, 0.3, 2400);

		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L2_hits[hitnum_L2].X(), EgGy - L2_hits[hitnum_L2].Y(), EgGz - L2_hits[hitnum_L2].Z());
                 FillHist(L2_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(L3_N > 0){
	       for(int hitnum_L3 = 0; hitnum_L3 <L3_N; hitnum_L3++){
		 L3_hit.SetXYZ( L3_hits[hitnum_L3].X(), L3_hits[hitnum_L3].Y(), L3_hits[hitnum_L3].Z() );
		 FillHist(L3_string + histname_phi, delta_phi( EgPhi, L3_hit.Phi() ), 1., -0.3, 0.3, 2400);
	       
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L3_hits[hitnum_L3].X(), EgGy - L3_hits[hitnum_L3].Y(), EgGz - L3_hits[hitnum_L3].Z());
                 FillHist(L3_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(L4_N > 0){
	       for(int hitnum_L4 = 0; hitnum_L4 <L4_N; hitnum_L4++){
		 L4_hit.SetXYZ( L4_hits[hitnum_L4].X(), L4_hits[hitnum_L4].Y(), L4_hits[hitnum_L4].Z() );
		 FillHist(L4_string + histname_phi, delta_phi( EgPhi, L4_hit.Phi() ), 1., -0.3, 0.3, 2400);

		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L4_hits[hitnum_L4].X(), EgGy - L4_hits[hitnum_L4].Y(), EgGz - L4_hits[hitnum_L4].Z());
                 FillHist(L4_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	   }
	   //Region 2, L123D1
	   else if( fabs(EgEta) > 1.3 && fabs(EgEta) < 1.6 ){
	     /*
	     L1_string.Insert(0, "R2_");
             L2_string.Insert(0, "R2_");
             L3_string.Insert(0, "R2_");
             D1_string.Insert(0, "R2_");
	     */
	     if(L1_N > 0){
	       for(int hitnum_L1 = 0; hitnum_L1 <L1_N; hitnum_L1++){
		 L1_hit.SetXYZ( L1_hits[hitnum_L1].X(), L1_hits[hitnum_L1].Y(), L1_hits[hitnum_L1].Z() );
		 FillHist(L1_string + histname_phi, delta_phi( EgPhi, L1_hit.Phi() ), 1., -0.3, 0.3, 2400);
	       
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L1_hits[hitnum_L1].X(), EgGy - L1_hits[hitnum_L1].Y(), EgGz - L1_hits[hitnum_L1].Z());
                 FillHist(L1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(L2_N > 0){
	       for(int hitnum_L2 = 0; hitnum_L2 <L2_N; hitnum_L2++){
		 L2_hit.SetXYZ( L2_hits[hitnum_L2].X(), L2_hits[hitnum_L2].Y(), L2_hits[hitnum_L2].Z() );
		 FillHist(L2_string + histname_phi, delta_phi( EgPhi, L2_hit.Phi() ), 1., -0.3, 0.3, 2400);
	       
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L2_hits[hitnum_L2].X(), EgGy - L2_hits[hitnum_L2].Y(), EgGz - L2_hits[hitnum_L2].Z());
                 FillHist(L2_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(L3_N > 0){
	       for(int hitnum_L3 = 0; hitnum_L3 <L3_N; hitnum_L3++){
		 L3_hit.SetXYZ( L3_hits[hitnum_L3].X(), L3_hits[hitnum_L3].Y(), L3_hits[hitnum_L3].Z() );
		 FillHist(L3_string + histname_phi, delta_phi( EgPhi, L3_hit.Phi() ), 1., -0.3, 0.3, 2400);
	
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L3_hits[hitnum_L3].X(), EgGy - L3_hits[hitnum_L3].Y(), EgGz - L3_hits[hitnum_L3].Z());
                 FillHist(L3_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D1_N > 0){
	       for(int hitnum_D1 = 0; hitnum_D1 <D1_N; hitnum_D1++){
		 D1_hit.SetXYZ( D1_hits[hitnum_D1].X(), D1_hits[hitnum_D1].Y(), D1_hits[hitnum_D1].Z() );
		 FillHist(D1_string + histname_phi, delta_phi( EgPhi, D1_hit.Phi() ), 1., -0.3, 0.3, 2400);
	
		 
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D1_hits[hitnum_D1].X(), EgGy - D1_hits[hitnum_D1].Y(), EgGz - D1_hits[hitnum_D1].Z());
                 FillHist(D1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	   }
	   //Region 3, L12D12
	   else if( fabs(EgEta) > 1.6 && fabs(EgEta) < 1.9 ){
	     /*
	     L1_string.Insert(0, "R3_");
             L2_string.Insert(0, "R3_");
             D1_string.Insert(0, "R3_");
             D2_string.Insert(0, "R3_");
	     */
	     if(L1_N > 0){
	       for(int hitnum_L1 = 0; hitnum_L1 <L1_N; hitnum_L1++){
		 L1_hit.SetXYZ( L1_hits[hitnum_L1].X(), L1_hits[hitnum_L1].Y(), L1_hits[hitnum_L1].Z() );
		 FillHist(L1_string + histname_phi, delta_phi( EgPhi, L1_hit.Phi() ), 1., -0.3, 0.3, 2400);
	    
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L1_hits[hitnum_L1].X(), EgGy - L1_hits[hitnum_L1].Y(), EgGz - L1_hits[hitnum_L1].Z());
                 FillHist(L1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(L2_N > 0){
	       for(int hitnum_L2 = 0; hitnum_L2 <L2_N; hitnum_L2++){
		 L2_hit.SetXYZ( L2_hits[hitnum_L2].X(), L2_hits[hitnum_L2].Y(), L2_hits[hitnum_L2].Z() );
		 FillHist(L2_string + histname_phi, delta_phi( EgPhi, L2_hit.Phi() ), 1., -0.3, 0.3, 2400);
	       
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L2_hits[hitnum_L2].X(), EgGy - L2_hits[hitnum_L2].Y(), EgGz - L2_hits[hitnum_L2].Z());
                 FillHist(L2_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D1_N > 0){
	       for(int hitnum_D1 = 0; hitnum_D1 <D1_N; hitnum_D1++){
		 D1_hit.SetXYZ( D1_hits[hitnum_D1].X(), D1_hits[hitnum_D1].Y(), D1_hits[hitnum_D1].Z() );
		 FillHist(D1_string + histname_phi, delta_phi( EgPhi, D1_hit.Phi() ), 1., -0.3, 0.3, 2400);

		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D1_hits[hitnum_D1].X(), EgGy - D1_hits[hitnum_D1].Y(), EgGz - D1_hits[hitnum_D1].Z());
                 FillHist(D1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D2_N > 0){
	       for(int hitnum_D2 = 0; hitnum_D2 <D2_N; hitnum_D2++){
		 D2_hit.SetXYZ( D2_hits[hitnum_D2].X(), D2_hits[hitnum_D2].Y(), D2_hits[hitnum_D2].Z() );
		 FillHist(D2_string + histname_phi, delta_phi( EgPhi, D2_hit.Phi() ), 1., -0.3, 0.3, 2400);
	  
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D2_hits[hitnum_D2].X(), EgGy - D2_hits[hitnum_D2].Y(), EgGz - D2_hits[hitnum_D2].Z());
                 FillHist(D2_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	   }
	   //Region 4, L1D123
	   else if( fabs(EgEta) > 1.9 && fabs(EgEta) < 2.5 ){
	     /*
	     L1_string.Insert(0, "R4_");
             D1_string.Insert(0, "R4_");
             D2_string.Insert(0, "R4_");
             D3_string.Insert(0, "R4_");
	     */
	     if(L1_N > 0){
	       for(int hitnum_L1 = 0; hitnum_L1 <L1_N; hitnum_L1++){
		 L1_hit.SetXYZ( L1_hits[hitnum_L1].X(), L1_hits[hitnum_L1].Y(), L1_hits[hitnum_L1].Z() );
		 FillHist(L1_string + histname_phi, delta_phi( EgPhi, L1_hit.Phi() ), 1., -0.3, 0.3, 2400);
	    
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L1_hits[hitnum_L1].X(), EgGy - L1_hits[hitnum_L1].Y(), EgGz - L1_hits[hitnum_L1].Z());
                 FillHist(L1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D1_N > 0){
	       for(int hitnum_D1 = 0; hitnum_D1 <D1_N; hitnum_D1++){
		 D1_hit.SetXYZ( D1_hits[hitnum_D1].X(), D1_hits[hitnum_D1].Y(), D1_hits[hitnum_D1].Z() );
		 FillHist(L2_string + histname_phi, delta_phi( EgPhi, L2_hit.Phi() ), 1., -0.3, 0.3, 2400);

		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D1_hits[hitnum_D1].X(), EgGy - D1_hits[hitnum_D1].Y(), EgGz - D1_hits[hitnum_D1].Z());
                 FillHist(D1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D2_N > 0){
	       for(int hitnum_D2 = 0; hitnum_D2 <L3_N; hitnum_D2++){
		 D2_hit.SetXYZ( D2_hits[hitnum_D2].X(), D2_hits[hitnum_D2].Y(), D2_hits[hitnum_D2].Z() );
		 FillHist(L3_string + histname_phi, delta_phi( EgPhi, L3_hit.Phi() ), 1., -0.3, 0.3, 2400);
		 
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D2_hits[hitnum_D2].X(), EgGy - D2_hits[hitnum_D2].Y(), EgGz - D2_hits[hitnum_D2].Z());
                 FillHist(D2_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D3_N > 0){
	       for(int hitnum_D3 = 0; hitnum_D3 <D3_N; hitnum_D3++){
		 D3_hit.SetXYZ( D3_hits[hitnum_D3].X(), D3_hits[hitnum_D3].Y(), D3_hits[hitnum_D3].Z() );
		 FillHist(D3_string + histname_phi, delta_phi( EgPhi, D3_hit.Phi() ), 1., -0.3, 0.3, 2400);
	       
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D3_hits[hitnum_D3].X(), EgGy - D3_hits[hitnum_D3].Y(), EgGz - D3_hits[hitnum_D3].Z());
                 FillHist(D3_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	   }
	   //Region 5, L1D123
	   else if( fabs(EgEta) > 2.5 && fabs(EgEta) < 2.8 ){
	     /*
	     L1_string.Insert(0, "R5_");
             D1_string.Insert(0, "R5_");
             D2_string.Insert(0, "R5_");
             D3_string.Insert(0, "R5_");
	     */
	     if(L1_N > 0){
	       for(int hitnum_L1 = 0; hitnum_L1 <L1_N; hitnum_L1++){
		 L1_hit.SetXYZ( L1_hits[hitnum_L1].X(), L1_hits[hitnum_L1].Y(), L1_hits[hitnum_L1].Z() );
		 FillHist(L1_string + histname_phi, delta_phi( EgPhi, L1_hit.Phi() ), 1., -0.3, 0.3, 2400);

		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - L1_hits[hitnum_L1].X(), EgGy - L1_hits[hitnum_L1].Y(), EgGz - L1_hits[hitnum_L1].Z());
                 FillHist(L1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D1_N > 0){
	       for(int hitnum_D1 = 0; hitnum_D1 <D1_N; hitnum_D1++){
		 D1_hit.SetXYZ( D1_hits[hitnum_D1].X(), D1_hits[hitnum_D1].Y(), D1_hits[hitnum_D1].Z() );
		 FillHist(L2_string + histname_phi, delta_phi( EgPhi, L2_hit.Phi() ), 1., -0.3, 0.3, 2400);
		 
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D1_hits[hitnum_D1].X(), EgGy - D1_hits[hitnum_D1].Y(), EgGz - D1_hits[hitnum_D1].Z());
                 FillHist(D1_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D2_N > 0){
	       for(int hitnum_D2 = 0; hitnum_D2 <L3_N; hitnum_D2++){
		 D2_hit.SetXYZ( D2_hits[hitnum_D2].X(), D2_hits[hitnum_D2].Y(), D2_hits[hitnum_D2].Z() );
		 FillHist(L3_string + histname_phi, delta_phi( EgPhi, L3_hit.Phi() ), 1., -0.3, 0.3, 2400);
	     
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D2_hits[hitnum_D2].X(), EgGy - D2_hits[hitnum_D2].Y(), EgGz - D2_hits[hitnum_D2].Z());
                 FillHist(D2_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	     if(D3_N > 0){
	       for(int hitnum_D3 = 0; hitnum_D3 <D3_N; hitnum_D3++){
		 D3_hit.SetXYZ( D3_hits[hitnum_D3].X(), D3_hits[hitnum_D3].Y(), D3_hits[hitnum_D3].Z() );
		 FillHist(D3_string + histname_phi, delta_phi( EgPhi, D3_hit.Phi() ), 1., -0.3, 0.3, 2400);
	  
		 TVector3 current_vector;
                 current_vector.SetXYZ( EgGx - D3_hits[hitnum_D3].X(), EgGy - D3_hits[hitnum_D3].Y(), EgGz - D3_hits[hitnum_D3].Z());
                 FillHist(D3_string + histname_eta, current_vector.Eta() - EgEta, 1., -0.3, 0.3, 2400);
	       }
	     }
	   }

	   
	 }//if EgEt range
	 
	 
       }//for range (Pt)
       
       
     }//EgN for
     
     
     
     
     
     
     
     
     if( jentry % 10000 == 0){
       cout << jentry << " / " << nentries << "  done" << endl;
     }
     
     
     
     // if (Cut(ientry) < 0) continue;
   }//end of nentries loop
   
   TFile *file = new TFile("First_result.root", "recreate");
   for(map<TString, TH1*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit ++){
     mapit->second->Write();
   }
   
   cout << "job is finished  "<< endl;
 

   TCanvas * c = new TCanvas("", "", 800, 600);
   c->cd();
   EtEt->Draw();
   c->Update();
   c->SaveAs("./EtEt.pdf");


  
}
