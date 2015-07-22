//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 10 12:54:36 2015 by ROOT version 5.34/14
// from TTree t/t
// found on file: Single_Ele_Eta3_with_Single_Crystal.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           bunchN;
   vector<int>     *pileup;
   Int_t           genPartN;
   vector<float>   *genPartE;
   vector<float>   *genPartPt;
   vector<float>   *genPartEta;
   vector<float>   *genPartPhi;
   vector<int>     *genPartCharge;
   vector<int>     *genPartId;
   Int_t           simTrkN;
   vector<float>   *simTrkPt;
   vector<float>   *simTrkEta;
   vector<float>   *simTrkPhi;
   vector<int>     *simTrkId;
   vector<int>     *simTrkType;
   vector<float>   *simTrkVx;
   vector<float>   *simTrkVy;
   vector<float>   *simTrkVz;
   vector<float>   *simVx;
   vector<float>   *simVy;
   vector<float>   *simVz;
   Int_t           egCrysN;
   vector<float>   *egCrysE;
   vector<float>   *egCrysEt;
   vector<float>   *egCrysEta;
   vector<float>   *egCrysPhi;
   vector<float>   *egCrysGx;
   vector<float>   *egCrysGy;
   vector<float>   *egCrysGz;
   vector<int>     *egCrysCharge;
   Int_t           egN;
   vector<float>   *egE;
   vector<float>   *egEt;
   vector<float>   *egEta;
   vector<float>   *egPhi;
   vector<float>   *egGx;
   vector<float>   *egGy;
   vector<float>   *egGz;
   Int_t           fHitN;
   vector<int>     *fHitDisk;
   vector<int>     *fHitBlade;
   vector<int>     *fHitSide;
   vector<float>   *fHitGx;
   vector<float>   *fHitGy;
   vector<float>   *fHitGz;
   vector<int>     *fClSize;
   vector<int>     *fClSizeX;
   vector<int>     *fClSizeY;
   Int_t           bHitN;
   vector<int>     *bHitLayer;
   vector<int>     *bHitLadder;
   vector<float>   *bHitGx;
   vector<float>   *bHitGy;
   vector<float>   *bHitGz;
   vector<int>     *bClSize;
   vector<int>     *bClSizeX;
   vector<int>     *bClSizeY;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchN;   //!
   TBranch        *b_pileup;   //!
   TBranch        *b_genPartN;   //!
   TBranch        *b_genPartE;   //!
   TBranch        *b_genPartPt;   //!
   TBranch        *b_genPartEta;   //!
   TBranch        *b_genPartPhi;   //!
   TBranch        *b_genPartCharge;   //!
   TBranch        *b_genPartId;   //!
   TBranch        *b_simTrkN;   //!
   TBranch        *b_simTrkPt;   //!
   TBranch        *b_simTrkEta;   //!
   TBranch        *b_simTrkPhi;   //!
   TBranch        *b_simTrkId;   //!
   TBranch        *b_simTrkType;   //!
   TBranch        *b_simTrkVx;   //!
   TBranch        *b_simTrkVy;   //!
   TBranch        *b_simTrkVz;   //!
   TBranch        *b_simVx;   //!
   TBranch        *b_simVy;   //!
   TBranch        *b_simVz;   //!
   TBranch        *b_egCrysN;   //!
   TBranch        *b_egCrysE;   //!
   TBranch        *b_egCrysEt;   //!
   TBranch        *b_egCrysEta;   //!
   TBranch        *b_egCrysPhi;   //!
   TBranch        *b_egCrysGx;   //!
   TBranch        *b_egCrysGy;   //!
   TBranch        *b_egCrysGz;   //!
   TBranch        *b_egCrysCharge;   //!
   TBranch        *b_egN;   //!
   TBranch        *b_egE;   //!
   TBranch        *b_egEt;   //!
   TBranch        *b_egEta;   //!
   TBranch        *b_egPhi;   //!
   TBranch        *b_egGx;   //!
   TBranch        *b_egGy;   //!
   TBranch        *b_egGz;   //!
   TBranch        *b_fHitN;   //!
   TBranch        *b_fHitDisk;   //!
   TBranch        *b_fHitBlade;   //!
   TBranch        *b_fHitSide;   //!
   TBranch        *b_fHitGx;   //!
   TBranch        *b_fHitGy;   //!
   TBranch        *b_fHitGz;   //!
   TBranch        *b_fClSize;   //!
   TBranch        *b_fClSizeX;   //!
   TBranch        *b_fClSizeY;   //!
   TBranch        *b_bHitN;   //!
   TBranch        *b_bHitLayer;   //!
   TBranch        *b_bHitLadder;   //!
   TBranch        *b_bHitGx;   //!
   TBranch        *b_bHitGy;   //!
   TBranch        *b_bHitGz;   //!
   TBranch        *b_bClSize;   //!
   TBranch        *b_bClSizeX;   //!
   TBranch        *b_bClSizeY;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Single_Ele_Eta3_with_Single_Crystal.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Single_Ele_Eta3_with_Single_Crystal.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Single_Ele_Eta3_with_Single_Crystal.root:/NtupleMaker");
      dir->GetObject("t",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pileup = 0;
   genPartE = 0;
   genPartPt = 0;
   genPartEta = 0;
   genPartPhi = 0;
   genPartCharge = 0;
   genPartId = 0;
   simTrkPt = 0;
   simTrkEta = 0;
   simTrkPhi = 0;
   simTrkId = 0;
   simTrkType = 0;
   simTrkVx = 0;
   simTrkVy = 0;
   simTrkVz = 0;
   simVx = 0;
   simVy = 0;
   simVz = 0;
   egCrysE = 0;
   egCrysEt = 0;
   egCrysEta = 0;
   egCrysPhi = 0;
   egCrysGx = 0;
   egCrysGy = 0;
   egCrysGz = 0;
   egCrysCharge = 0;
   egE = 0;
   egEt = 0;
   egEta = 0;
   egPhi = 0;
   egGx = 0;
   egGy = 0;
   egGz = 0;
   fHitDisk = 0;
   fHitBlade = 0;
   fHitSide = 0;
   fHitGx = 0;
   fHitGy = 0;
   fHitGz = 0;
   fClSize = 0;
   fClSizeX = 0;
   fClSizeY = 0;
   bHitLayer = 0;
   bHitLadder = 0;
   bHitGx = 0;
   bHitGy = 0;
   bHitGz = 0;
   bClSize = 0;
   bClSizeX = 0;
   bClSizeY = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bunchN", &bunchN, &b_bunchN);
   fChain->SetBranchAddress("pileup", &pileup, &b_pileup);
   fChain->SetBranchAddress("genPartN", &genPartN, &b_genPartN);
   fChain->SetBranchAddress("genPartE", &genPartE, &b_genPartE);
   fChain->SetBranchAddress("genPartPt", &genPartPt, &b_genPartPt);
   fChain->SetBranchAddress("genPartEta", &genPartEta, &b_genPartEta);
   fChain->SetBranchAddress("genPartPhi", &genPartPhi, &b_genPartPhi);
   fChain->SetBranchAddress("genPartCharge", &genPartCharge, &b_genPartCharge);
   fChain->SetBranchAddress("genPartId", &genPartId, &b_genPartId);
   fChain->SetBranchAddress("simTrkN", &simTrkN, &b_simTrkN);
   fChain->SetBranchAddress("simTrkPt", &simTrkPt, &b_simTrkPt);
   fChain->SetBranchAddress("simTrkEta", &simTrkEta, &b_simTrkEta);
   fChain->SetBranchAddress("simTrkPhi", &simTrkPhi, &b_simTrkPhi);
   fChain->SetBranchAddress("simTrkId", &simTrkId, &b_simTrkId);
   fChain->SetBranchAddress("simTrkType", &simTrkType, &b_simTrkType);
   fChain->SetBranchAddress("simTrkVx", &simTrkVx, &b_simTrkVx);
   fChain->SetBranchAddress("simTrkVy", &simTrkVy, &b_simTrkVy);
   fChain->SetBranchAddress("simTrkVz", &simTrkVz, &b_simTrkVz);
   fChain->SetBranchAddress("simVx", &simVx, &b_simVx);
   fChain->SetBranchAddress("simVy", &simVy, &b_simVy);
   fChain->SetBranchAddress("simVz", &simVz, &b_simVz);
   fChain->SetBranchAddress("egCrysN", &egCrysN, &b_egCrysN);
   fChain->SetBranchAddress("egCrysE", &egCrysE, &b_egCrysE);
   fChain->SetBranchAddress("egCrysEt", &egCrysEt, &b_egCrysEt);
   fChain->SetBranchAddress("egCrysEta", &egCrysEta, &b_egCrysEta);
   fChain->SetBranchAddress("egCrysPhi", &egCrysPhi, &b_egCrysPhi);
   fChain->SetBranchAddress("egCrysGx", &egCrysGx, &b_egCrysGx);
   fChain->SetBranchAddress("egCrysGy", &egCrysGy, &b_egCrysGy);
   fChain->SetBranchAddress("egCrysGz", &egCrysGz, &b_egCrysGz);
   fChain->SetBranchAddress("egCrysCharge", &egCrysCharge, &b_egCrysCharge);
   fChain->SetBranchAddress("egN", &egN, &b_egN);
   fChain->SetBranchAddress("egE", &egE, &b_egE);
   fChain->SetBranchAddress("egEt", &egEt, &b_egEt);
   fChain->SetBranchAddress("egEta", &egEta, &b_egEta);
   fChain->SetBranchAddress("egPhi", &egPhi, &b_egPhi);
   fChain->SetBranchAddress("egGx", &egGx, &b_egGx);
   fChain->SetBranchAddress("egGy", &egGy, &b_egGy);
   fChain->SetBranchAddress("egGz", &egGz, &b_egGz);
   fChain->SetBranchAddress("fHitN", &fHitN, &b_fHitN);
   fChain->SetBranchAddress("fHitDisk", &fHitDisk, &b_fHitDisk);
   fChain->SetBranchAddress("fHitBlade", &fHitBlade, &b_fHitBlade);
   fChain->SetBranchAddress("fHitSide", &fHitSide, &b_fHitSide);
   fChain->SetBranchAddress("fHitGx", &fHitGx, &b_fHitGx);
   fChain->SetBranchAddress("fHitGy", &fHitGy, &b_fHitGy);
   fChain->SetBranchAddress("fHitGz", &fHitGz, &b_fHitGz);
   fChain->SetBranchAddress("fClSize", &fClSize, &b_fClSize);
   fChain->SetBranchAddress("fClSizeX", &fClSizeX, &b_fClSizeX);
   fChain->SetBranchAddress("fClSizeY", &fClSizeY, &b_fClSizeY);
   fChain->SetBranchAddress("bHitN", &bHitN, &b_bHitN);
   fChain->SetBranchAddress("bHitLayer", &bHitLayer, &b_bHitLayer);
   fChain->SetBranchAddress("bHitLadder", &bHitLadder, &b_bHitLadder);
   fChain->SetBranchAddress("bHitGx", &bHitGx, &b_bHitGx);
   fChain->SetBranchAddress("bHitGy", &bHitGy, &b_bHitGy);
   fChain->SetBranchAddress("bHitGz", &bHitGz, &b_bHitGz);
   fChain->SetBranchAddress("bClSize", &bClSize, &b_bClSize);
   fChain->SetBranchAddress("bClSizeX", &bClSizeX, &b_bClSizeX);
   fChain->SetBranchAddress("bClSizeY", &bClSizeY, &b_bClSizeY);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
