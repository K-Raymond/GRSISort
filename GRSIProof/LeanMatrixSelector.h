//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 25 13:18:27 2016 by ROOT version 5.34/24
// from TTree FragmentTree/FragmentTree
// found on file: fragment07844_000.root
//////////////////////////////////////////////////////////

#ifndef LeanMatrixSelector_h
#define LeanMatrixSelector_h

#include "TChain.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"

// Header file for the classes stored in the TTree if any.
#include "TGRSISelector.h"
#include "TGriffin.h"
#include "TSceptar.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class LeanMatrixSelector : public TGRSISelector {
  private:
    double ggTlow = 0.;
    double ggThigh = 316.;
    double gbTlow = -75.;
    double gbThigh = 350.;
    double ggTunknownlow = 500;
    double ggTunknownhigh = 900;

    double ggBGlow = 1000.;
    double ggBGhigh = 2000.;
    double gbBGlow = 600.;
    double gbBGhigh = 2000.;
    double ggBGScale = (ggThigh - ggTlow) / (ggBGhigh - ggBGlow);
    double gbBGScale = (gbThigh - gbTlow) / (gbBGhigh - gbBGlow);

    // Define duty cycle

    double CycleTapeStart = 0;
    double CycleBgStart = 0;
    double CycleBeamOnStart = 0;
    double CycleBeamOffStart = 0;
    long CycleEnd = 0;

    // Support for residuals
    std::vector<TGraph*> gEngResidualVec;

  public:
    TGriffin *fGrif; // Pointers to spot that events will be
    TSceptar *fScep;

    LeanMatrixSelector(TTree * /*tree*/ = 0)
        : TGRSISelector(), fGrif(0), fScep(0) {
        SetOutputPrefix("ExampleEvent"); // Changes prefix of output file
    }
    // These functions are expected to exist
    virtual ~LeanMatrixSelector() {}
    virtual Int_t Version() const { return 2; }
    void CreateHistograms();
    void FillHistograms();
    void InitializeBranches(TTree *tree);

    ClassDef(LeanMatrixSelector, 2); // Makes ROOT happier
};

#endif

#ifdef LeanMatrixSelector_cxx
void LeanMatrixSelector::InitializeBranches(TTree *tree) {
    
    if (!tree)
        return;
    if (tree->SetBranchAddress("TGriffin", &fGrif) == TTree::kMissingBranch) {
        fGrif = new TGriffin;
    }
    if (tree->SetBranchAddress("TSceptar", &fScep) == TTree::kMissingBranch) {
        fScep = new TSceptar;
    }
       // Add in residuals
    printf("Attemping to load residuals...");
    if( gFile->cd("Energy_Residuals") ) {
        printf(" found, loading\n");
        TGraph* TempResidual;
        for( int i = 0; i < 64 ; i++ ) {
            gDirectory->GetObject(Form("Graph;%d", i), TempResidual );
            //fGrif->LoadEnergyResidual(i,TempResidual);
        }
        gFile->cd();
        printf("Residuals loaded!\n");
    } else {
        printf(" no residuals found.\n");
    }
}

#endif // #ifdef LeanMatrixSelector_cxx
