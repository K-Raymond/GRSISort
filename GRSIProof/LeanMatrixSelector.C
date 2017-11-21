#define LeanMatrixSelector_cxx
// The class definition in LeanMatrixSelector.h has been generated automatically
#include "LeanMatrixSelector.h"

void LeanMatrixSelector::CreateHistograms() {

    // Define histogram limits

    double low = 0;
    double high = 8192;
    double nofBins = 8192;

    // Define coincidence parameters

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

    // Define Histograms

    // spectra
    // singles
    fH1["gE"] = new TH1D("gE", "#gamma Singles", nofBins, low, high);
    fH1["gE_b"] = new TH1D("gE_b", "#gamma Singles in rough #beta coincidence",
                           nofBins, low, high);
    /*  fH1["gE_b-prompt-all"] = new TH1D("gE_b-prompt", "#gamma Singles in
    prompt
    #beta coincidence", nofBins, low, high);
    fH1["gE_b-prompt-implant"] = new TH1D("gE_b-prompt-implant","#gamma Singles
    in
    prompt #beta coincidence - implant window", nofBins, low, high);
    fH1["gE_b-prompt-decay"] = new TH1D("gE_b-prompt-decay","#gamma Singles in
    prompt #beta coincidence - decay window", nofBins, low, high);
    fH1["gE_b-tbg-implant"] = new TH1D("gE_b-tbg-implant","#gamma Singles in
    time
    random  #beta coincidence correction  - implant window", nofBins, low,
    high);
    fH1["gE_b-tbg-decay"] = new TH1D("gE_b-tbg-decay","#gamma Singles in prompt
    #beta coincidence - decay window", nofBins, low, high);
    */
    // what are the tcorr matrices good for??
    // What are the timing spectra good for??

    // addback
    fH1["aE"] = new TH1D("aE", "Addback Singles", nofBins, low, high);

    // matrices
    fH2["ggE"] = new TH2F("ggE", "#gamma #gamma Coincidence", 6000, 0, 3000,
                          6000, 0, 3000);
    fH2["gg_prompt"] =
        new TH2F("gg_prompt", "#gamma #gamma Coincidence - prompt window",
                 nofBins, low, high, nofBins, low, high);
    // Send histograms to Output list to be added and written.
    for (auto it : fH1)
        GetOutputList()->Add(it.second);
    for (auto it : fH2)
        GetOutputList()->Add(it.second);
    for (auto it : fHSparse)
        GetOutputList()->Add(it.second);
}

bool PromptCoincidence(TGriffinHit *g, TSceptarHit *s) {
    // Check if hits are less then 300 ns apart.
    return std::fabs(g->GetTime() - s->GetTime()) < 300.;
}
/*
bool PromtCoincidenceImplant(TGriffinHit *g, TSceptarHit *s){
  return std::fabs((g->GetTime() - s->GetTime())) < 300) && ()

}

bool ImplantWindow(TGriffinHit *g){
  return std::fabs(g->GetTime() => CycleBeamOnStart && g->GetTime() =<
CycleBeamOffStart);

}
*/

bool PromptCoincidence(TGriffinHit *g1, TGriffinHit *g2) {
    // Check if hits are less then 500 ns apart.
    return std::fabs(g1->GetTime() - g2->GetTime()) < 500.;
}

void LeanMatrixSelector::FillHistograms() {

    // define time within cycle, since PPG was not being saved

    // Loop over all Griffin Hits
    for (int i = 0; i < fGrif->GetMultiplicity(); ++i) {
        //  gTimeInCycle = (long long)((long)g->GetHit(i)->GetTime())%CycleEnd);
        fH1.at("gE")->Fill(fGrif->GetGriffinHit(i)->GetEnergy());

        for (int j = 0; j < fGrif->GetMultiplicity(); ++j) {
            if (i == j)
                continue;
            if (PromptCoincidence(fGrif->GetGriffinHit(i),
                                  fGrif->GetGriffinHit(j))) {
                if (ggTlow <= std::fabs(fGrif->GetGriffinHit(j)->GetTime() -
                                        fGrif->GetGriffinHit(i)->GetTime()) &&
                    std::fabs(fGrif->GetGriffinHit(j)->GetTime() -
                              fGrif->GetGriffinHit(i)->GetTime()) < ggThigh) {
                    fH2.at("gg_prompt")
                        ->Fill(fGrif->GetGriffinHit(i)->GetEnergy(),
                               fGrif->GetGriffinHit(j)->GetEnergy());
                }
            }

            //  fH2.at("ggE")->Fill(fGrif->GetGriffinHit(i)->GetEnergy(),
            //  fGrif->GetGriffinhit()
        }
        // Loop over all sceptar hits
        for (auto j = 0; j < fScep->GetMultiplicity(); ++j) {

            //	  bTimeInCycle = (long
            // long)((long)s->GetHit(j)->GetTime())%CycleEnd);
            if (PromptCoincidence(fGrif->GetGriffinHit(i),
                                  fScep->GetSceptarHit(j))) {
                fH1.at("gE_b")->Fill(fGrif->GetGriffinHit(i)->GetEnergy());
            }
        }
    }
}
