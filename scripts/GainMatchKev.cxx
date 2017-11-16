// To compile:
// Note: GRSISort looks for .cxx extensions when compiling (for example it looks
// in the myAnalysis directory)
// Alternatively you may use the following to compile:
// g++ GainMatch.C -o MyGainMatch -std=c++0x -I$GRSISYS/GRSISort/include/
// `grsi-config --cflags --all-libs --root`

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <math.h> // round, floor, ceil, trunc
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
//#include "TFileIter.h"
#include "TApplication.h"
#include "TF1.h"
#include "TGainMatch.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TPPG.h"
#include "TPeak.h"
#include "TROOT.h"
#include "TScaler.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include "TTree.h"

#ifndef __CINT__
#include "TGriffin.h"
#include "TSceptar.h"
#endif

///// GLOBAL VARIABLES /////

// When using fragment trees, set ISCALIBRATION = 1
// For analysis trees, set ISCALIBRATION = 0
const bool IS_CALIBRATION = 0;
// Set 1 for Quadradic peak fit centroids
// Set 0 for linear fit centroids
const bool IS_QUADRADIC = 0;

// Input the peaks that are needed for gain matching
/*
const std::vector<double_t> CALIBRATION_PEAKS = {121.7817, 244.6974, 964.057,
                                                 1085.837, 1112.076, 1408.013};

const std::vector<double_t> WIDTH_PEAKS = {16, 20, 20, 20, 20, 30};
*/

const std::vector<double_t> CALIBRATION_PEAKS = 
    {138.2,
    208.52,
    445.98,
    560.21,
    683.06,
    813.06,
    813.20,
    1050,65,
    1173.52,
    1229.64,
    1504.10,
    2042.77,
    2327.82,
    2475.06,
    2677.4 };

const std::vector<double_t> WIDTH_PEAKS = 
    {20, 
    20,  
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    20,
    20};

//// END GLOBAL VARIABLES ////

void CalibrateData(TTree *intree);

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage is: %s <fragment or analysis tree file> ).\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Setup the variables and object that will be used
    TFile *infile = new TFile(argv[1], "UPDATE");

    // We check if the file is open two ways. This might be exessive.
    if (infile == nullptr) {
        printf("Failed to open file '%s'!\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    if (!infile->IsOpen()) {
        printf("Failed to open file '%s'!\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    TTree *intree = nullptr;
    if (IS_CALIBRATION)
        intree = (TTree *)infile->Get("FragmentTree");
    else
        intree = (TTree *)infile->Get("AnalysisTree");

    if (!intree) {
        printf("Failed to find correct TTree,"
               "Check if IS_CALIBRATION is set for ether"
               " fragment or analysis tree.\n");
        exit(EXIT_FAILURE);
    }

    CalibrateData(intree);

    delete infile;
    exit(EXIT_SUCCESS);
}

void CalibrateData(TTree *intree) {
    // Check if global variables for peak fitting are sane
    if (CALIBRATION_PEAKS.size() != WIDTH_PEAKS.size()) {
        printf("Global variables that describe peaks are not sane\n"
               "Check CALIBRATION_PEAKS or WIDTH_PEAKS\n");
        std::cout << "Note: CALIBRATION_PEAKS is " << CALIBRATION_PEAKS.size() << " elements.\n" ;
        std::cout << "and WIDTH_PEAKS is " <<  WIDTH_PEAKS.size() << " elements.\n";
        exit(EXIT_FAILURE);
    }

    int CaliPeakQnty = CALIBRATION_PEAKS.size();

    std::vector<double_t> PeakCentroids;

    TChannel *channel = nullptr;
    TChannel::ReadCalFromTree(intree);

    // Keep statistics at 1 keV/bin
    TH2D *mat_en = new TH2D("mat_en", "", 64, 0, 64, 5000, 0, 5000);

    if (IS_CALIBRATION)
        intree->Project("mat_en",
                        "TFragment.GetEnergy():TFragment.GetChannelNumber()");
    else
        // intree->Project("mat_en","TGriffin.fGriffinLowGainHits.GetEnergy():TGriffin.fGriffinHits.GetChannel().fNumber","TSceptar.GetMultiplicity()>0");
        // for tin data
        intree->Project("mat_en", "TGriffin.fGriffinLowGainHits.GetEnergy():"
                                  "TGriffin.fGriffinLowGainHits.GetChannel()."
                                  "fNumber"); // for cal files

    for (int i = 0; i < 64; i++) { // 64 individual clovers in the whole array

        // Load the current channel/calibration
        channel = TChannel::GetChannelByNumber(i);

        // Project a matrix for each channel
        TH1D *h_en = mat_en->ProjectionY(Form("h_%.2i", i), i + 1, i + 1);

        for (int k = 0; k < CaliPeakQnty; k++) {
            TSpectrum s;
            h_en->GetXaxis()->SetRangeUser(
                CALIBRATION_PEAKS[k] - WIDTH_PEAKS[k],
                CALIBRATION_PEAKS[k] + WIDTH_PEAKS[k]);
            s.Search(h_en, 2, "", 0.25);
            double_t DataPeak = s.GetPositionX()[0];
            h_en->GetXaxis()->UnZoom();

            TPeak *TempPeak = new TPeak(DataPeak, DataPeak - WIDTH_PEAKS[k],
                                        DataPeak + WIDTH_PEAKS[k]);
            TempPeak->Fit(h_en, "MQ+");
            PeakCentroids.push_back(TempPeak->GetCentroid());
            delete TempPeak;
        }

        // Output status of peak finding
        /*
        printf("\nFor channel %d, centroids found at ", i);
        for (int k = 0; k < CaliPeakQnty - 1; k++)
                        printf("%.2f, ", PeakCentroids[k]);
        if( PeakCentroids.size() == 1 ) // Only one centroid
                printf( "%.2f\n", PeakCentroids.back() );
        else
                printf( "and %.2f\n", PeakCentroids.back() );
        */

        TF1 *fit = nullptr;
        TGraph *graph = nullptr;
        TFitResultPtr fitptr;

        // new coefficients are convolutions of existing coefficients with new
        // fit parameters
        double a = 0; // offset
        double b = 0; // linear coefficient
        double c = 0; // quadratic coefficient

        if (!IS_QUADRADIC) // linear fit
        {
            fit = new TF1("fit", "pol1", 0, 5000);
            graph = new TGraph(CaliPeakQnty, CALIBRATION_PEAKS.data(),
                               PeakCentroids.data());
            fitptr = graph->Fit(fit, "IEMNCFQS");

            double slope = fitptr->Parameter(1);
            double offset = fitptr->Parameter(0);
            double old_offset = channel->GetENGCoeff()[0];
            double old_slope = channel->GetENGCoeff()[1];
            a = (old_offset - offset)/slope;
            b = (old_slope/slope);
            printf("slope: %.2f, offset %.2f\n", slope, offset);
            /*
            a = fitptr->Parameter(0) +
                fitptr->Parameter(1) * channel->GetENGCoeff()[0];
            b = fitptr->Parameter(1) * channel->GetENGCoeff()[1];
            */
        } else // quadratic fit
        {
            fit = new TF1("fit", "pol2", 0, 5000);
            graph = new TGraph(CaliPeakQnty, PeakCentroids.data(),
                               CALIBRATION_PEAKS.data());
            fitptr = graph->Fit(fit, "QS");

            a = fitptr->Parameter(0) +
                fitptr->Parameter(1) * channel->GetENGCoeff()[0] +
                fitptr->Parameter(2) * pow(channel->GetENGCoeff()[0], 2);
            b = fitptr->Parameter(1) * channel->GetENGCoeff()[1] +
                2 * fitptr->Parameter(2) * channel->GetENGCoeff()[0] *
                    channel->GetENGCoeff()[1];
            c = fitptr->Parameter(2) * pow(channel->GetENGCoeff()[1], 2);
        }

        if (channel->GetENGCoeff().size() == 2) // existing linear calibration
            printf("Old coefficients are\t%f\t%f\n", channel->GetENGCoeff()[0],
                   channel->GetENGCoeff()[1]);
        else // existing quadratic calibration
            printf("Old coefficients are\t%f\t%f\t%f\n",
                   channel->GetENGCoeff()[0], channel->GetENGCoeff()[1],
                   channel->GetENGCoeff()[2]);

        if (!IS_QUADRADIC) // linear fit
            printf("New coefficients are\t%f\t%f\n", a, b);
        else // quadratic fit
            printf("New coefficients are\t%f\t%f\t%.7f\n", a, b, c);

        // save new coefficients
        channel->DestroyENGCal();
        channel->AddENGCoefficient((Float_t)a);     // offset
        channel->AddENGCoefficient((Float_t)b);     // slope
        if (IS_QUADRADIC)                           // quadratic fit
            channel->AddENGCoefficient((Float_t)c); // quadratic coefficient

        // Cleanup from loop
        delete fit;
        delete graph;
    }

    // If using individual files, save the coefficients directly in the analysis
    // file, and save time later
    // For grouped files, there are multiple TChannels, causing errors;
    // the cal file will have to be read in individually to the analysis files,
    // and each analysis file re-written with the new coefficients
    TChannel::WriteToRoot();
    //TChannel::WriteCalFile("./test.cal");

    delete mat_en;
    printf("\n");

    return;
}
