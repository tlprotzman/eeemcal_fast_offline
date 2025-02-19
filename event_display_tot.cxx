#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TError.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLatex.h>

#include <iostream>
#include <vector>

const int NUM_SAMPLES = 20;

// EEEMCal mapping - instead of "layers", we have a single plane, where each crystal is one connector
// FPGA IP | ID
// 208     | 0
// 209     | 1
// 210     | 2
// 211     | 3
int eeemcal_fpga_map[25] = {0, 3, 3, 0, 3,
                            2, 1, 1, 1, 2,
                            2, 1, 1, 1, 3,
                            2, 2, 1, 2, 3,
                            2, 0, 0, 1, 2};

// ASIC | ID
// 0    | 0
// 1    | 1
int eeemcal_asic_map[25] = { 1, 1, 1, 0, 0,
                             1, 1, 1, 1, 1,
                             1, 0, 0, 0, 0,
                             1, 0, 1, 0, 0,
                             0, 1, 1, 0, 0};

// Connector | ID
// A        | 0
// B        | 1
// C        | 2
// D        | 3
int eeemcal_connector_map[25] = { 2,  0,  1,  0,  1,
                                  0,  2,  0,  3,  3,
                                  1,  2,  0,  3,  0,
                                  2,  0,  1,  1,  2,
                                  3,  1,  1,  1,  2};

int eeemcal_16i_channel_a_map[16] = { 0,  1,  2,  3,  4,  5,  6,  7,
              9, 10, 11, 12, 13, 14, 15, 16};

int eeemcal_16i_channel_b_map[16] = {19, 20, 21, 22, 23, 24, 25, 26,
             27, 28, 29, 30, 31, 32, 33, 34};

int eeemcal_16i_channel_c_map[16] = {55, 56, 57, 58, 59, 60, 61, 62,
             63, 64, 65, 66, 67, 68, 69, 70};

int eeemcal_16i_channel_d_map[16] = {36, 37, 38, 39, 40, 41, 42, 43,
             45, 46, 47, 48, 49, 50, 51, 52};

int *eeemcal_16i_channel_map[4] = {eeemcal_16i_channel_a_map, eeemcal_16i_channel_b_map, eeemcal_16i_channel_c_map, eeemcal_16i_channel_d_map};

int eeemcal_4x4_channel_a_map[4] = {0, 4, 9, 12};
int eeemcal_4x4_channel_b_map[4] = {20, 24, 27, 31};
int eeemcal_4x4_channel_c_map[4] = {58, 62, 65, 69};
int eeemcal_4x4_channel_d_map[4] = {38, 42, 48, 52};
int *eeemcal_4x4_channel_map[4] = {eeemcal_4x4_channel_a_map, eeemcal_4x4_channel_b_map, eeemcal_4x4_channel_c_map, eeemcal_4x4_channel_d_map};

int eeemcal_16p_channel_map[4] = {6, 26, 63, 46};

// Mode | Readout
// 0    | 16i
// 1    | 4x4
// 2    | 16p
void event_display_tot(int run, int event=0, int mode=0) {
    gErrorIgnoreLevel = kWarning;
    // Read in the waveforms
    auto path = getenv("OUTPUT_PATH");
    TFile *file = new TFile(Form("%s/run%03d.root", path, run));
    if (!file) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Get the tree from the file
    TTree *tree;
    file->GetObject("events", tree);
    if (!tree) {
        std::cerr << "Error getting tree from file" << std::endl;
        return;
    }

    // Get number of events
    int n_events = tree->GetEntries();
    if (event == -1 || event >= n_events) {
        std::cerr << "Run has " << n_events << " events" << std::endl;
        return;
    }
    
    TCanvas *c = new TCanvas("c", "c", 1600, 900);
    bool open = false;
    // Set the branch addresses
    uint waveform[576][NUM_SAMPLES];
    tree->SetBranchAddress("tot", &waveform);
    for (event = 0; event < 10; event++) {
        std::cout << "\rEvent " << event << std::flush;
        c = new TCanvas(Form("c_%d", event), "c", 1600, 1200);
        tree->GetEntry(event);
        
        gStyle->SetOptStat(0);

        // Draw the event
        std::vector<TGraph*> graphs;
        c->cd(0);
        auto label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.05);
        label->DrawLatex(0.05, 0.9, Form("TOT: Run %d, Event %d", run, event));
        auto pad = new TPad("pad", "pad", 0.05, 0.05, 0.95, 0.85);
        pad->Draw();
        // Add text to the top of the pad with the run and event number
        pad->cd();
        pad->Divide(5, 5, 0.001, 0.001);
        for (int crystal = 0; crystal < 25; crystal++) {
            if (mode == 0) {    // 16i
                pad->cd(crystal+1);
                gPad->Divide(4, 4, 0, 0);
                for (int sipm = 0; sipm < 16; sipm++) {
                    pad->cd(crystal+1);
                    gPad->cd(sipm+1);
                    int channel_number = eeemcal_fpga_map[crystal]*144 + eeemcal_asic_map[crystal]*72 +  eeemcal_16i_channel_map[eeemcal_connector_map[crystal]][sipm];
                    TGraph *g = new TGraph(NUM_SAMPLES-1);
                    g->SetTitle(Form("crystal_%d_sipm_%d_event_%d_ch_%d", crystal, sipm, event, channel_number));
                    g->GetXaxis()->SetTitle("Sample");
                    g->GetYaxis()->SetTitle("TOT Counts");
                    g->GetXaxis()->SetRange(0, NUM_SAMPLES-1);
                    // g->SetLineColor(sipm+1); // Different color for each SiPM
                    // g->SetMarkerColor(sipm+1);
                    g->SetMarkerStyle(20);
                    g->SetMarkerSize(0.5);
                    graphs.push_back(g);
                    double max_signal = 0;
                    for (int sample = 0; sample < NUM_SAMPLES; sample++) {
                        double signal = waveform[channel_number][sample];
                        g->SetPoint(sample, sample + 0.5, signal);
                        if (signal > max_signal) {
                            max_signal = signal;
                        }
                    }
                    g->SetMinimum(0);
                    g->SetMaximum(4096);
                    g->Draw("APL");
                    // Set background color based on max signal
                    int red = (int)(255 * max_signal / 4096);
                    int green = (int)(255 * (1 - abs(max_signal - 2048) / 2048));
                    int blue = (int)(255 * (1 - max_signal / 4096));
                    gPad->SetFillColorAlpha(TColor::GetColor(red, green, blue), 0.2);
                }
                // Draw borders between 4x4 groups
                if (crystal % 5 == 4) {
                    gPad->SetFrameLineWidth(2);
                    gPad->SetFrameLineColor(kBlack);
                }
            }
        }
        c->Draw();
        if (open) {
            c->SaveAs(Form("output/Run%03d_event_display_tot.pdf", run));
        } else {
            c->SaveAs(Form("output/Run%03d_event_display_tot.pdf(", run));
            open = true;
        }
    }
    std::cout << std::endl;
    c->SaveAs(Form("output/Run%03d_event_display_tot.pdf)", run));
}