#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TError.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLatex.h>
#include <TF1.h>

#include <iostream>
#include <vector>

const int NUM_SAMPLES = 20;

void adc_tot_correlation(int run) {
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);
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

    // Set the branch addresses
    uint waveform[576][NUM_SAMPLES];
    uint tot[576][NUM_SAMPLES];
    tree->SetBranchAddress("adc", &waveform);
    tree->SetBranchAddress("tot", &tot);

    TCanvas *c = new TCanvas("c", "c", 1600, 900);
    std::vector<TH2F*> hists;
    for (int channel = 0; channel < 576; channel++) {
        TH2F *hist = new TH2F(Form("adc_tot_ch%d", channel), "ADC vs TOT;Max ADC;Max TOT", 1024/8, 0, 1024, 4096/32, 0, 4096);
        hists.push_back(hist);
    }

    int n_events = tree->GetEntries();
    for (int event = 0; event < n_events; event++) {
        tree->GetEntry(event);

        for (int channel = 0; channel < 576; channel++) {
            // std::cerr << "\rEvent " << event << ", Channel " << channel << std::flush;
            double max_adc = 0;
            double max_tot = 0;
            for (int sample = 0; sample < NUM_SAMPLES; sample++) {
                if (waveform[channel][sample] > max_adc) {
                    max_adc = waveform[channel][sample];
                }
                if (tot[channel][sample] > max_tot) {
                    max_tot = tot[channel][sample];
                }
            }
            if (max_adc > 150 && max_adc < 950) {
                hists[channel]->Fill(max_adc, max_tot);
            }
        }
    }

    c->cd();
    bool open = false;
    TF1 *fit = new TF1("fit", "[0]*x+[1]", 150, 950);
    for (int channel = 0; channel < 576; channel++) {
        hists[channel]->Fit(fit, "QR");
        hists[channel]->Draw("COLZ");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.12, 0.85, Form("Channel %d", channel));
        latex.DrawLatex(0.12, 0.81, Form("Slope: %.2f", fit->GetParameter(0)));
        latex.DrawLatex(0.12, 0.77, Form("Fit #chi^{2}/NDF: %.2f", fit->GetChisquare()/fit->GetNDF()));
        if (open) {
            c->SaveAs(Form("output/adc_tot_correlation_run%03d.pdf", run));
        } else {
            c->SaveAs(Form("output/adc_tot_correlation_run%03d.pdf(", run));
            open = true;
        }
    }
    c->SaveAs(Form("output/adc_tot_correlation_run%03d.pdf)", run));
    std::cout << "done" << std::endl;
}
