// Idea: Load each run, calculate the sum of the max sample in the
// central crystal.  When the central crystal is best centered, it 
// should have the largest sum.  We might need to incorperate the 
// left and right crystal too, let's see if we can make it easy first 

#include <TROOT.h>
#include <TH1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include <iostream>
#include <vector>

const int run_m04_pos = 39;
const int run_m02_pos = 38;
const int run_0_pos = 37;
const int run_p02_pos = 40;
const int run_p04_pos = 41;

const int center_fpga = 1;
const int center_asic = 0;
const int center_connector = 0;

int eeemcal_16i_channel_a_map[16] = { 0,  1,  2,  3,  4,  5,  6,  7,
                                      9, 10, 11, 12, 13, 14, 15, 16};

const int NUM_SAMPLES = 20;

void position_scan() {
    gStyle->SetOptStat(0);
    std::vector<TH1*> position_hists;
    std::vector<int> positions = {-4, -2, 0, 2, 4};
    std::vector<int> runs = {run_m04_pos, run_m02_pos, run_0_pos, run_p02_pos, run_p04_pos};

    TGraphErrors *mean_vs_position = new TGraphErrors(positions.size());

    int i = 0;
    for (auto run : runs) {
        auto path = getenv("OUTPUT_PATH");
        TFile *file = new TFile(Form("%s/run%03d.root", path, run));
        if (!file) {
            std::cerr << "Error opening file" << std::endl;
            return;
        }

        TTree *tree;
        file->GetObject("events", tree);
        if (!tree) {
            std::cerr << "Error getting tree from file" << std::endl;
            return;
        }

        uint waveform[576][NUM_SAMPLES];
        uint tot[576][NUM_SAMPLES];
        tree->SetBranchAddress("adc", &waveform);
        tree->SetBranchAddress("tot", &tot);

        TH1D *adc_sum_hist = new TH1D(Form("adc_sum_hist_run%03d", run), "ADC Sum;ADC;Counts", 500, 0, 8000);
        position_hists.push_back(adc_sum_hist);

        int n_events = tree->GetEntries();
        for (int event = 0; event < n_events; event++) {
            tree->GetEntry(event);
            int adc_sum = 0;
            for (int channel = 0; channel < 16; channel++) {
                double max_adc = 0;
                int actual_channel = eeemcal_16i_channel_a_map[center_connector] + center_fpga*144 + center_asic*72;
                for (int sample = 0; sample < NUM_SAMPLES; sample++) {
                    int sample_adc = waveform[actual_channel][sample] - waveform[actual_channel][0];
                    if (sample_adc > max_adc) {
                        max_adc = sample_adc;
                    }
                }
                adc_sum += max_adc;
            }
            if (adc_sum > 0) {
                adc_sum_hist->Fill(adc_sum);
            }
        }
        TF1 *fit = new TF1("fit", "gaus", 3500, 5000);
        adc_sum_hist->Fit(fit, "QR");
        mean_vs_position->SetPoint(i, positions[i], fit->GetParameter(1));
        mean_vs_position->SetPointError(i, 0, fit->GetParError(1));
        i++;
    }

    TCanvas *c = new TCanvas("c", "c", 1600, 500);
    c->Divide(5, 1);
    for (int i = 0; i < position_hists.size(); i++) {
        c->cd(i+1);
        position_hists[i]->SetTitle(Form("Horizontal Position: %d", positions[i]));
        position_hists[i]->Draw();
        double mean = position_hists[i]->GetFunction("fit")->GetParameter(1);
        double stddev = position_hists[i]->GetFunction("fit")->GetParameter(2);
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.20, 0.85, Form("Mean = %.2f", mean));
        latex.DrawLatex(0.20, 0.80, Form("StdDev = %.2f", stddev));
    }
    c->SaveAs("adc_sum_hist.png");

    for (int i = 0; i < position_hists.size(); i++) {
    }



    
    TF1 *fit = new TF1("fit", "gaus", -4, 4);
    
    mean_vs_position->Fit(fit, "QR");
    
    c = new TCanvas("c2", "c2", 1000, 800);
    mean_vs_position->SetTitle("Mean vs Horizontal Position;Horizontal Position (mm);Mean (ADC)");
    // mean_vs_position->SetMinimum(3500);
    // mean_vs_position->SetMaximum(6000);
    mean_vs_position->SetMarkerStyle(20);
    mean_vs_position->Draw("AP");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    double a = fit->GetParameter(2);
    double ea = fit->GetParError(2);
    double b = fit->GetParameter(1);
    double eb = fit->GetParError(1);
    latex.DrawLatexNDC(0.15, 0.85, Form("Center of fit: %.03f#pm %.03f mm", fit->GetParameter(1), fit->GetParError(1)));
    c->SaveAs("mean_vs_position.png");
}