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
#include <string>

const int center_fpga = 1;
const int center_asic = 0;
const int center_connector = 0;

int eeemcal_16i_channel_a_map[16] = { 0,  1,  2,  3,  4,  5,  6,  7,
                                      9, 10, 11, 12, 13, 14, 15, 16};

const int NUM_SAMPLES = 20;

void position_scan() {
    int mode = 0;   // 0 horizontal, 1 vertical
    gStyle->SetOptStat(0);
    std::vector<TH1*> position_hists;
    std::vector<int> h_positions = {-4, -2, 0, 2, 4};
    std::vector<int> h_runs = {39, 38, 37, 40, 41};

    std::vector<int> v_positions = {-4, -2, 0, 2, 4, 6, 8};
    std::vector<int> v_runs = {52, 51, 47, 48, 50, 54, 55};

    std::vector<int> positions;
    std::vector<int> runs;

    std::string axis;

    if (mode == 0) {
        positions = h_positions;
        runs = h_runs;
        axis = "Horizontal";
    } else {
        positions = v_positions;
        runs = v_runs;
        axis = "Vertical";
    }

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
        // TF1 *fit = new TF1("fit", "gaus", 200, 325);
        adc_sum_hist->Fit(fit, "QR");
        mean_vs_position->SetPoint(i, positions[i], fit->GetParameter(1));
        mean_vs_position->SetPointError(i, 0, fit->GetParError(1));
        i++;
    }

    TCanvas *c = new TCanvas("c", "c", 1600, 500);
    c->Divide(position_hists.size(), 1);
    for (int i = 0; i < position_hists.size(); i++) {
        c->cd(i+1);
        position_hists[i]->SetTitle(Form("%s Position: %d", axis.c_str(), positions[i]));
        position_hists[i]->Draw();
        double mean = position_hists[i]->GetFunction("fit")->GetParameter(1);
        double stddev = position_hists[i]->GetFunction("fit")->GetParameter(2);
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.20, 0.85, Form("Mean = %.2f", mean));
        latex.DrawLatex(0.20, 0.80, Form("StdDev = %.2f", stddev));
        latex.DrawLatex(0.20, 0.75, Form("StdDev/Mean = %.2f", stddev/mean));
    }
    c->SaveAs("adc_sum_hist.png");
    // return;

    for (int i = 0; i < position_hists.size(); i++) {
    }


    float min = *std::min_element(positions.begin(), positions.end());
    float max = *std::max_element(positions.begin(), positions.end());
    
    TF1 *fit = new TF1("fit", "gaus", min, max);
    
    mean_vs_position->Fit(fit, "QR");
    
    c = new TCanvas("c2", "c2", 1000, 800);
    mean_vs_position->SetTitle(Form("Mean vs %s Position;%s Position (mm);Mean (ADC)", axis.c_str(), axis.c_str()));
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
    if (mode == 0) {
        c->SaveAs("horizontal_position.png");
    } else if (mode == 1) {
        c->SaveAs("vertical_position.png");
    }
}