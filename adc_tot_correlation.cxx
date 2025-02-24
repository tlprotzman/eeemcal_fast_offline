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
#include <TLine.h>

#include <iostream>
#include <vector>
#include <ostream>

const int NUM_SAMPLES = 20;

int lfhcal_channel_map[72] = {64, 63, 66, 65, 69, 70, 67, 68,
    55, 56, 57, 58, 62, 61, 60, 59,
    45, 46, 47, 48, 52, 51, 50, 49,
    37, 36, 39, 38, 42, 43, 40, 41,
    34, 33, 32, 31, 27, 28, 29, 30,
    25, 26, 23, 24, 20, 19, 22, 21,
    16, 14, 15, 12,  9, 11, 10, 13,
     7,  6,  5,  4,  0,  1,  2,  3,
     -1, -1, -1, -1, -1, -1, -1, -1};

// EEEMCal mapping - instead of "layers", we have a single plane, where each crystal is one connector
// FPGA IP | ID
// 208     | 0
// 209     | 1
// 210     | 2
// 211     | 3
int eeemcal_fpga_map[25] = {0, 3, 3, 0, 3,
         2, 1, 1, 1, 2,
         2, 1, 1, 1, 3,
         2, 2, 1, 1, 3,
         2, 0, 0, 1, 2};

// ASIC | ID
// 0    | 0
// 1    | 1
int eeemcal_asic_map[25] = { 1, 1, 1, 0, 0,
          1, 1, 1, 1, 1,
          1, 0, 0, 0, 0,
          1, 0, 1, 1, 0,
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
               3,  1,  3,  1,  2};

int eeemcal_16i_channel_a_map[16] = { 2,  6, 11, 15,  0,  4,  9, 13,
    1,  5, 10, 14,  3,  7, 12, 16};

int eeemcal_16i_channel_b_map[16] = {20, 24, 29, 33, 18, 22, 27, 31,
   19, 23, 28, 32, 21, 25, 30, 34};

int eeemcal_16i_channel_c_map[16] = {67, 63, 59, 55, 69, 65, 61, 57,
   70, 66, 60, 56, 68, 64, 58, 54};
     
int eeemcal_16i_channel_d_map[16] = {50, 46, 40, 36, 52, 48, 42, 38,
   51, 47, 43, 39, 49, 45, 41, 37};

int *eeemcal_16i_channel_map[4] = {eeemcal_16i_channel_a_map, eeemcal_16i_channel_b_map, eeemcal_16i_channel_c_map, eeemcal_16i_channel_d_map};

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
            int tot_val = 0;
            int tot_sample = 0;
            int adc_val = 0;
            for (int sample = 0; sample < NUM_SAMPLES; sample++) {
                if (tot[channel][sample] > tot_val) {
                    tot_val = tot[channel][sample];
                    tot_sample = sample;
                }
                if (waveform[channel][sample] > adc_val) {
                    adc_val = waveform[channel][sample];
                }
            }
            if (tot_val > 5) {
                adc_val -= waveform[channel][0];
                if (adc_val > 200 && waveform[channel][tot_sample] < 1000) {
                    hists[channel]->Fill(adc_val, tot_val);
                }
            }
        }
    }
    
    bool open = false;
    std::vector<float> slopes(576);
    std::vector<float> slope_errors(576);
    std::vector<float> intercepts(576);
    std::vector<float> intercept_errors(576);
    for (int i = 0; i < 576; i++) {
        slopes[i] = 0;
        slope_errors[i] = 0;
    }
    for (int crystal = 0; crystal < 25; crystal++) {
        TCanvas *c = new TCanvas("c", "c", 1600, 900);
        c->cd();
        auto label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.05);
        label->DrawLatex(0.05, 0.9, Form("Crystal %d ADC TOT Correlation", crystal));
    
    
        auto pad = new TPad("pad", "pad", 0.05, 0.05, 0.95, 0.85);
        pad->Draw();
        // Add text to the top of the pad with the run and event number
        pad->cd();
        pad->Divide(4, 4, 0.000, 0.000);

        for (int sipm = 0; sipm < 16; sipm++) {
            pad->cd(sipm+1);
            int fit_start = 700;
            int fit_end = 900;
            TF1 *fit = new TF1("fit", "[0]*x+[1]", fit_start, fit_end);
            
            int channel_fpga = eeemcal_fpga_map[crystal];
            int channel_asic = eeemcal_asic_map[crystal];
            int channel_connector = eeemcal_connector_map[crystal];
            int channel = 144 * channel_fpga + 72 * channel_asic + eeemcal_16i_channel_map[channel_connector][sipm];
            
            hists[channel]->Fit(fit, "QR");
            slopes[144 * eeemcal_fpga_map[crystal] + 72 * eeemcal_asic_map[crystal] + eeemcal_16i_channel_map[eeemcal_connector_map[crystal]][sipm]] = fit->GetParameter(0);
            slope_errors[144 * eeemcal_fpga_map[crystal] + 72 * eeemcal_asic_map[crystal] + eeemcal_16i_channel_map[eeemcal_connector_map[crystal]][sipm]] = fit->GetParError(0);
            intercepts[144 * eeemcal_fpga_map[crystal] + 72 * eeemcal_asic_map[crystal] + eeemcal_16i_channel_map[eeemcal_connector_map[crystal]][sipm]] = fit->GetParameter(1);
            intercept_errors[144 * eeemcal_fpga_map[crystal] + 72 * eeemcal_asic_map[crystal] + eeemcal_16i_channel_map[eeemcal_connector_map[crystal]][sipm]] = fit->GetParError(1);
            hists[channel]->Draw("COLZ");
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.03);
            latex.DrawLatex(0.12, 0.85, Form("Channel %d", channel));
            latex.DrawLatex(0.12, 0.81, Form("Slope: %.2f#pm%.4f", fit->GetParameter(0), fit->GetParError(0)));
            latex.DrawLatex(0.12, 0.77, Form("Intercept: %.2f#pm%.4f", fit->GetParameter(1), fit->GetParError(1)));
            TLine *line1 = new TLine(fit_start, 0, fit_start, hists[channel]->GetYaxis()->GetXmax());
            line1->SetLineColor(kRed);
            line1->SetLineStyle(2);
            line1->Draw();
    
            TLine *line2 = new TLine(fit_end, 0, fit_end, hists[channel]->GetYaxis()->GetXmax());
            line2->SetLineColor(kRed);
            line2->SetLineStyle(2);
            line2->Draw();
    
    
        }
        // latex.DrawLatex(0.12, 0.77, Form("Fit #chi^{2}/NDF: %.2f", fit->GetChisquare()/fit->GetNDF()));
        if (open) {
            c->SaveAs(Form("output/Run%03d_adc_tot_correlation.pdf", run));
        } else {
            c->SaveAs(Form("output/Run%03d_adc_tot_correlation.pdf(", run));
            open = true;
        }
    }
    TCanvas *c = new TCanvas("c", "c", 1600, 900);
    c->SaveAs(Form("output/Run%03d_adc_tot_correlation.pdf)", run));
    std::cout << "done" << std::endl;

    // Write slopes and intercepts to root file
    auto slopes_histogram = new TH1F("adc_tot_slope", "Slopes;Channel;Slope", 576, 0, 576);
    auto intercepts_histogram = new TH1F("adc_tot_intercept", "Intercepts;Channel;Intercept", 576, 0, 576);
    for (int i = 0; i < 576; i++) {
        slopes_histogram->SetBinContent(i, slopes[i]);
        slopes_histogram->SetBinError(i, slope_errors[i]);
        intercepts_histogram->SetBinContent(i, intercepts[i]);
        intercepts_histogram->SetBinError(i, intercept_errors[i]);
    }
    TFile *output_file = new TFile(Form("output/Run%03d_adc_tot_correlation.root.new", run), "RECREATE");
    slopes_histogram->Write();
    intercepts_histogram->Write();
    output_file->Close();
}
