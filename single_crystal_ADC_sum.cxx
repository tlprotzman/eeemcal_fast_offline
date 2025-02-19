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

int eeemcal_16i_channel_a_map[16] = { 2,  6, 11, 15,  0,  4,  9, 13,
                                      1,  5, 10, 14,  3,  7, 12, 16};

int eeemcal_16i_channel_b_map[16] = {20, 24, 29, 33, 18, 22, 27, 31,
                                     19, 23, 28, 32, 21, 25, 30, 34};

int eeemcal_16i_channel_c_map[16] = {67, 63, 59, 55, 69, 65, 61, 57,
                                     70, 66, 60, 56, 68, 64, 58, 54};
             
int eeemcal_16i_channel_d_map[16] = {50, 46, 40, 36, 52, 48, 42, 38,
                                     51, 47, 43, 39, 49, 45, 41, 37};
          
int *eeemcal_16i_channel_map[4] = {eeemcal_16i_channel_a_map, eeemcal_16i_channel_b_map, eeemcal_16i_channel_c_map, eeemcal_16i_channel_d_map};

int eeemcal_4x4_channel_a_map[4] = {0, 4, 9, 12};
int eeemcal_4x4_channel_b_map[4] = {20, 24, 27, 31};
int eeemcal_4x4_channel_c_map[4] = {58, 62, 65, 69};
int eeemcal_4x4_channel_d_map[4] = {38, 42, 48, 52};
int *eeemcal_4x4_channel_map[4] = {eeemcal_4x4_channel_a_map, eeemcal_4x4_channel_b_map, eeemcal_4x4_channel_c_map, eeemcal_4x4_channel_d_map};

int eeemcal_16p_channel_map[4] = {6, 26, 63, 46};

int sipms_per_crystal[3] = {16, 4, 1};
int crystal_ID[25] = {5, 10, 15, 20, 25,
                      4, 9, 14, 19, 24,
                      3, 8, 13, 18, 23,
                      2, 7, 12, 17, 22,
                      1, 6, 11, 16, 21};


void single_crystal_ADC_sum(int run_number) {
    int readout = 0;
    gStyle->SetOptStat(0);
    auto path = getenv("OUTPUT_PATH");
    TFile *file = new TFile(Form("%s/run%03d.root", path, run_number));
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

    uint adc[576][20];
    uint tot[576][20];
    tree->SetBranchAddress("adc", &adc);
    tree->SetBranchAddress("tot", &tot);

    std::vector <TH1D*> crystal_sums;
    for (int crystal = 0; crystal < 25; crystal++) {
        TH1D *adc_total_sum = new TH1D(Form("crystal_%02d_sum", crystal), Form("Crystal %d ADC Sum;ADC;Counts", crystal), 256 * sipms_per_crystal[readout], 0, 1024 * sipms_per_crystal[readout]);
        crystal_sums.push_back(adc_total_sum);
    }
    TH1D *center_calo_sum = new TH1D("center_calo_sum", "Center Calorimeter ADC Sum;ADC;Counts", 256 * sipms_per_crystal[readout], 0, 1024 * sipms_per_crystal[readout]);
    TH1D *full_calo_sum = new TH1D("full_calo_sum", "Full Calorimeter ADC Sum;ADC;Counts", 256 * sipms_per_crystal[readout], 0, 1024 * sipms_per_crystal[readout]);

    for (int event = 0; event < tree->GetEntries(); event++) {
        tree->GetEntry(event);
        
        int center_sum = 0;
        int event_sum = 0;
        for (int crystal = 0; crystal < 25; crystal++) {
            int adc_sum = 0;
            int crystal_fpga = eeemcal_fpga_map[crystal];
            int crystal_asic = eeemcal_asic_map[crystal];
            int crystal_connector = eeemcal_connector_map[crystal];

            for (int channel = 0; channel < 16; channel++) {
                int crystal_channel = 144 * crystal_fpga + 72 * crystal_asic + eeemcal_16i_channel_map[crystal_connector][channel];

                int max_adc = 0;
                int max_sample = 0;
                for (int sample = 0; sample < 20; sample++) {
                    int sample_adc = adc[crystal_channel][sample] - adc[crystal_channel][0];
                    if (sample_adc > max_adc) {
                        max_adc = sample_adc;
                        max_sample = sample;
                    }
                }
                // add the neighboring two samples to max_adc
                // if (max_sample > 0) {
                //     max_adc += adc[crystal_channel][max_sample-1] - adc[crystal_channel][0];
                // }
                // if (max_sample < 19) {
                    //     max_adc += adc[crystal_channel][max_sample+1] - adc[crystal_channel][0];
                    // }
                adc_sum += max_adc;
                if (crystal == 6 || crystal == 7 || crystal == 8 || crystal == 11 || crystal == 12 || crystal == 13 || crystal == 16 || crystal == 17 || crystal == 18) {
                    center_sum += max_adc;
                }
                event_sum += max_adc;
            }
            crystal_sums[crystal]->Fill(adc_sum);
        }
        center_calo_sum->Fill(center_sum);
        full_calo_sum->Fill(event_sum);
    }
    

    
    int lower_range = 200 * sipms_per_crystal[readout];
    int upper_range = 900 * sipms_per_crystal[readout];

    double max_value = 0;
    for (int crystal = 0; crystal < 25; crystal++) {
        TF1 *fit = new TF1("fit", "gaus", lower_range, upper_range);
        auto result = crystal_sums[crystal]->Fit("fit", "QR");
        
        if (fit->Eval(fit->GetParameter(1)) > max_value) {
            max_value = fit->Eval(fit->GetParameter(1));
        }
    }
    std::cout << "max is " << max_value << std::endl;

    // Draw individual crystal sums
    TCanvas *c = new TCanvas("c", "c", 1600, 1200);
    c->cd(0);
    auto label = new TLatex();
    label->SetNDC();
    label->SetTextSize(0.05);
    label->DrawLatex(0.05, 0.9, Form("Crystal ADC Sums Run %d", run_number));


    auto pad = new TPad("pad", "pad", 0.05, 0.05, 0.95, 0.85);
    pad->Draw();
    // Add text to the top of the pad with the run and event number
    pad->cd();
    pad->Divide(5, 5, 0.000, 0.000);
    
    for (int crystal = 0; crystal < 25; crystal++) {
        pad->cd(crystal+1);
        auto fit = crystal_sums[crystal]->GetFunction("fit");
        
        crystal_sums[crystal]->SetTitle("");
        crystal_sums[crystal]->Draw("e");
        crystal_sums[crystal]->SetMaximum(max_value * 5);
        crystal_sums[crystal]->GetXaxis()->SetLabelSize(0.06);
        crystal_sums[crystal]->GetYaxis()->SetLabelSize(0.06);
        crystal_sums[crystal]->GetXaxis()->SetTitle("");
        crystal_sums[crystal]->GetYaxis()->SetTitle("");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(33);
        latex.DrawLatex(0.95, 0.95, Form("Crystal %d", crystal_ID[crystal]));
        if (fit) {
            int entries_in_range = crystal_sums[crystal]->Integral(crystal_sums[crystal]->FindBin(lower_range), crystal_sums[crystal]->FindBin(upper_range));
            double mean = fit->GetParameter(1);
            double stddev = fit->GetParameter(2);
            double mean_error = fit->GetParError(1);
            double stddev_error = fit->GetParError(2);
            double stddev_over_mean = stddev/mean;
            double stddev_over_mean_error = stddev_over_mean * sqrt(pow(mean_error/mean, 2) + pow(stddev_error/stddev, 2));
            latex.DrawLatex(0.95, 0.85, Form("Mean = %.2f#pm%.2f", mean, mean_error));
            latex.DrawLatex(0.95, 0.75, Form("StdDev = %.2f#pm%.2f", stddev, stddev_error));
            latex.DrawLatex(0.95, 0.65, Form("StdDev/Mean = %.2f#pm%.4f", stddev_over_mean, stddev_over_mean_error));
            latex.DrawLatex(0.95, 0.55, Form("Entries in range = %d", entries_in_range));
        }
    }

    c->cd(0);
    label->SetTextSize(0.03);
    label->SetTextAlign(33);
    label->SetTextAngle(90);
    label->DrawLatex(0.04, 0.85, "Counts/4 ADC");
    label->SetTextAngle(0);
    label->DrawLatex(0.925, 0.05, "ADC");
    c->SaveAs(Form("output/run%03d_adc_sum.pdf(", run_number));

    // Draw center 9 crystal sum
    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
    auto fit = new TF1("fit", "gaus", 6000, 10000);
    center_calo_sum->Fit("fit", "QR");
    center_calo_sum->SetTitle("Central 9 Crystals");
    center_calo_sum->Draw("e");
    double mean = fit->GetParameter(1);
    double stddev = fit->GetParameter(2);
    double mean_error = fit->GetParError(1);
    double stddev_error = fit->GetParError(2);
    double stddev_over_mean = stddev/mean;
    double stddev_over_mean_error = stddev_over_mean * sqrt(pow(mean_error/mean, 2) + pow(stddev_error/stddev, 2));

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(33);
    latex.DrawLatex(0.89, 0.85, Form("Mean = %.2f#pm%.2f", mean, mean_error));
    latex.DrawLatex(0.89, 0.8, Form("StdDev = %.2f#pm%.2f", stddev, stddev_error));
    latex.DrawLatex(0.89, 0.75, Form("StdDev/Mean = %.2f#pm%.4f", stddev_over_mean, stddev_over_mean_error));
    c2->SaveAs(Form("output/run%03d_adc_sum.pdf", run_number));

    // Draw full calo sum
    TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1200);
    fit = new TF1("fit", "gaus", 8000, 16000);
    full_calo_sum->Fit("fit", "QR");
    full_calo_sum->SetTitle("Full Calorimeter");
    full_calo_sum->Draw("e");
    mean = fit->GetParameter(1);
    stddev = fit->GetParameter(2);
    mean_error = fit->GetParError(1);
    stddev_error = fit->GetParError(2);
    stddev_over_mean = stddev/mean;
    stddev_over_mean_error = stddev_over_mean * sqrt(pow(mean_error/mean, 2) + pow(stddev_error/stddev, 2));

    latex.DrawLatex(0.89, 0.85, Form("Mean = %.2f#pm%.2f", mean, mean_error));
    latex.DrawLatex(0.89, 0.8, Form("StdDev = %.2f#pm%.2f", stddev, stddev_error));
    latex.DrawLatex(0.89, 0.75, Form("StdDev/Mean = %.2f#pm%.4f", stddev_over_mean, stddev_over_mean_error));
    c3->SaveAs(Form("output/run%03d_adc_sum.pdf)", run_number));
    
}