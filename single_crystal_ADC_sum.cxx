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


double get_max_ADC(uint adc[576][20], int channel) {
    int single_adc = 0;
    int max_sample = 0;
    for (int sample = 0; sample < 20; sample++) {
        int sample_adc = adc[channel][sample] - adc[channel][0];
        if (sample_adc > single_adc) {
            single_adc = sample_adc;
            max_sample = sample;
        }
    }
    return  single_adc;
}

double get_full_waveform_sum(uint adc[576][20], int channel) {
    double sum = 0;
    for (int sample = 2; sample < 10; sample++) {
        double this_sample = (double)adc[channel][sample] - (double)adc[channel][0];
        if (this_sample < 0) {
            this_sample = 0;
        }
        sum += this_sample;
    }
    // std::cout << sum << std::endl;
    return sum;
}

double decode_toa_sample(uint adc[576][20], uint toa[576][20], int channel) {
    int toa_sample = 0;
    int toa_found = 0;
    for (int sample = 0; sample < 20; sample++) {
        if (toa[channel][sample] > 0) {
            toa_sample = sample;
            toa_found++;
        }
    }
    if (toa_sample == 0) {
        return 0;
    }
    std::cout << "toa sample: " << toa_sample << std::endl;
    if (toa_found > 1) {
        std::cerr << "MORE THAN ONE TOA FOUND" << std::endl;
    }
    return adc[channel][toa_sample] - adc[channel][0];
}

double decode_tot_sample(uint adc[576][20], uint tot[576][20], int channel) {
    int tot_sample = 0;
    int tot_found = 0;
    for (int sample = 0; sample < 20; sample++) {
        if (tot[channel][sample] > 0) {
            tot_sample = sample;
            tot_found++;
        }
    }
    if (tot_sample == 0) {
        return 0;
    }
    std::cout << "ToT sample: " << tot_sample << std::endl;
    if (tot_found > 1) {
        std::cerr << "MORE THAN ONE TOT FOUND" << std::endl;
    }
    return adc[channel][tot_sample] - adc[channel][0];
}

void single_crystal_ADC_sum(int run_number) {
    int readout = 0;
    gStyle->SetOptStat(0);
    auto path = getenv("OUTPUT_PATH");
    TFile *file = new TFile(Form("%s/Run%03d.root", path, run_number));
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
    uint toa[576][20]; 
    tree->SetBranchAddress("adc", &adc);
    tree->SetBranchAddress("tot", &tot);
    tree->SetBranchAddress("toa", &toa);

    std::vector<TH1D*> sipm_single_sums;
    std::vector<TH1D*> sipm_full_sums;
    for (int crystal = 0; crystal < 25; crystal++) {
        for (int sipm = 0; sipm < sipms_per_crystal[readout]; sipm++) {
            TH1D *sipm_sum = new TH1D(Form("crystal_%02d_sipm_%02d_sum_single", crystal, sipm), Form("Crystal %d SiPM %d ADC Sum;ADC;Counts", crystal, sipm), 256, 0, 1024);
            sipm_single_sums.push_back(sipm_sum);
            sipm_sum = new TH1D(Form("crystal_%02d_sipm_%02d_sum_full", crystal, sipm), Form("Crystal %d SiPM %d ADC Sum;ADC;Counts", crystal, sipm), 256, 0, 3500);
            sipm_full_sums.push_back(sipm_sum);
        }
    }


    std::vector <TH1D*> crystal_single_sums;
    std::vector <TH1D*> crystal_full_sums;
    for (int crystal = 0; crystal < 25; crystal++) {
        TH1D *adc_total_sum = new TH1D(Form("crystal_%02d_sum_single", crystal), Form("Crystal %d ADC Sum;ADC;Counts", crystal), 256 * sipms_per_crystal[readout], 0, 1024 * sipms_per_crystal[readout]);
        crystal_single_sums.push_back(adc_total_sum);
        adc_total_sum = new TH1D(Form("crystal_%02d_sum_full", crystal), Form("Crystal %d ADC Sum;ADC;Counts", crystal), 25 * sipms_per_crystal[readout], 0, 2500 * sipms_per_crystal[readout]);
        crystal_full_sums.push_back(adc_total_sum);
    }
    TH1D *center_calo_single_sum = new TH1D("center_calo_single_sum_single", "Center Calorimeter ADC Sum;ADC;Counts", 256 * sipms_per_crystal[readout], 0, 1024 * sipms_per_crystal[readout]);
    TH1D *center_calo_full_sum = new TH1D("center_calo_full_sum_single", "Center Calorimeter ADC Sum;ADC;Counts", 25 * sipms_per_crystal[readout], 0, 4000 * sipms_per_crystal[readout]);
    TH1D *full_calo_single_sum = new TH1D("full_calo_single_sum_single", "Full Calorimeter ADC Sum;ADC;Counts", 256 * sipms_per_crystal[readout], 0, 1024 * sipms_per_crystal[readout]);
    TH1D *full_calo_full_sum = new TH1D("full_calo_full_sum_single", "Full Calorimeter ADC Sum;ADC;Counts", 25 * sipms_per_crystal[readout], 0, 4000 * sipms_per_crystal[readout]);

    for (int event = 0; event < tree->GetEntries(); event++) {
        tree->GetEntry(event);
        
        int center_single_sum = 0;
        int event_single_sum = 0;
        double center_full_sum = 0;
        double event_full_sum = 0;
        for (int crystal = 0; crystal < 25; crystal++) {
            int crystal_single_sum = 0;
            int crystal_full_sum = 0;
            int crystal_fpga = eeemcal_fpga_map[crystal];
            int crystal_asic = eeemcal_asic_map[crystal];
            int crystal_connector = eeemcal_connector_map[crystal];

            for (int channel = 0; channel < 16; channel++) {
                int crystal_channel = 144 * crystal_fpga + 72 * crystal_asic + eeemcal_16i_channel_map[crystal_connector][channel];
                double single_adc = get_max_ADC(adc, crystal_channel);
                // decode_toa_sample(adc, toa, crystal_channel);
                // decode_tot_sample(adc, tot, crystal_channel);
                double full_adc = get_full_waveform_sum(adc, crystal_channel);
                crystal_single_sum += single_adc;
                crystal_full_sum += full_adc;
                if (crystal == 6 || crystal == 7 || crystal == 8 || crystal == 11 || crystal == 12 || crystal == 13 || crystal == 16 || crystal == 17 || crystal == 18) {
                    center_single_sum += single_adc;
                    center_full_sum += full_adc;
                }
                event_single_sum += single_adc;
                event_full_sum += full_adc;
                sipm_single_sums[crystal * sipms_per_crystal[readout] + channel]->Fill(single_adc);
                sipm_full_sums[crystal * sipms_per_crystal[readout] + channel]->Fill(full_adc);
            }
            crystal_single_sums[crystal]->Fill(crystal_single_sum);
            crystal_full_sums[crystal]->Fill(crystal_full_sum);
        }
        center_calo_single_sum->Fill(center_single_sum);
        center_calo_full_sum->Fill(center_full_sum);
        // std::cout << center_full_sum << std::endl;
        full_calo_single_sum->Fill(event_single_sum);
        full_calo_full_sum->Fill(event_full_sum);
        // std::cout << event_full_sum << std::endl;
    }
    

    
    int lower_range = 200 * sipms_per_crystal[readout];
    int upper_range = 900 * sipms_per_crystal[readout];

    double max_value = 0;
    for (int crystal = 0; crystal < 25; crystal++) {
        TF1 *fit = new TF1("fit", "gaus", lower_range, upper_range);
        auto result = crystal_single_sums[crystal]->Fit("fit", "QR");
        
        if (fit->Eval(fit->GetParameter(1)) > max_value) {
            max_value = fit->Eval(fit->GetParameter(1));
        }
    }
    std::cout << "max is " << max_value << std::endl;


    // SINGLE ADC SUMS
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
        auto fit = crystal_single_sums[crystal]->GetFunction("fit");
        
        crystal_single_sums[crystal]->SetTitle("");
        crystal_single_sums[crystal]->Draw("e");
        crystal_single_sums[crystal]->SetMaximum(max_value * 5);
        crystal_single_sums[crystal]->GetXaxis()->SetLabelSize(0.06);
        crystal_single_sums[crystal]->GetYaxis()->SetLabelSize(0.06);
        crystal_single_sums[crystal]->GetXaxis()->SetTitle("");
        crystal_single_sums[crystal]->GetYaxis()->SetTitle("");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(33);
        latex.DrawLatex(0.95, 0.95, Form("Crystal %d", crystal_ID[crystal]));
        if (fit) {
            int entries_in_range = crystal_single_sums[crystal]->Integral(crystal_single_sums[crystal]->FindBin(lower_range), crystal_single_sums[crystal]->FindBin(upper_range));
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
    c->SaveAs(Form("output/Run%03d_adc_single_sum.pdf(", run_number));

    // Draw center 9 crystal sum
    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
    auto fit = new TF1("fit", "gaus", 6000, 10000);
    center_calo_single_sum->Fit("fit", "QR");
    center_calo_single_sum->SetTitle("Central 9 Crystals");
    center_calo_single_sum->Draw("e");
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
    c2->SaveAs(Form("output/Run%03d_adc_single_sum.pdf", run_number));

    // Draw full calo sum
    TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1200);
    fit = new TF1("fit", "gaus", 8000, 16000);
    full_calo_single_sum->Fit("fit", "QR");
    full_calo_single_sum->SetTitle("Full Calorimeter");
    full_calo_single_sum->Draw("e");
    mean = fit->GetParameter(1);
    stddev = fit->GetParameter(2);
    mean_error = fit->GetParError(1);
    stddev_error = fit->GetParError(2);
    stddev_over_mean = stddev/mean;
    stddev_over_mean_error = stddev_over_mean * sqrt(pow(mean_error/mean, 2) + pow(stddev_error/stddev, 2));

    latex.DrawLatex(0.89, 0.85, Form("Mean = %.2f#pm%.2f", mean, mean_error));
    latex.DrawLatex(0.89, 0.8, Form("StdDev = %.2f#pm%.2f", stddev, stddev_error));
    latex.DrawLatex(0.89, 0.75, Form("StdDev/Mean = %.2f#pm%.4f", stddev_over_mean, stddev_over_mean_error));
    c3->SaveAs(Form("output/Run%03d_adc_single_sum.pdf", run_number));


    // Each crystal, per SiPM
    latex.SetTextSize(0.06);   
    for (int crystal = 0; crystal < 25; crystal++) {
        TCanvas *canvas = new TCanvas(Form("crystal_%02d_sipm_sum", crystal), Form("crystal_%02d_sipm_sum", crystal_ID[crystal]), 1600, 1200);
        canvas->cd(0);
        auto label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.05);
        label->DrawLatex(0.05, 0.9, Form("Crystal %d SiPM ADC Sums Run %d", crystal_ID[crystal], run_number));
        auto pad = new TPad("pad", "pad", 0.05, 0.05, 0.95, 0.85);
        pad->Draw();
        pad->cd();
        pad->Divide(4, 4, 0.000, 0.000);
        for (int sipm = 0; sipm < sipms_per_crystal[readout]; sipm++) {
            pad->cd(sipm+1);
            auto fit = new TF1("fit", "gaus", 200, 900);
            sipm_single_sums[crystal * sipms_per_crystal[readout] + sipm]->Fit("fit", "QR");
            sipm_single_sums[crystal * sipms_per_crystal[readout] + sipm]->Draw("e");
            int entries_in_range = sipm_single_sums[crystal * sipms_per_crystal[readout] + sipm]->Integral(sipm_single_sums[crystal * sipms_per_crystal[readout] + sipm]->FindBin(200), sipm_single_sums[crystal * sipms_per_crystal[readout] + sipm]->FindBin(900));
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
        canvas->SaveAs(Form("output/Run%03d_adc_single_sum.pdf", run_number));
    }
    auto end_page = new TCanvas("end_page", "end_page", 1600, 1200);
    end_page->SaveAs(Form("output/Run%03d_adc_single_sum.pdf)", run_number));



    lower_range = 1150 * sipms_per_crystal[readout];
    upper_range = 1800 * sipms_per_crystal[readout];

    max_value = 0;
    for (int crystal = 0; crystal < 25; crystal++) {
        TF1 *fit = new TF1("fit", "gaus", lower_range, upper_range);
        auto result = crystal_full_sums[crystal]->Fit("fit", "QR");
        
        if (fit->Eval(fit->GetParameter(1)) > max_value) {
            max_value = fit->Eval(fit->GetParameter(1));
        }
    }
    std::cout << "max is " << max_value << std::endl;


    // Full ADC plots
    // Draw individual crystal sums
    c = new TCanvas("c5", "c", 1600, 1200);
    c->cd(0);
    label->SetNDC();
    label->SetTextSize(0.05);
    label->DrawLatex(0.05, 0.9, Form("Crystal ADC Sums Run %d", run_number));


    pad = new TPad("pad2", "pad", 0.05, 0.05, 0.95, 0.85);
    pad->Draw();
    // Add text to the top of the pad with the run and event number
    pad->cd();
    pad->Divide(5, 5, 0.000, 0.000);
    
    for (int crystal = 0; crystal < 25; crystal++) {
        pad->cd(crystal+1);
        auto fit = crystal_full_sums[crystal]->GetFunction("fit");
        
        crystal_full_sums[crystal]->SetTitle("");
        crystal_full_sums[crystal]->Draw("e");
        crystal_full_sums[crystal]->SetMaximum(max_value * 5);
        crystal_full_sums[crystal]->GetXaxis()->SetLabelSize(0.06);
        crystal_full_sums[crystal]->GetYaxis()->SetLabelSize(0.06);
        crystal_full_sums[crystal]->GetXaxis()->SetTitle("");
        crystal_full_sums[crystal]->GetYaxis()->SetTitle("");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(33);
        latex.DrawLatex(0.95, 0.95, Form("Crystal %d", crystal_ID[crystal]));
        if (fit) {
            int entries_in_range = crystal_full_sums[crystal]->Integral(crystal_full_sums[crystal]->FindBin(lower_range), crystal_full_sums[crystal]->FindBin(upper_range));
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
    c->SaveAs(Form("output/Run%03d_adc_full_sum.pdf(", run_number));

    // Draw center 9 crystal sum
    c2 = new TCanvas("c6", "c2", 1600, 1200);
    fit = new TF1("fit", "gaus", 25000, 38000);
    center_calo_full_sum->Fit("fit", "QR");
    center_calo_full_sum->SetTitle("Central 9 Crystals");
    center_calo_full_sum->Draw("e");
    mean = fit->GetParameter(1);
    stddev = fit->GetParameter(2);
    mean_error = fit->GetParError(1);
    stddev_error = fit->GetParError(2);
    stddev_over_mean = stddev/mean;
    stddev_over_mean_error = stddev_over_mean * sqrt(pow(mean_error/mean, 2) + pow(stddev_error/stddev, 2));

    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(33);
    latex.DrawLatex(0.89, 0.85, Form("Mean = %.2f#pm%.2f", mean, mean_error));
    latex.DrawLatex(0.89, 0.8, Form("StdDev = %.2f#pm%.2f", stddev, stddev_error));
    latex.DrawLatex(0.89, 0.75, Form("StdDev/Mean = %.2f#pm%.4f", stddev_over_mean, stddev_over_mean_error));
    c2->SaveAs(Form("output/Run%03d_adc_full_sum.pdf", run_number));

    // Draw full calo sum
    c3 = new TCanvas("c7", "c3", 1600, 1200);
    fit = new TF1("fit", "gaus", 30000, 45000);
    full_calo_full_sum->Fit("fit", "QR");
    full_calo_full_sum->SetTitle("Full Calorimeter");
    full_calo_full_sum->Draw("e");
    mean = fit->GetParameter(1);
    stddev = fit->GetParameter(2);
    mean_error = fit->GetParError(1);
    stddev_error = fit->GetParError(2);
    stddev_over_mean = stddev/mean;
    stddev_over_mean_error = stddev_over_mean * sqrt(pow(mean_error/mean, 2) + pow(stddev_error/stddev, 2));

    latex.DrawLatex(0.89, 0.85, Form("Mean = %.2f#pm%.2f", mean, mean_error));
    latex.DrawLatex(0.89, 0.8, Form("StdDev = %.2f#pm%.2f", stddev, stddev_error));
    latex.DrawLatex(0.89, 0.75, Form("StdDev/Mean = %.2f#pm%.4f", stddev_over_mean, stddev_over_mean_error));
    c3->SaveAs(Form("output/Run%03d_adc_full_sum.pdf", run_number));


    // Each crystal, per SiPM
    latex.SetTextSize(0.06);   
    for (int crystal = 0; crystal < 25; crystal++) {
        TCanvas *canvas = new TCanvas(Form("crystal_%02d_sipm_full_sum", crystal), Form("crystal_%02d_sipm_sum", crystal_ID[crystal]), 1600, 1200);
        canvas->cd(0);
        auto label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.05);
        label->DrawLatex(0.05, 0.9, Form("Crystal %d SiPM ADC Sums Run %d", crystal_ID[crystal], run_number));
        auto pad = new TPad("pad", "pad", 0.05, 0.05, 0.95, 0.85);
        pad->Draw();
        pad->cd();
        pad->Divide(4, 4, 0.000, 0.000);
        for (int sipm = 0; sipm < sipms_per_crystal[readout]; sipm++) {
            pad->cd(sipm+1);
            auto fit = new TF1("fit", "gaus", 600, 3500);
            sipm_full_sums[crystal * sipms_per_crystal[readout] + sipm]->Fit("fit", "QR");
            sipm_full_sums[crystal * sipms_per_crystal[readout] + sipm]->Draw("e");
            int entries_in_range = sipm_full_sums[crystal * sipms_per_crystal[readout] + sipm]->Integral(sipm_full_sums[crystal * sipms_per_crystal[readout] + sipm]->FindBin(200), sipm_single_sums[crystal * sipms_per_crystal[readout] + sipm]->FindBin(900));
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
        canvas->SaveAs(Form("output/Run%03d_adc_full_sum.pdf", run_number));
    }
    end_page->SaveAs(Form("output/Run%03d_adc_full_sum.pdf)", run_number));
}