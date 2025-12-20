#include <map>
#include <vector>
#include <iosfwd>
#include <istream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLatex.h>

float sigma_cut = 1;
float energy_fraction_cut = 0.3;
long n_events = 1000000;
// n_events = 1000;

float center_x = 1.98022;
float sigma_x = 0.193676 * sigma_cut;
float center_y = 1.97743;
float sigma_y = 0.193863 * sigma_cut;

float tot_min_cut = 10000;

float adc_calib = 32444.1;  // 32444.1 Signal_ADC = 1 GeV 

std::vector<int> run_numbers = {321, 326, 331, 336};
std::vector<float> energies  = {2.0, 3.0, 4.0, 5.0};

// int run_number = 321;
// int run_number = 412;

std::map<int, std::vector<int>> read_mapping(const std::string& filename) {
    std::map<int, std::vector<int>> mapping;
    std::ifstream file(filename);
    std::string line;
    
    // Skip header
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int fpga, asic, connector, crystal, sipm, channel;
        char comma;
        
        // Parse CSV: FPGA,ASIC,Connector,Crystal,SiPM,Channel
        if (iss >> fpga >> comma >> asic >> comma >> connector >> comma >> crystal >> comma >> sipm >> comma >> channel) {
            // Ensure the vector is large enough for all SiPMs (16 per crystal)
            if (mapping[crystal].size() <= sipm) {
                mapping[crystal].resize(sipm + 1, -1);
            }
            mapping[crystal][sipm] = channel;
        }
    }
    
    file.close();
    return mapping;
}

float calculate_signal(uint32_t *adc_values, float gain) {
    // Pedestal is the mean of the first three samples
    float pedestal = (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;
    
    // Signal is the sum of the three greatest values
    std::vector<float> samples;
    for (int i = 3; i < 20; ++i) {
        samples.push_back(adc_values[i] - pedestal);
    }
    std::sort(samples.begin(), samples.end(), std::greater<float>());
    float signal = samples[0] + samples[1] + samples[2] + samples[3];
    return signal * gain;
}

bool calculate_signal(uint32_t *adc_values, uint32_t *tot_values, float gain, float &signal) {
    // Check if there is a ToT value
    uint32_t tot = 0;
    for (int i = 0; i < 20; ++i) {
        if (tot_values[i] > 0) {
            tot = tot_values[i];
            break;
        }
    }
    // If there is no ToT value, return the single SiPM ADC signal
    if (tot == 0) {
        signal = calculate_signal(adc_values, gain);
        return false;
    }

    // Otherwise, just return the first ToT value
    signal = tot;
    return true;
}

bool is_tot(uint32_t *tot_values) {
    for (int i = 0; i < 20; ++i) {
        // std::cout << "TOT Value[" << i << "]: " << tot_values[i] << std::endl;
        if (tot_values[i] > 50) {
            return true;
        }
    }
    return false;
}

bool position_cut(float x, float y) {
    // Cut anything outside of the center
    if (std::abs(x - center_x) > sigma_x || std::abs(y - center_y) > sigma_y) {
        return false;
    }
    return true;
}

bool calculate_cog(TH2* distribution, float *values) {
    float total_signal = 0;
    float x_weighted_sum = 0;
    float y_weighted_sum = 0;
    float x_cog = 0;
    float y_cog = 0;
    float w = 4;

    for (int i = 0; i < 25; i++) {
        total_signal += values[i];
    }
    
    float total_weight = 0;
    for (int i = 0; i < 25; i++) {
        int x = i % 5;
        int y = i / 5;
        float signal = values[i];
        float weight = w + std::log(signal / total_signal); // Avoid log(0)
        if (weight < 0) {
            weight = 0;
        }
        total_weight += weight;
        x_weighted_sum += x * weight;
        y_weighted_sum += y * weight;

    }
    if (total_weight > 0) {
        x_cog = x_weighted_sum / total_weight;
        y_cog = y_weighted_sum / total_weight;
        distribution->Fill(x_cog, y_cog);
    }
    return position_cut(x_cog, y_cog);
}

void print_progress(int progress) {
    std::cout << " [";
    for (int i = 0; i < 25; i++) {
        if (i < progress) {
            std::cout << "*";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "]\r" << std::flush;
}

void tot_calibration() {
    run_numbers.clear();
    energies.clear();
    for (int r = 319; r <= 338; r++) {
        run_numbers.push_back(r);
        energies.push_back(1.6 + 0.2 * (r - 319));
    }

    gStyle->SetOptStat(0);
    auto mapping = read_mapping("eeemcal_desy_dec2025_mapping.csv");

    TH1* gain_factor = nullptr;
    TH1* crystal_gain_factor;
    TFile* gain_file = TFile::Open("output/gain_factors.root");
    if (gain_file && !gain_file->IsZombie()) {
        gain_factor = (TH1*)gain_file->Get("gain_factors");
        crystal_gain_factor = (TH1*)gain_file->Get("crystal_factor");
        std::cout << "Loaded gain factors from file." << std::endl;
    }
    if (!gain_factor) {
        gain_factor = new TH1F("gain_factors", "Gain Factor", 400, 0, 400);
        for (int i = 1; i <= 400; i++) {
            gain_factor->SetBinContent(i, 1.0);
        }
    }
    if (!crystal_gain_factor) {
        crystal_gain_factor = new TH1F("crystal_factor", "Crystal Gain", 25, 0, 25);
        for (int i = 1; i <= 25; i++) {
            crystal_gain_factor->SetBinContent(i, 1);
        }
    }

    std::vector<TH1*> missing_energies;
    std::vector<TH1*> tot_distributions;
    std::vector<TH1*> tot_energies;
    std::vector<TH2*> cog_distributions;

    TH2 *total_distribution =new TH2F("run%d_tot_energy", "ToT Energy;Energy (GeV);ToT Signal", 100, 0, 5, 1000, 0, 60000);
    TH2 *total_distribution_invt =new TH2F("run%d_tot_energy_invt", "Energy as a function of ToT;ToT Signal;Energy (GeV);", 1000, 0, 60000, 100, 0, 5);



    for (int run = 0; run < run_numbers.size(); run++) {
        int run_number = run_numbers[run];
        float energy = energies[run];
        // Process data file
        char file_path[256];
        sprintf(file_path, "/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number);
        TFile* root_file = TFile::Open(file_path);
        TTree* tree = (TTree*)root_file->Get("events");
        uint32_t adc[576][20];
        uint32_t tot[576][20];
        int tot_events = 0;
        tree->SetBranchAddress("adc", &adc);
        tree->SetBranchAddress("tot", &tot);
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("adc", 1);
        tree->SetBranchStatus("tot", 1);

        
        missing_energies.push_back(new TH1F(Form("run%d_missing_energy", run_number), Form("Run %d ADC portion of energy", run_number), 100, 0, 6));
        tot_distributions.push_back(new TH1F(Form("run%d_tot_distribution", run_number), Form("Run %d ToT;ToT;Events", run_number), 1024, 0, 16 * 4096));
        tot_energies.push_back(new TH2F(Form("run%d_tot_energy", run_number), Form("Run %d ToT Energy;Energy (GeV);(ToT Signal)", run_number), 100, 0, 5, 1000, 0, 60000));
        cog_distributions.push_back(new TH2F(Form("run%d_cog_distribution", run_number), Form("Run %d Center of Gravity Distribution;X (# Crystals);Y (# Crystals)", run_number), 100, -0.5, 4.5, 100, -0.5, 4.5));

        Long64_t nentries = tree->GetEntries();
        if (n_events < nentries) {
            nentries = n_events;
        }
        std::cout << "Run " << run_number << std::endl;
        std::cout << "Processing " << nentries << " events" << std::endl;
        int complete = 0;
        for (Long64_t entry = 0; entry < nentries; ++entry) {
            if (entry * 25 / nentries > complete) {
                complete = entry * 25 / nentries;
                print_progress(complete);
            }

            tree->GetEntry(entry);
            float central_signal = 0.0f;
            float central_nine_signal = 0.0f;
            float total_signal = 0.0f;
            bool is_tot_event = false;

            float signals[25];
            for (int crystal = 0; crystal < 25; crystal++) {
                if (crystal == 9) {
                    continue;
                }
                if (crystal == 12) {
                    continue;
                }
                float crystal_signal = 0.0f;
                for (int sipm = 0; sipm < 16; sipm++) {
                    int channel = mapping[crystal][sipm];
                    float gain = gain_factor->GetBinContent(crystal * 16 + sipm + 1);
                    float channel_signal = calculate_signal(adc[channel], gain);
                    crystal_signal += channel_signal;
                    if (!is_tot_event && is_tot(tot[channel])) {
                        is_tot_event = true;
                    }
                }
                signals[crystal] = crystal_signal * crystal_gain_factor->GetBinContent(crystal + 1);
            }
            if (is_tot_event) {
                tot_events++;
                // continue;
            }

            // Get the ToT
            float center_signal = 0;
            for (int sipm = 0; sipm < 16; sipm++) {
                float signal = 0;
                int channel = mapping[12][sipm];
                float gain = gain_factor->GetBinContent(12 * 16 + sipm + 1);
                calculate_signal(adc[channel], tot[channel], gain, signal);
                center_signal += signal * gain;
            }
            if (center_signal < tot_min_cut) {
                center_signal = 0;
            }


            float x_cog, y_cog;
            bool keep = calculate_cog(cog_distributions[run], signals);
            keep &= (!is_tot_event);    // Get rid of events with a tot in a non-central channel
            if (!keep) {
                continue;
            }

            tot_distributions[run]->Fill(center_signal);

            // Sum up the total energy found in the non-central crystals
            float non_central_energy = 0;
            for (int crystal = 0; crystal < 25; crystal++) {
                if (crystal == 12) {continue;}
                non_central_energy += signals[crystal];
            }
            // Calibrate to GeV
            non_central_energy /= adc_calib;
            missing_energies[run]->Fill(non_central_energy);
            tot_energies[run]->Fill(energy - non_central_energy, center_signal);
            total_distribution->Fill(energy - non_central_energy, center_signal);
            total_distribution_invt->Fill(center_signal, energy - non_central_energy);
        }
        std::cout << std::endl;
    }


    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    text->SetTextFont(42);

    TCanvas* canvas = new TCanvas("tot_calib", "", 800, 600);
    canvas->SaveAs("output/tot_calib.pdf(");
    for (int run = 0; run < run_numbers.size(); run++) {
        canvas->Clear();
        canvas->Divide(2, 2);
        canvas->cd(1);
        missing_energies[run]->Draw();
        text->SetTextAlign(31);
        text->DrawLatex(0.85, 0.85, Form("Mean: %.2f", missing_energies[run]->GetMean()));
        text->DrawLatex(0.85, 0.80, Form("ADC percentage: %.2f \%", 100 * missing_energies[run]->GetMean() / energies[run]));

        canvas->cd(2);
        tot_distributions[run]->Draw();
        canvas->cd(3);
        tot_energies[run]->Draw("colz");
        canvas->cd(4);
        cog_distributions[run]->Draw("colz");
        canvas->SaveAs("output/tot_calib.pdf");
    }
    canvas->Clear();
    total_distribution->Draw("colz");
    canvas->SaveAs("output/tot_calib.pdf");


    TF1 *range_one_fit = new TF1("range_one", "pol1", 35000, 50000);
    total_distribution_invt->Fit(range_one_fit, "R");
    range_one_fit->SetLineColor(kRed);

    TF1 *range_two_fit = new TF1("range_two", "pol1", 20000, 26000);
    total_distribution_invt->Fit(range_two_fit, "R");
    range_two_fit->SetLineColor(kGreen + 2);



    total_distribution_invt->Draw("colz");
    range_one_fit->Draw("same");
    range_two_fit->Draw("same");
    canvas->SaveAs("output/tot_calib.pdf)");
}