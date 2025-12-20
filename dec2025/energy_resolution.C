#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <istream>
#include <iosfwd>
#include <sstream>

#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph.h>
#include <TLegend.h>

float sigma_cut = 2;
float energy_fraction_cut = 0.3;
long n_events = 1000000;
// n_events = 1000;
float tot_min_cut = 10000;

float center_x = 1.98022;
float sigma_x = 0.193676 * sigma_cut;
float center_y = 1.97743;
float sigma_y = 0.193863 * sigma_cut;

float adc_calib = 32444.1;  // 32444.1 Signal_ADC = 1 GeV 



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

void fit_peak(TH1* hist) {
    float mean = hist->GetMean();
    float rms = hist->GetRMS();
    float fit_min = mean - 1.5 * rms;
    float fit_max = mean + 1.5 * rms;
    TF1* rough_fit = new TF1("rough_fit", "gaus", fit_min, fit_max);
    hist->Fit(rough_fit, "RQ");
    float peak = rough_fit->GetParameter(1);
    float sigma = rough_fit->GetParameter(2);
    rough_fit->SetLineColor(kBlue);
    rough_fit->SetLineStyle(2);

    TF1* second_fit = new TF1("second_fit", "gaus", peak - sigma, peak + sigma);
    hist->Fit(second_fit, "RQ");
    peak = second_fit->GetParameter(1);
    sigma = second_fit->GetParameter(2);

    TF1* final_fit = new TF1("final_fit", "crystalball", peak - sigma, peak + sigma);
    final_fit->SetParameters(second_fit->GetParameter(0), peak, sigma, 1.5, 2.0);
    // TF1* final_fit = new TF1("final_fit", "gaus", peak - sigma, peak + sigma);
    hist->Fit(final_fit, "R");
}

void draw_text(TF1* fit, int run_number, float beam_energy) {
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    text->SetTextFont(42);
    text->SetTextAlign(31);

    float peak = fit->GetParameter(1);;
    float sigma = fit->GetParameter(2);
    float resolution = (sigma / peak) * 100.0f;
    text->DrawLatex(0.93, 0.85, Form("%.01f GeV Electrons", beam_energy));
    text->DrawLatex(0.93, 0.80, Form("Run %d", run_number));
    text->DrawLatex(0.93, 0.75, Form("Peak: %.03f", peak));
    text->DrawLatex(0.93, 0.70, Form("Sigma: %.03f", sigma));
    text->DrawLatex(0.93, 0.65, Form("Resolution: %.2f%%", resolution));
    text->DrawLatex(0.93, 0.60, "Signal method 2");

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

void energy_resolution() {
    gStyle->SetOptStat(0);
    
    std::vector<int> run_numbers = {316, 321, 326, 331, 336};
    std::vector<float> energies  = {1.0, 2.0, 3.0, 4.0, 5.0};
    // run_numbers.clear();
    // energies.clear();
    // for (int r = 316; r <= 338; r++) {
    //     run_numbers.push_back(r);
    //     energies.push_back(1. + 0.2 * (r - 316));
    // }


    int central_crystal_index = 12;
    int center_nine_indexes[8] = {7, 8, 9, 11, 13, 17, 18, 19};
    int remaining_indexes[16] = {0, 1, 2, 3, 4, 5, 6, 10, 11, 15, 16, 20, 21, 22, 23, 24};

    std::vector<TH1*> central_crystal_energy_vec;
    std::vector<TH1*> central_nine_energy_vec;
    std::vector<TH1*> total_energy_vec;
    std::vector<TH2*> cog_distribution_vec;

    std::vector<std::vector<TH1*>> crystal_energy;
    std::vector<std::vector<TH1*>> crystal_energy_shares;

    // Set up all the histograms we need
    for (int run = 0; run < run_numbers.size(); run++) {
        int run_number = run_numbers[run];
        central_crystal_energy_vec.push_back(new TH1F(Form("run%d_central_crystal_energy", run_number), Form("Run %d Central Crystal Energy;Energy (ADC);Events", run_number), 500, 0, 8));
        central_nine_energy_vec.push_back(new TH1F(Form("run%d_central_nine_energy", run_number), Form("Run %d Central 3x3 Energy;Energy (ADC);Events", run_number), 500, 0, 8));
        total_energy_vec.push_back(new TH1F(Form("run%d_total_energy", run_number), Form("Run %d Total Energy;Energy (ADC);Events", run_number), 500, 0, 8));
        cog_distribution_vec.push_back(new TH2F(Form("run%d_cog_distribution", run_number), Form("Run %d Center of Gravity Distribution;X (# Crystals);Y (# Crystals)", run_number), 100, -0.5, 4.5, 100, -0.5, 4.5));
    
        crystal_energy.push_back(std::vector<TH1*>());
        crystal_energy_shares.push_back(std::vector<TH1*>());
        for (int i = 0; i < 25; i++) { 
            crystal_energy[run].push_back(new TH1F(Form("run%d_crystal_%02d_energy", run_number, i), Form("Run %d Crystal %02d Energy;Energy (ADC);Events", run_number, i), 500, 0, 8));
            crystal_energy_shares[run].push_back(new TH1F(Form("run%d_crystal_%02d_energy_share", run_number, i), Form("Run %d Crystal %02d Energy Share;Energy Share;Events", run_number, i), 10000, 0, 1));
        }
    }

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

    for (int run = 0; run < run_numbers.size(); run++) {
        // Process data file
        int run_number = run_numbers[run];
        float energy = energies[run];
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

        Long64_t nentries = tree->GetEntries();
        if (n_events < nentries) {
            nentries = n_events;
        }
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

            for (int i = 0; i < 25; i++) {
                signals[i] /= adc_calib;
            }
    

            // Convert the ToT to energy - these are _incredibly_ rough
            float center_energy = 0;
            if (center_signal > 30000) {
                center_energy = -0.862098 + 9.68314e-05 * center_signal;
            } else {
                center_energy = -0.1 + 4.29422e-05 * center_signal;
            }
            signals[12] = center_energy;

            // std::cout << "Center energy: " << center_energy << std::endl;




            float x_cog, y_cog;
            bool keep = calculate_cog(cog_distribution_vec[run], signals);
            keep &= (!is_tot_event);
            // if (!keep) {
            //     std::cout << "skipping?" << std::endl;
            //     continue;
            // }


            // Populate energies
            central_signal = signals[central_crystal_index];
            for (int i = 0; i < 8; i++) {
                int idx = center_nine_indexes[i];
                central_nine_signal += signals[idx];
            }
            for (int i = 0; i < 16; i++) {
                int idx = remaining_indexes[i];
                total_signal += signals[idx];
            }
            if (signals[12] / total_signal < energy_fraction_cut) { // Skip events where the majority of the energy share is not in the central crystal;
                continue;
            }

            central_nine_signal += central_signal;
            total_signal += central_nine_signal;
            central_crystal_energy_vec[run]->Fill(central_signal);
            central_nine_energy_vec[run]->Fill(central_nine_signal);
            total_energy_vec[run]->Fill(total_signal);
            // std::cout << total_signal << std::endl;

            // Energy share fill
            for (int i = 0; i < 25; i++) {
                crystal_energy[run][i]->Fill(signals[i]);
                crystal_energy_shares[run][i]->Fill(signals[i] / total_signal);
            }

        }

        fit_peak(central_crystal_energy_vec[run]);
        fit_peak(central_nine_energy_vec[run]);
        fit_peak(total_energy_vec[run]);
    }


    int crystal_mapping[25] = {
        4, 9, 14, 19, 24,
        3, 8, 13, 18, 23,
        2, 7, 12, 17, 22,
        1, 6, 11, 16, 21,
        0, 5, 10, 15, 20
    };

    TH1 *dummy = new TH1F("energy_resolution_center", "Energy Resolution;Energy (GeV);#sigma/mean", 1, 0, 6);
    TGraph *res_center = new TGraph(run_numbers.size());
    TGraph *res_middle = new TGraph(run_numbers.size());
    TGraph *res_all    = new TGraph(run_numbers.size());
    res_center->SetMarkerColor(kRed);
    res_center->SetLineColor(kRed);
    res_center->SetLineStyle(kDashed);
    res_center->SetMarkerStyle(21);
    res_center->SetMarkerSize(2);

    res_middle->SetMarkerColor(kBlue);
    res_middle->SetLineColor(kBlue);
    res_middle->SetLineStyle(kDashed);
    res_middle->SetMarkerStyle(22);
    res_middle->SetMarkerSize(2);

    res_all->SetMarkerColor(kMagenta);
    res_all->SetLineColor(kMagenta);
    res_all->SetLineStyle(kDashed);
    res_all->SetMarkerStyle(23);
    res_all->SetMarkerSize(2);

    for (int run = 0; run < run_numbers.size(); run++) {
        float res = central_crystal_energy_vec[run]->GetFunction("final_fit")->GetParameter(2) / central_crystal_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
        res_center->SetPoint(run, energies[run], 100.0 * res);
        res = central_nine_energy_vec[run]->GetFunction("final_fit")->GetParameter(2) / central_nine_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
        res_middle->SetPoint(run, energies[run], 100.0 * res);
        res = total_energy_vec[run]->GetFunction("final_fit")->GetParameter(2) / total_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
        res_all->SetPoint(run, energies[run], 100.0 * res);
    }

    TCanvas* canvas = new TCanvas("gain_matching", "", 800, 600);
    dummy->Draw();
    dummy->SetMinimum(0);
    dummy->SetMaximum(10);
    res_center->Draw("same lp");
    res_middle->Draw("same lp");
    res_all->Draw("same lp");

    TLegend *l = new TLegend(0.6, 0.6, 0.89, 0.89);
    l->SetLineWidth(0);
    l->AddEntry(res_center, "Center Crystal");
    l->AddEntry(res_middle, "Central 9 Crystals");
    l->AddEntry(res_all, "All Crystals");
    l->Draw();



    
    canvas->SaveAs("output/energy_resolution.pdf(");
    for (int run = 0; run < run_numbers.size(); run++) {
        canvas->Clear();
        int run_number = run_numbers[run];
        float energy = energies[run];

        canvas->SetRightMargin(0.05);
        central_crystal_energy_vec[run]->Draw("HIST e");
        auto fit = central_crystal_energy_vec[run]->GetFunction("final_fit");
        fit->Draw("same");
        draw_text(fit, run_number, energy);
        canvas->SaveAs("output/energy_resolution.pdf");

        central_nine_energy_vec[run]->Draw("HIST e");
        fit = central_nine_energy_vec[run]->GetFunction("final_fit");
        fit->Draw("same");
        draw_text(fit, run_number, energy);
        canvas->SaveAs("output/energy_resolution.pdf");

        total_energy_vec[run]->Draw("HIST e");
        fit = total_energy_vec[run]->GetFunction("final_fit");
        fit->Draw("same");
        draw_text(fit, run_number, energy);
        canvas->SaveAs("output/energy_resolution.pdf");


        canvas->SetRightMargin(0.1);
        cog_distribution_vec[run]->Draw("COLZ");
        gPad->SetLogz();

        std::vector<TLine*> lines;
        for (int i = 1; i < 5; i++) {
            float loc = i - 0.5f;
            TLine* line_x = new TLine(loc, -0.5, loc, 4.5);
            line_x->SetLineStyle(2);
            line_x->Draw();
            TLine* line_y = new TLine(-0.5, loc, 4.5, loc);
            line_y->SetLineStyle(2);
            line_y->Draw();
        }
        TLine* line_x = new TLine(2, -0.5, 2, 4.5);
        line_x->SetLineStyle(2);
        line_x->SetLineColor(kRed);
        line_x->Draw();

        TLine* line_y = new TLine(-0.5, 2, 4.5, 2);
        line_y->SetLineStyle(2);
        line_y->SetLineColor(kRed);
        line_y->Draw();

        TEllipse* circle = new TEllipse(center_x, center_y, sigma_x, sigma_y);
        circle->SetLineColor(kRed);
        circle->SetLineWidth(2);
        circle->SetFillStyle(0);
        circle->Draw();

        canvas->SaveAs("output/energy_resolution.pdf");

        canvas->Clear();
        canvas->Divide(5, 5);
        for (int i = 0; i < 25; i++) {
            canvas->cd(i + 1);
            crystal_energy[run][crystal_mapping[i]]->Draw("HIST e");
            // gPad->SetLogx();
        }

        canvas->SaveAs("output/energy_resolution.pdf");
        

        std::vector<float> gev_calib;
        canvas->Clear();
        canvas->Divide(5, 5);
        float mean_calib = 0;
        for (int i = 0; i < 25; i++) {
            canvas->cd(i + 1);
            crystal_energy_shares[run][crystal_mapping[i]]->Draw("HIST e");
            crystal_energy_shares[run][crystal_mapping[i]]->GetXaxis()->SetRangeUser(0.001, 1);
            gPad->SetLogx();
            float signal_for_1gev = crystal_energy[run][crystal_mapping[i]]->GetMean() / crystal_energy_shares[run][crystal_mapping[i]]->GetMean();
            float x_coord = 0.4;
            if (i == 12) {
                x_coord = 0.15;
            }
            TLatex *text = new TLatex();
            text->SetNDC();
            text->SetTextSize(0.04);
            if (crystal_mapping[i] == 9) {
                continue;
            }
        }
        canvas->SaveAs("output/energy_resolution.pdf");
        mean_calib /= 24;   // Since we are excluding crystal 9

    }
    canvas->Clear();
    canvas->SaveAs("output/energy_resolution.pdf)");

    // float mean_x = cog_distribution->GetMean(1);
    // float mean_y = cog_distribution->GetMean(2);
    // float sigma_x = cog_distribution->GetStdDev(1);
    // float sigma_y = cog_distribution->GetStdDev(2);

    // std::cout << "Mean X: " << mean_x << " +/- " << sigma_x << std::endl;
    // std::cout << "Mean Y: " << mean_y << " +/- " << sigma_y << std::endl;

    // std::cout << "Total TOT events: " << tot_events << " out of " << nentries << std::endl;

}

