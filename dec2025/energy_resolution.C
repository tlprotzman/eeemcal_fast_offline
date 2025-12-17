#include <map>
#include <vector>

float sigma_cut = 2;

float center_x = 1.97505;
float sigma_x = 0.192906 * sigma_cut;
float center_y = 1.97451;
float sigma_y = 0.196824 * sigma_cut;
float radius = 0.1;

int run_number = 315;
float beam_energy = 1;

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

// float calculate_signal(int *adc_values, float gain) {
//     // Pedestal is the mean of the first three samples
//     float pedestal = (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;
//     // Signal is the sum of samples 4 through 8 minus pedestals
//     float signal = 0.0f;
//     for (int i = 3; i < 8; ++i) {
//         signal += (adc_values[i] - pedestal);
//     }
//     return signal * gain;
// }

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

float calculate_signal(uint32_t *adc_values, uint32_t *tot_values, float gain) {
    // Start with the ADC portion of the signal
    float signal = calculate_signal(adc_values, gain);
    // Check if there is a ToT value
    uint32_t tot = 0;
    for (int i = 0; i < 20; ++i) {
        if (tot_values[i] > 0) {
            tot = tot_values[i];
            break;
        }
    }
    // If there is a ToT value, then we discard the ADC and convert the ToT to energy
    if (tot == 0) {
        return signal;
    }
    float adc_equivalent = 0.0f;
    return 0;

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

void draw_text(TF1* fit) {
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
    text->DrawLatex(0.93, 0.75, Form("Peak: %.0f", peak));
    text->DrawLatex(0.93, 0.70, Form("Sigma: %.0f", sigma));
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

void energy_resolution() {
    gStyle->SetOptStat(0);
    float energy = 1;

    int central_crystal_index = 12;
    int center_nine_indexes[8] = {7, 8, 9, 11, 13, 17, 18, 19};
    int remaining_indexes[16] = {0, 1, 2, 3, 4, 5, 6, 10, 11, 15, 16, 20, 21, 22, 23, 24};


    TH1* central_crystal_energy = new TH1F("central_crystal_energy", "Central Crystal Energy;Energy (ADC);Events", 500, 0, 75000);
    TH1* central_nine_energy    = new TH1F("central_nine_energy", "Central 3x3 Energy;Energy (ADC);Events", 500, 0, 75000);
    TH1* total_energy           = new TH1F("total_energy", "Total Energy;Energy (ADC);Events", 500, 0, 75000);
    TH2* cog_distribution       = new TH2F("cog_distribution", "Center of Gravity Distribution;X (# Crystals);Y (# Crystals)", 100, -0.5, 4.5, 100, -0.5, 4.5);

    std::vector<std::vector<TH1*>> sipm_energy;
    std::vector<TH1*> crystal_energy;
    std::vector<TH1*> crystal_energy_shares;
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 16; j++) {
            sipm_energy.push_back(std::vector<TH1*>());
            sipm_energy[i].push_back(new TH1F(Form("crystal_%02d_sipm_%02d_energy", i, j), Form("Crystal %02d SiPM %02d Energy;Energy (ADC);Events", i, j), 500, 0, 4000));
        }   
        crystal_energy.push_back(new TH1F(Form("crystal_%02d_energy", i), Form("Crystal %02d Energy;Energy (ADC);Events", i), 500, 0, 40000));
        crystal_energy_shares.push_back(new TH1F(Form("crystal_%02d_energy_share", i), Form("Crystal %02d Energy Share;Energy Share;Events", i), 100, 0, 1));
    }

    TH2* pedestals = new TH2F("pedestals", "Pedestals;Channel;Pedestal (ADC)", 72*8, 0, 72*8, 1024, 0, 1024);

    auto mapping = read_mapping("eeemcal_desy_dec2025_mapping.csv");

    TH1* gain_factor = nullptr;
    TFile* gain_file = TFile::Open("output/gain_factors.root");
    if (gain_file && !gain_file->IsZombie()) {
        gain_factor = (TH1*)gain_file->Get("gain_factors");
        std::cout << "Loaded gain factors from file." << std::endl;
    }
    if (!gain_factor) {
        gain_factor = new TH1F("gain_factors", "Gain Factor", 400, 0, 400);
        for (int i = 1; i <= 400; ++i) {
            gain_factor->SetBinContent(i, 1.0);
        }
    }


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

    Long64_t nentries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nentries; ++entry) {
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
            float crystal_signal = 0.0f;
            for (int sipm = 0; sipm < 16; sipm++) {
                int channel = mapping[crystal][sipm];
                float gain = gain_factor->GetBinContent(crystal * 16 + sipm + 1);
                float channel_signal = calculate_signal(adc[channel], gain);
                sipm_energy[crystal][sipm]->Fill(channel_signal);
                pedestals->Fill(channel, (adc[channel][0] + adc[channel][1] + adc[channel][2]) / 3.0f);
                crystal_signal += channel_signal;
                if (!is_tot_event && is_tot(tot[channel])) {
                    is_tot_event = true;
                }
            }
            signals[crystal] = crystal_signal;
        }
        if (is_tot_event) {
            tot_events++;
            // continue;
        }

        float x_cog, y_cog;
        bool keep = calculate_cog(cog_distribution, signals);
        keep &= (!is_tot_event);
        if (!keep) {
            continue;
        }


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
        central_nine_signal += central_signal;
        total_signal += central_nine_signal;
        central_crystal_energy->Fill(central_signal);
        central_nine_energy->Fill(central_nine_signal);
        total_energy->Fill(total_signal);

        // Energy share fill
        for (int i = 0; i < 25; i++) {
            crystal_energy[i]->Fill(signals[i]);
            crystal_energy_shares[i]->Fill(signals[i] / total_signal);
        }

    }

    fit_peak(central_crystal_energy);
    fit_peak(central_nine_energy);
    fit_peak(total_energy);


    int crystal_mapping[25] = {
        4, 9, 14, 19, 24,
        3, 8, 13, 18, 23,
        2, 7, 12, 17, 22,
        1, 6, 11, 16, 21,
        0, 5, 10, 15, 20
    };

    TCanvas* canvas = new TCanvas("gain_matching", "", 800, 600);
    canvas->SetRightMargin(0.05);
    central_crystal_energy->Draw("HIST e");
    auto fit = central_crystal_energy->GetFunction("final_fit");
    fit->Draw("same");
    draw_text(fit);
    canvas->SaveAs("output/energy_resolution.pdf(");

    central_nine_energy->Draw("HIST e");
    fit = central_nine_energy->GetFunction("final_fit");
    fit->Draw("same");
    draw_text(fit);
    canvas->SaveAs("output/energy_resolution.pdf");

    total_energy->Draw("HIST e");
    fit = total_energy->GetFunction("final_fit");
    fit->Draw("same");
    draw_text(fit);
    canvas->SaveAs("output/energy_resolution.pdf");


    canvas->SetRightMargin(0.1);
    cog_distribution->Draw("COLZ");
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
    pedestals->Draw("COLZ");
    canvas->SaveAs("output/energy_resolution.pdf");

    canvas->Clear();
    canvas->Divide(5, 5);
    for (int i = 0; i < 25; i++) {
        canvas->cd(i + 1);
        crystal_energy[crystal_mapping[i]]->Draw("HIST e");
    }

    canvas->SaveAs("output/energy_resolution.pdf");
    
    canvas->Clear();
    canvas->Divide(5, 5);
    for (int i = 0; i < 25; i++) {
        canvas->cd(i + 1);
        crystal_energy_shares[crystal_mapping[i]]->Draw("HIST e");
    }
    canvas->SaveAs("output/energy_resolution.pdf");




    for (int i = 0; i < 25; i++) {
        canvas->Clear();
        canvas->Divide(4, 4);
        for (int j = 0; j < 16; j++) {
            canvas->cd(j + 1);
            sipm_energy[i][j]->Draw("HIST e");
        }
        canvas->SaveAs("output/energy_resolution.pdf");
    }
    
    


    canvas->SaveAs("output/energy_resolution.pdf)");

    float mean_x = cog_distribution->GetMean(1);
    float mean_y = cog_distribution->GetMean(2);
    float sigma_x = cog_distribution->GetStdDev(1);
    float sigma_y = cog_distribution->GetStdDev(2);

    std::cout << "Mean X: " << mean_x << " +/- " << sigma_x << std::endl;
    std::cout << "Mean Y: " << mean_y << " +/- " << sigma_y << std::endl;

    std::cout << "Total TOT events: " << tot_events << " out of " << nentries << std::endl;

}

