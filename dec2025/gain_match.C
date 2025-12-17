#include <map>
#include <vector>


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

// float calculate_signal(int *adc_values) {
//     // Pedestal is the mean of the first three samples
//     float pedestal = (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;
//     // Signal is the sum of samples 4 through 8 minus pedestals
//     float signal = 0.0f;
//     for (int i = 3; i < 8; ++i) {
//         signal += (adc_values[i] - pedestal);
//     }
//     return signal;
// }

float calculate_signal(int *adc_values) {
    // Pedestal is the mean of the first three samples
    float pedestal = (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;
    
    // Signal is the sum of the three greatest values
    std::vector<float> samples;
    for (int i = 3; i < 20; ++i) {
        samples.push_back(adc_values[i] - pedestal);
    }
    std::sort(samples.begin(), samples.end(), std::greater<float>());
    float signal = samples[0] + samples[1] + samples[2] + samples[3];
    return signal;
}

void gain_match() {
    gStyle->SetOptStat(0);
    std::vector<int> run_numbers = {383, 385, /*386, 387,*/ 388, 389, 390, 391, 392, 393, 394, /*395,*/ 396, 397, 398, 399, /*400,*/ 401, 402, 403};//, 404, 405, 406};
    std::vector<int> run_crystal = { 14,  19, /* 24,  23,*/  18,  13,  14,   8,   3,   4,   9, /*  2,*/   7,  12,  17,  22, /* 21,*/  16,  11,   6};//,   1,   0,   5};
    std::vector<TH1*> histograms(400); // 25 crystals * 16 SiPMs

    // Set up histograms
    auto mapping = read_mapping("eeemcal_desy_dec2025_mapping.csv");
    for (int crystal = 0; crystal < 25; crystal++) {
        for (int sipm = 0; sipm < 16; sipm++) {
            int channel = mapping[crystal][sipm];
            auto hist = new TH1F(Form("crystal_%d_sipm_%d", crystal, sipm),
                                 Form("Crystal %d SiPM %d Signal;Signal (ADC);Counts", crystal, sipm),
                                 200, 0, 2000);
            histograms[crystal * 16 + sipm] = hist;
        }
    }


    // Process each run
    for (int i = 0; i < run_numbers.size(); i++) {
        int run_number = run_numbers[i];
        int crystal = run_crystal[i];
        char file_path[256];
        sprintf(file_path, "/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number);
        TFile* root_file = TFile::Open(file_path);
        TTree* tree = (TTree*)root_file->Get("events");
        int adc[576][20];
        tree->SetBranchAddress("adc", &adc);
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("adc", 1);

        // Fill histograms from tree
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            for (int sipm = 0; sipm < 16; sipm++) {
                int channel = mapping[crystal][sipm];
                float signal = calculate_signal(adc[channel]);
                histograms[crystal * 16 + sipm]->Fill(signal);
            }
        }
    }

    std::vector<float> peak_locations(400);

    // Fit each histogram to find the peak
    for (int channel = 0; channel < histograms.size(); channel++) {
        TH1* hist = histograms[channel];
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

        TF1* final_fit = new TF1("final_fit", "gaus", peak - sigma, peak + sigma);
        hist->Fit(final_fit, "RQ");

        peak_locations[channel] = final_fit->GetParameter(1);
    }

    // Print the average peak location
    float total_peak = 0.0f;
    int count = 0;
    for (float peak : peak_locations) {
        if (peak > 0) {
            total_peak += peak;
            count++;
        }
    }
    float average_peak = total_peak / count;
    std::cout << "Average Peak Location: " << average_peak << std::endl;

    // Calculate the gain factors
    float target = 1300;//average_peak;
    std::vector<float> gain_factors(400);
    for (int i = 0; i < peak_locations.size(); i++) {
        if (peak_locations[i] > 0) {
            gain_factors[i] = target / peak_locations[i];
        } else {
            gain_factors[i] = 1.0f; // Default factor if no peak
        }
    }

    // Save all histograms to a single PDF
    char output_file[256] = "output/gain_matching.pdf";
    TCanvas* canvas = new TCanvas("gain_matching", "", 800, 600);
    for (int crystal = 0; crystal < run_crystal.size(); crystal++) {
        canvas->Clear();
        canvas->Divide(4, 4);
        for (int sipm = 0; sipm < 16; sipm++) {
            canvas->cd(sipm + 1);
            auto hist = histograms[run_crystal[crystal] * 16 + sipm];
            hist->Draw();
        }
        if (crystal == 0) {
            canvas->SaveAs(Form("%s(", output_file));
        } else {
            canvas->SaveAs(output_file);
        }
    }

    canvas->Clear();


    // Make a histogram of gain factors
    TH1F* gain_hist = new TH1F("gain_factors", "Gain Factors;Channel;Gain Factor", 400, 0, 400);
    for (int channel = 0; channel < gain_factors.size(); channel++) {
        gain_hist->SetBinContent(channel + 1, gain_factors[channel]);
        gain_hist->SetBinError(channel + 1, 0.0001);
        std::cout << "Channel " << channel << ": Gain Factor = " << gain_factors[channel] << std::endl;
    }
    gain_hist->SetMinimum(0);
    gain_hist->SetMaximum(2);
    gain_hist->Draw("e");
    canvas->SaveAs(Form("%s)", output_file));
    TFile* out_root = new TFile("output/gain_factors.root", "RECREATE");
    gain_hist->Write();
    out_root->Close();
}