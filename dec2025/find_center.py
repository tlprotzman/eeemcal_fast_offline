import ROOT
from array import array

# Vertical
# runs = [68, 69, 70, 71/, 72, 75, 76, 78]
# position = [0, 2, 4, -2, -4, 6, 2, 0]
# runs =     [96, 97, 98,  99, 101, 102]
# position = [ 0, -5,  5, -10, -12, -15]

# Horizontal
# runs =     [83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 96]
# position = [0,  -4, -2,  2,  4,  6,  8, 10, 12, 14,  4, 4.7]

# New vertical
runs = [234, 235, 236, 237, 239]
position = [0, 3, -3, -6, 6]

# New horizontal
runs = [234, 241, 242, 243, 244, 245, 246]
position = [0, -3, 3, 6, 9, 1.7, 4.1]

# Corner
runs = [248, 249, 250, 255]
position = [42.5, 44.7, 40.4, 38.5]

runs = [251, 252, 253, 254]
position = [-44.5, -40.5, -38.6, -46.4]

# New new (sigh) vertical
runs = [303, 304, 305, 306, 307]
position = [0, -3.2, -5.1, 3.1, 5]

runs = [303, 308, 309, 310, 311]
position = [0, 3.2, 5.2, -2.9, -5.2]

# new corner vert
runs = [355, 356, 357, 358, 359]
position = [-42.5, -44.5, -46.5, -40.3, -38.3]

# new corner horiz
runs = [360, 361, 362, 363, 364]
positions = [39.5, 37.5, 44.5, 47.3, -42.5]

# crystal 20
runs = [365, 366, 367, 368, 369]
position = [-42.5, -39.5, -37.5, -45.4, -47.5]

peak = []
peak_err = []

for i, run in enumerate(runs):
    file = ROOT.TFile(f"/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run{run:03d}.root")
    hist_diff = ROOT.TH1F("{run}_adc_pedestal_diff", f"Run{run} ADC - Pedestal", 1000, 0, 10000)
    tree = file.Get("events")
    eeemcal_16i_channel_a_map = [2,  6, 11, 15,  0,  4,  9, 13,
                                1,  5, 10, 14,  3,  7, 12, 16]

    eeemcal_16i_channel_b_map = [20, 24, 29, 33, 18, 22, 27, 31,
                                19, 23, 28, 32, 21, 25, 30, 34]

    eeemcal_16i_channel_c_map = [67, 63, 59, 55, 69, 65, 61, 57,
                                70, 66, 60, 56, 68, 64, 58, 54]

    eeemcal_16i_channel_d_map = [50, 46, 40, 36, 52, 48, 42, 38,
                                51, 47, 43, 39, 49, 45, 41, 37]
    # center_channels = [2,  6, 11, 15,  0,  4,  9, 13,
                                    #   1,  5, 10, 14,  3,  7, 12, 16]
  
    # center_channels = [ch + (144 * 2) + (72 * 1) for ch in center_channels]
    center_channels = [ch + (144 * 1) + (72 * 1) for ch in eeemcal_16i_channel_c_map]
    for event in tree:
        tot = 0
        for ch in center_channels:
            tot += event.hit_max[ch] - event.hit_pedestal[ch]
        hist_diff.Fill(tot)

    fit_result = hist_diff.Fit("gaus", "", "", 6000, 8500)
    fit_func = hist_diff.GetFunction("gaus")
    mean = fit_func.GetParameter(1)
    mean_err = fit_func.GetParError(1)

    canvas = ROOT.TCanvas()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)

    fit_func = hist_diff.GetFunction("gaus")
    mean = fit_func.GetParameter(1)
    std_dev = fit_func.GetParameter(2)
    text = f"Mean: {mean:.2f}, Std Dev: {std_dev:.2f}"


    # Display histogram
    # canvas.SetLogy()
    hist_diff.Draw()
    latex.DrawLatex(0.15, 0.8, text)
    text = f"resolution: {std_dev/mean*100:.2f} %"
    latex.DrawLatex(0.15, 0.7, text)

    canvas.SaveAs(f"output/{run}.png")

    
    peak.append(mean)
    peak_err.append(mean_err)
    file.Close()
    # break

print(peak)
print(peak_err)

canvas = ROOT.TCanvas()
graph = ROOT.TGraphErrors(len(runs),
                          array('d', position),
                          array('d', peak),
                          array('d', [0]*len(runs)),
                          array('d', peak_err))
graph.SetTitle("Center Crystal Response vs Position; Position (mm); ADC Peak Value")
graph.SetMarkerStyle(21)
graph.Draw("AP")

# Fit the graph with a Gaussian
fit_func = ROOT.TF1("fit_func", "gaus", min(position)-1, max(position)+1)
fit_func.SetParameters(8000, 6, 20.7785)
graph.Fit(fit_func, "")  # "Q" for quiet mode

# Optionally, draw the fit function on the canvas
fit_func.SetLineColor(ROOT.kRed)
fit_func.Draw("same")

# Write the fit mean and uncertainty on the canvas
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
fit_mean = fit_func.GetParameter(1)
fit_mean_err = fit_func.GetParError(1)
text = f"Fit Mean: {fit_mean:.2f} #pm {fit_mean_err:.2f} mm"
latex.DrawLatex(0.15, 0.8, text)

# Fit with polynomial degree 2
fit_func_pol2 = ROOT.TF1("fit_func_pol2", "pol2", min(position)-1, max(position)+1)
graph.Fit(fit_func_pol2, "")
fit_func_pol2.SetLineColor(ROOT.kBlue)
fit_func_pol2.Draw("same")

# Find peak location for pol2 (vertex of parabola: -b/2a)
a = fit_func_pol2.GetParameter(2)
b = fit_func_pol2.GetParameter(1)
pol2_peak = -b / (2 * a)
text = f"Pol2 Peak: {pol2_peak:.2f} mm"
latex.DrawLatex(0.15, 0.75, text)
print(f"Gaussian peak location: {fit_mean:.2f} mm")
print(f"Pol2 peak location: {pol2_peak:.2f} mm")

canvas.SaveAs("output/center_crystal_response.png")