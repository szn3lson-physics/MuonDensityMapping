//====================================
// Kacper Dorszewski                 ||
// plot_rate_analysis_poisson.C      ||
//====================================

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <sstream>
#include <TLine.h>
#include <TLatex.h>
#include <vector>

const char* INPUT_FILE = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_rot_series/series_900.txt";
const char* OUTPUT_PLOT = "/home/kacper/MuonDensityMapping/figures/czestosc_zliczen_detektor_1.pdf";

void plot_rate_analysis_poisson() {
    // -----------------------------
    //   Apply global style
    // -----------------------------
    gStyle->SetOptStat(0);       // No statistics box
    gStyle->SetPadGridX(kTRUE);  // Grid (X)
    gStyle->SetPadGridY(kTRUE);  // Grid (Y)
    gStyle->SetPadTickX(1);      // Ticks on all sides
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kLightTemperature);

    int bin_width_deg = 5;
    int n_bins = 180 / bin_width_deg; 
    double x_min = 0.0;               
    double x_max = 180.0;

    // -----------------------------
    //   1. Create histogram and fill
    // -----------------------------
    TH1F *h_dat = new TH1F("h_dat", "Analiza czestosci zliczen mionow;Kat;Liczba czastek", n_bins, x_min, x_max);
    
    std::vector<double> sum_counts_dat(n_bins, 0.0);
    std::vector<double> sum_time_min(n_bins, 0.0);

    std::string line;

    // --- Wczytywanie danych detektora ---
    std::ifstream f_dat_txt(INPUT_FILE);
    if (f_dat_txt.is_open()) {
        while (std::getline(f_dat_txt, line)) {
            if (line.empty() || line.find("Angle") != std::string::npos) continue;
            
            std::stringstream ss(line);
            std::string item;
            int angle = 0; double time_min = 0.0; double counts = 0.0;

            if (std::getline(ss, item, ';')) angle = std::stoi(item);
            if (std::getline(ss, item, ';')) time_min = std::stod(item);
            if (std::getline(ss, item, ';')) counts = std::stod(item);

            int bin_idx = (angle - 1) / bin_width_deg;
            if (bin_idx >= 0 && bin_idx < n_bins) {
                sum_counts_dat[bin_idx] += counts;
                sum_time_min[bin_idx] += time_min;
            }
        }
        f_dat_txt.close();
    } else {
        std::cerr << "[BLAD] Nie mozna otworzyc pliku: " << INPUT_FILE << std::endl;
        return;
    }

    // --- Obliczanie wartosci srednich i wypelnianie histogramu ---
    double global_counts_dat = 0.0;
    double global_time_hours = 0.0;

    for (int i = 0; i < n_bins; ++i) {
        global_counts_dat += sum_counts_dat[i];
        if (sum_time_min[i] > 0) {
            double time_hours = sum_time_min[i] / 60.0;
            global_time_hours += time_hours;
            
            double rate_dat = sum_counts_dat[i] / time_hours;
            h_dat->SetBinContent(i + 1, rate_dat);
        }
    }

    // -----------------------------
    //   CALCULATION OF MEAN_Y AND SIGMA_Y
    // -----------------------------
    double mean_rate = 0.0;
    double global_sigma = 0.0;
    
    // Set manually the sigma of the histogram
    double manual_sigma = 2.0;
    
    if (global_time_hours > 0) {
        double avg_counts_per_bin = global_counts_dat / n_bins;
        double avg_time_h_per_bin = global_time_hours / n_bins;
        
        mean_rate = avg_counts_per_bin / avg_time_h_per_bin;
        global_sigma = TMath::Sqrt(avg_counts_per_bin) / avg_time_h_per_bin;
    }

    // Significance level
    double sigma = manual_sigma * global_sigma;

    // -----------------------------
    //   Histogram style
    // -----------------------------
    h_dat->SetLineColor(kRed);
    h_dat->SetLineWidth(2);
    h_dat->SetFillColorAlpha(kRed, 0.3); // Transparent fill
    h_dat->SetStats(0);

    // -----------------------------
    //   Canvas
    // -----------------------------
    TCanvas *c1 = new TCanvas("c1", "Histogram", 900, 600);

    // -----------------------------
    //   Draw (first to establish the axes)
    // -----------------------------
    h_dat->Draw("HIST");

    // Ensure the Y range is sufficient to display the lines
    double y_max_needed = TMath::Max(h_dat->GetMaximum() * 1.05, mean_rate + sigma * 1.1);
    h_dat->SetMaximum(y_max_needed);
    h_dat->SetMinimum(0); // Histograms usually start at 0

    // -----------------------------
    //   Mean Y line (Baseline)
    // -----------------------------
    // Draw a horizontal line at Mean_Y height
    TLine *mean_y_line = new TLine(x_min, mean_rate, x_max, mean_rate);
    mean_y_line->SetLineColor(kBlue);
    mean_y_line->SetLineWidth(3);
    mean_y_line->SetLineStyle(2); // Dotted line
    mean_y_line->Draw("SAME");
    
    // Text showing the mean Y value
    TLatex *mean_y_text = new TLatex();
    mean_y_text->SetTextFont(62); // 62 = bold w ROOT
    mean_y_text->SetTextSize(0.04);
    mean_y_text->SetTextColor(kBlue);
    mean_y_text->SetTextAlign(31); // Align Right, Bottom
    mean_y_text->DrawLatex(x_max * 0.95, mean_rate * 1.05, Form("Srednia (Y): %.2f", mean_rate));


    // -----------------------------
    //   Significance lines (+- sigma Y)
    // -----------------------------
    
    // Upper significance level
    double y_high = mean_rate + sigma;
    TLine *sig_high = new TLine(x_min, y_high, x_max, y_high);
    sig_high->SetLineColor(kBlack);
    sig_high->SetLineWidth(2);
    sig_high->SetLineStyle(3); // Dashed line
    sig_high->Draw("SAME");

    // Lower significance level
    double y_low = mean_rate - sigma;
    // Limit the bottom line to 0, as counts cannot be negative
    if (y_low < 0) y_low = 0; 
    TLine *sig_low = new TLine(x_min, y_low, x_max, y_low);
    sig_low->SetLineColor(kBlack);
    sig_low->SetLineWidth(2);
    sig_low->SetLineStyle(3); // Dashed line
    sig_low->Draw("SAME");

    // Text for the sigma range
    TLatex *sigma_text = new TLatex();
    sigma_text->SetTextFont(62); // 62 = bold w ROOT
    sigma_text->SetTextSize(0.035);
    sigma_text->SetTextColor(kBlack);
    sigma_text->SetTextAlign(12); // Align Left, Center
    
    // +sigma label
    sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01, y_high * 0.98, Form("+%.1f#sigma", manual_sigma));
    // -sigma label
    if (y_low > 0) {
        sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01, y_low * 0.9, Form("-%.1f#sigma", manual_sigma));
    }

    c1->SaveAs(OUTPUT_PLOT);
    std::cout << "[SUKCES] Wygenerowano wykres: " << OUTPUT_PLOT << std::endl;
}