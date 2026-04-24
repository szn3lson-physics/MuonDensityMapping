#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TMath.h>

void His_1D_Rates(
    const std::string& input_file,
    const std::string& output_pdf,
    const std::string& title_text,
    int bin_width_deg,
    double sigma_input
) {
    // -----------------------------
    //   1. Data accumulation
    // -----------------------------
    int n_bins = 180 / bin_width_deg;

    std::vector<double> sum_counts(n_bins, 0.0);
    std::vector<double> sum_time_min(n_bins, 0.0);

    std::fstream file(input_file, std::ios::in);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << input_file << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Skip the header: Angle;Time_Minutes;Counts_34f

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string item;
        int angle = 0;
        double time_min = 0.0, counts = 0.0;

        if (std::getline(ss, item, ';')) angle = std::stoi(item);
        if (std::getline(ss, item, ';')) time_min = std::stod(item);
        if (std::getline(ss, item, ';')) counts = std::stod(item);

        // Convert 1-180 angle to 0-based bin index
        int bin_idx = (angle - 1) / bin_width_deg;
        if (bin_idx >= 0 && bin_idx < n_bins) {
            sum_counts[bin_idx] += counts;
            sum_time_min[bin_idx] += time_min;
        }
    }
    file.close();

    // -----------------------------
    //   2. Create histogram and fill
    // -----------------------------
    TH1F *hist = new TH1F(
        "hist",
        Form("%s (Bin: %d#circ);Zenith Angle [#circ];Muon Rate [Counts / Hour]", title_text.c_str(), bin_width_deg),
        n_bins, 0, 180
    );

    double total_rate = 0.0;
    double global_counts = 0.0;
    double global_time_hours = 0.0;

    for (int i = 0; i < n_bins; ++i) {
        double rate = 0.0;
        if (sum_time_min[i] > 0) {
            double time_hours = sum_time_min[i] / 60.0;
            rate = sum_counts[i] / time_hours;
            
            global_counts += sum_counts[i];
            global_time_hours += time_hours;
        }
        hist->SetBinContent(i + 1, rate); // ROOT bins are 1-based
        total_rate += rate;
    }

    // -----------------------------
    //   3. CALCULATION OF MEAN_Y AND SIGMA_Y
    // -----------------------------
    double mean_y = total_rate / n_bins;

    // Statistically correct Poisson error for a Rate
    double avg_counts = global_counts / n_bins;
    double avg_time_h = global_time_hours / n_bins;
    
    // 1 Sigma = sqrt(N) / T
    double sigma_1 = 0.0;
    if (avg_time_h > 0) {
        sigma_1 = TMath::Sqrt(avg_counts) / avg_time_h;
    }

    double sigma = sigma_input * sigma_1;

    // -----------------------------
    //   4. Apply global style
    // -----------------------------
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kLightTemperature);

    // -----------------------------
    //   5. Histogram style
    // -----------------------------
    hist->SetLineColor(kRed);
    hist->SetLineWidth(2);
    hist->SetFillColorAlpha(kRed, 0.3);
    hist->SetStats(0);

    // -----------------------------
    //   6. Canvas and Drawing
    // -----------------------------
    TCanvas *c1 = new TCanvas("c1", "Histogram", 900, 600);

    // Dynamic Y-axis scaling to fit the +Sigma line comfortably
    double y_max_needed = TMath::Max(hist->GetMaximum() * 1.15, mean_y + sigma * 1.2);
    hist->SetMaximum(y_max_needed);
    hist->SetMinimum(0);

    hist->Draw("HIST");

    double x_min = hist->GetXaxis()->GetXmin();
    double x_max = hist->GetXaxis()->GetXmax();

    // -----------------------------
    //   7. Mean Y line
    // -----------------------------
    TLine *mean_y_line = new TLine(x_min, mean_y, x_max, mean_y);
    mean_y_line->SetLineColor(kBlue);
    mean_y_line->SetLineWidth(3);
    mean_y_line->SetLineStyle(2);
    mean_y_line->Draw("SAME");

    TLatex *mean_y_text = new TLatex();
    mean_y_text->SetTextSize(0.04);
    mean_y_text->SetTextColor(kBlue);
    mean_y_text->SetTextAlign(31);
    mean_y_text->DrawLatex(x_max * 0.95, mean_y * 1.05,
                           Form("Mean Rate: %.1f / h", mean_y));

    // -----------------------------
    //   8. Significance lines
    // -----------------------------
    double y_high = mean_y + sigma;
    double y_low  = TMath::Max(0.0, mean_y - sigma);

    TLine *sig_high = new TLine(x_min, y_high, x_max, y_high);
    sig_high->SetLineColor(kBlack);
    sig_high->SetLineWidth(2);
    sig_high->SetLineStyle(3);
    sig_high->Draw("SAME");

    TLine *sig_low = new TLine(x_min, y_low, x_max, y_low);
    sig_low->SetLineColor(kBlack);
    sig_low->SetLineWidth(2);
    sig_low->SetLineStyle(3);
    sig_low->Draw("SAME");

    TLatex *sigma_text = new TLatex();
    sigma_text->SetTextSize(0.035);
    sigma_text->SetTextColor(kBlack);
    sigma_text->SetTextAlign(12);

    sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01,
                          y_high * 1.02,
                          Form("+%.1f#sigma", sigma_input));

    if (y_low > 0) {
        sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01,
                              y_low * 0.92,
                              Form("-%.1f#sigma", sigma_input));
    }

    // -----------------------------
    //   9. SAVE PDF
    // -----------------------------
    std::string dir = output_pdf.substr(0, output_pdf.find_last_of('/'));
    gSystem->mkdir(dir.c_str(), kTRUE);
    c1->SaveAs(output_pdf.c_str());

    delete c1;
    delete hist;
}