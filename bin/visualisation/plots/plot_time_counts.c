#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLine.h>
#include <TLatex.h>
#include <TSystem.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

void plot_rate_1D(
    const std::string& input_file,
    const std::string& output_pdf,
    const std::string& title_text,
    int bin_width_deg,
    double sigma_input
) {
    int n_bins = 180 / bin_width_deg;

    // -----------------------------
    //   1. Create histograms
    // -----------------------------
    // Create working histograms for Counts and Time
    TH1F *hist_counts = new TH1F("h_counts", "", n_bins, 0, 180);
    TH1F *hist_time   = new TH1F("h_time", "", n_bins, 0, 180);
    
    // Main histogram for the Result (Counts / Minute)
    TH1F *hist_rate = new TH1F(
        "hist_rate",
        Form("%s (Bin: %d#circ);Angle [#circ];Normalized Rate [Counts / Minute]", title_text.c_str(), bin_width_deg),
        n_bins, 0, 180
    );

    // -----------------------------
    //   2. Loading from CSV file
    // -----------------------------
    std::ifstream file(input_file);
    if (!file.is_open()) {
        std::cerr << "[ERROR] Cannot open: " << input_file << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Skip header (Angle;Time_Minutes;Counts_34f)

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string item;
        
        std::getline(ss, item, ';'); int angle = std::stoi(item);
        std::getline(ss, item, ';'); double time = std::stod(item);
        std::getline(ss, item, ';'); double counts = std::stod(item);

        // Bin center for a given angle (so ROOT places it correctly)
        double bin_center = angle + 0.5;

        // Fill with values (weight = time or counts)
        hist_counts->Fill(bin_center, counts);
        hist_time->Fill(bin_center, time);
    }
    file.close();

    // DIVISION: Rate = Counts / Time
    hist_rate->Divide(hist_counts, hist_time);

    // -----------------------------
    //   3. CALCULATION OF MEAN_Y and SIGMA_Y
    // -----------------------------
    double manual_sigma = sigma_input;

    // Calculate mean only from non-empty bins
    double sum_of_content = 0.0;
    int valid_bins = 0;
    for (int i = 1; i <= n_bins; ++i) {
        double val = hist_rate->GetBinContent(i);
        if (hist_time->GetBinContent(i) > 0) { // If the detector was looking there at all
            sum_of_content += val;
            valid_bins++;
        }
    }
    
    double mean_y = (valid_bins > 0) ? (sum_of_content / valid_bins) : 0.0;
    
    // According to your original code: sigma = N * sqrt(mean)
    double sigma = manual_sigma * TMath::Sqrt(mean_y);

    // -----------------------------
    //   4. Apply global style
    // -----------------------------
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // -----------------------------
    //   5. Histogram style
    // -----------------------------
    hist_rate->SetLineColor(kRed);
    hist_rate->SetLineWidth(2);
    hist_rate->SetFillColorAlpha(kRed, 0.3);

    // -----------------------------
    //   6. Canvas
    // -----------------------------
    TCanvas *c1 = new TCanvas("c1", "Histogram", 900, 600);
    c1->SetBottomMargin(0.12);
    c1->SetLeftMargin(0.12);

    hist_rate->GetYaxis()->SetTitleOffset(1.2);
    hist_rate->Draw("HIST");

    double x_min = hist_rate->GetXaxis()->GetXmin();
    double x_max = hist_rate->GetXaxis()->GetXmax();

    double y_max_needed = TMath::Max(hist_rate->GetMaximum() * 1.05, mean_y + sigma * 1.1);
    if (y_max_needed <= 0) y_max_needed = 1.0; // Safety margin
    
    hist_rate->SetMaximum(y_max_needed);
    hist_rate->SetMinimum(0);

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
    mean_y_text->DrawLatex(x_max * 0.98, mean_y * 1.05,
                           Form("Mean (Y): %.3f", mean_y));

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
                          Form("+%.1f#sigma", manual_sigma));

    if (y_low > 0) {
        sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01,
                              y_low * 1.02,
                              Form("-%.1f#sigma", manual_sigma));
    }

    // -----------------------------
    //   9. SAVE PDF
    // -----------------------------
    std::string dir = output_pdf.substr(0, output_pdf.find_last_of('/'));
    gSystem->mkdir(dir.c_str(), kTRUE); // Automatically creates folders if they don't exist
    c1->SaveAs(output_pdf.c_str());

    // Memory cleanup
    delete c1;
    delete hist_rate;
    delete hist_counts;
    delete hist_time;
}

// =============================================================================
// MAIN FUNCTION TO RUN
// =============================================================================
void plot_time_counts() {
    // Which datasets we are grouping (from 1 to 10)
    std::vector<int> datasets = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; 
    // Bin widths
    std::vector<int> bin_widths = {3, 6, 9, 12, 15, 18};

    for (int x : datasets) {
        for (int bw : bin_widths) {
            
            // Input file (generated by the previous C++ script)
            std::string input = "D:/MuonDensityMapping/output/all/" + std::to_string(x) + "_zussamen/time_all.txt";

            // Where to save the generated PDF
            std::string output = "D:/MuonDensityMapping/results/all/" + std::to_string(x) + "_zussamen/time_analysis/hist_rate_" + std::to_string(bw) + "deg.pdf";
            
            // Call with a hardcoded sigma of 3 (like in your code: bw, 3)
            plot_rate_1D(
                input, 
                output, 
                "Normalized Rate (34f) / Data " + std::to_string(x), 
                bw, 
                3.0
            );
        }
    }
    std::cout << "\nAll plots generated successfully!" << std::endl;
}