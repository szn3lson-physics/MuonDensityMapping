#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLine.h>

void plot_rate_per_minute() {

    // ================= PARAMETERS =================
    const int START = 1;
    const int STOP  = 10; // Set to 10 for all files
    const double bin_width_minutes = 1.0; // Bin width = 1 minute
    // ==============================================

    std::vector<double> global_times;
    std::vector<double> boundaries; // Cut points where red lines will be drawn
    
    double global_time_offset = 0.0;

    std::cout << "Starting to load and merge files..." << std::endl;

    for (int i = START; i <= STOP; ++i) {

        // Dynamic file path
        std::string filename = "D:\\MuonDensityMapping\\data\\dane_" + std::to_string(i) + ".log";
        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "[WARNING] Cannot open file: " << filename << ". Skipping..." << std::endl;
            continue;
        }

        std::string line;
        bool is_first_line = true;
        double local_t0 = 0.0;
        double max_t_this_file = 0.0;
        int count_this_file = 0;

        while (std::getline(file, line)) {

            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string token;
            int col = 0;
            double t_raw = -1.0;

            // Parse semicolon-separated values (column index 1 contains time)
            while (std::getline(ss, token, ';')) {
                if (col == 1) {
                    try {
                        t_raw = std::stod(token);
                    } catch (...) {
                        // Ignore text headers
                    }
                    break;
                }
                col++;
            }

            // If time parsed correctly
            if (t_raw >= 0) {

                if (is_first_line) {
                    local_t0 = t_raw; // Reset clock for this file
                    is_first_line = false;
                }

                // Time in minutes for this file (from 0 to its end)
                double t_min = (t_raw - local_t0) / 60000.0;

                // Global time (shifted by previous files' durations)
                double t_global = global_time_offset + t_min;
                global_times.push_back(t_global);
                count_this_file++;

                if (t_min > max_t_this_file) {
                    max_t_this_file = t_min;
                }
            }
        }

        file.close();

        std::cout << "File " << i << ": Loaded " << count_this_file
                  << " counts. Series duration: "
                  << max_t_this_file << " min." << std::endl;

        // Update global offset so next file starts where this one ended
        global_time_offset += max_t_this_file;

        // Store series boundary (to draw red separator line)
        // Do not draw after last file
        if (i < STOP && count_this_file > 0) {
            boundaries.push_back(global_time_offset);
        }
    }

    if (global_times.empty()) {
        std::cerr << "[ERROR] No valid data to plot!" << std::endl;
        return;
    }

    double total_time = global_time_offset;
    int n_bins = std::max(1, (int)std::ceil(total_time / bin_width_minutes));

    // Create main histogram
    TH1F *hist = new TH1F(
        "hist",
        "Counting Frequency Over Time (Series 1-10);Total Experiment Time [minutes];Counts per minute",
        n_bins, 0, total_time
    );

    for (double t : global_times) {
        hist->Fill(t);
    }

    // ================= PLOT STYLING =================
    gStyle->SetOptStat(0); 
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(2);
    hist->SetFillStyle(0);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->SetMinimum(0);

    // Determine Y max for full-height separator lines
    double y_max = hist->GetMaximum() * 1.05;
    hist->SetMaximum(y_max);

    TCanvas *c1 = new TCanvas("c1", "Merged Series", 1200, 600);
    c1->SetLeftMargin(0.10);
    c1->SetBottomMargin(0.12);
    c1->SetRightMargin(0.05);

    hist->Draw("HIST");

    // ================= DRAW RED SEPARATOR LINES =================
    for (double b_time : boundaries) {
        TLine *line = new TLine(b_time, 0, b_time, y_max);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->SetLineStyle(2); // 2 = dashed line
        line->Draw("SAME");
    }

    c1->SaveAs("D:\\MuonDensityMapping\\figures\\frequency_over_time.png");

    std::cout << "\nPlot successfully generated." << std::endl;
}
