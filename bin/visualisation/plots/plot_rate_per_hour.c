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
#include <TLegend.h>
#include <TLatex.h>

void plot_rate_per_hour() {
    // ================= PARAMETERS =================
    const int START = 1;
    const int END   = 12; // Scans files from dane_1.log to dane_12.log
    
    const char* INPUT_FILE ="/home/kacper/MuonDensityMapping/bin/data_analysis/data/dane_";
    std::string path = "/home/kacper/MuonDensityMapping/figures/rate_per_hour/";

    // Bin width in hours. 
    const double bin_width_hours = 1.0; 
    // ==============================================

    std::vector<double> t_all, t_4, t_3, t_3f, t_34, t_34f;
    std::vector<double> boundaries;      // Czerwone linie (granice plików)
    std::vector<double> internal_resets; // Zielone linie (restarty mikrokontrolera)
    std::vector<double> file_centers; 
    std::vector<int> file_numbers;    
    
    double global_time_offset = 0.0;

    std::cout << "Starting data extraction and coincidence segregation..." << std::endl;

    for (int i = START; i <= END; ++i) {
        std::string filename = INPUT_FILE + std::to_string(i) + ".log";
        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "[WARNING] Cannot open file: " << filename << ". Skipping..." << std::endl;
            continue;
        }

        std::string line;
        bool is_first_line = true;
        
        double local_t0 = 0.0;
        double prev_t_raw = 0.0;
        double local_segment_offset = 0.0;
        
        double max_t_this_file = 0.0;
        int count_this_file = 0;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            if (line.rfind("Adafruit", 0) == 0) continue;
            if (line.rfind("Pod", 0) == 0) continue; 

            std::stringstream ss(line);
            std::string token;
            int col = 0;
            double t_raw = -1.0;
            std::string coin_str = "";
            
            while (std::getline(ss, token, ';')) {
                if (col == 1) {
                    try { t_raw = std::stod(token); } catch (...) {}
                } else if (col == 2) {
                    coin_str = token;
                    break;
                }
                col++;
            }

            if (t_raw >= 0 && coin_str.length() == 6) {
                if (is_first_line) {
                    local_t0 = t_raw; 
                    prev_t_raw = t_raw;
                    local_segment_offset = 0.0;
                    is_first_line = false;
                }
                
                // WYKRYCIE RESETU MIKROKONTROLERA
                if (t_raw < prev_t_raw) {
                    local_segment_offset += (prev_t_raw - local_t0);
                    local_t0 = t_raw;
                    
                    // Zapisujemy globalny czas wystąpienia tego resetu do narysowania zielonej linii
                    double reset_global_time = global_time_offset + (local_segment_offset / 3600000.0);
                    internal_resets.push_back(reset_global_time);
                }
                
                double t_hours = (local_segment_offset + (t_raw - local_t0)) / 3600000.0;
                double t_global = global_time_offset + t_hours;
                
                prev_t_raw = t_raw;

                bool is_4   = (coin_str == "000000");
                bool is_3   = (coin_str == "001011" || coin_str == "111000");
                bool is_3f  = (is_3 || coin_str == "010101" || coin_str == "100110");
                bool is_34  = (is_3 || is_4);
                bool is_34f = (is_3f || is_4);

                t_all.push_back(t_global);
                if (is_4)   t_4.push_back(t_global);
                if (is_3)   t_3.push_back(t_global);
                if (is_3f)  t_3f.push_back(t_global);
                if (is_34)  t_34.push_back(t_global);
                if (is_34f) t_34f.push_back(t_global);

                count_this_file++;
                if (t_hours > max_t_this_file) max_t_this_file = t_hours;
            }
        }
        file.close();

        std::cout << "File " << i << ": Events: " << count_this_file 
                  << " | Duration: " << max_t_this_file << " hours." << std::endl;

        if (count_this_file > 0) {
            double start_of_file = global_time_offset;
            global_time_offset += max_t_this_file;
            double end_of_file = global_time_offset;
            
            file_centers.push_back(start_of_file + (end_of_file - start_of_file) / 2.0);
            file_numbers.push_back(i);
            
            if (i < END) boundaries.push_back(end_of_file);
        }
    }

    if (t_all.empty()) {
        std::cerr << "[ERROR] No valid data found to plot!" << std::endl;
        return;
    }

    double total_time = global_time_offset;
    int n_bins = std::max(1, (int)std::ceil(total_time / bin_width_hours));

    // ================= CREATE HISTOGRAMS =================
    TH1F *h_all = new TH1F("h_all", "", n_bins, 0, total_time);
    TH1F *h_34f = new TH1F("h_34f", "", n_bins, 0, total_time);
    TH1F *h_34  = new TH1F("h_34", "", n_bins, 0, total_time);
    TH1F *h_3f  = new TH1F("h_3f", "", n_bins, 0, total_time);
    TH1F *h_3   = new TH1F("h_3", "", n_bins, 0, total_time);
    TH1F *h_4   = new TH1F("h_4", "", n_bins, 0, total_time);

    for (double t : t_all) h_all->Fill(t);
    for (double t : t_34f) h_34f->Fill(t);
    for (double t : t_34)  h_34->Fill(t);
    for (double t : t_3f)  h_3f->Fill(t);
    for (double t : t_3)   h_3->Fill(t);
    for (double t : t_4)   h_4->Fill(t);

    double scale_factor = 1.0 / bin_width_hours;
    h_all->Scale(scale_factor);
    h_34f->Scale(scale_factor);
    h_34->Scale(scale_factor);
    h_3f->Scale(scale_factor);
    h_3->Scale(scale_factor);
    h_4->Scale(scale_factor);

    // ================= STYLING =================
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    auto style_hist = [](TH1F* h, int color) {
        h->SetLineColor(color);
        h->SetLineWidth(2);
        h->SetFillStyle(0);
        h->GetXaxis()->SetTitle("Total Experiment Time [hours]");
        h->GetYaxis()->SetTitle("Count Rate [Counts / hour]");
        h->GetYaxis()->SetTitleOffset(1.0);
    };

    style_hist(h_all, kBlack);
    style_hist(h_34f, kRed);
    style_hist(h_34,  kBlue);
    style_hist(h_3f,  kGreen + 2);
    style_hist(h_3,   kMagenta);
    style_hist(h_4,   kCyan + 1);

    TCanvas *c1 = new TCanvas("c1", "Plots", 1200, 650);
    c1->SetLeftMargin(0.08);
    c1->SetBottomMargin(0.18); 
    c1->SetRightMargin(0.05);

    // ================= DRAWING FUNCTION =================
    auto save_plot = [&](std::vector<TH1F*> hists, std::vector<std::string> names, std::string filename, std::string title) {
        c1->Clear(); 
        
        double y_max = hists[0]->GetMaximum();
        double y_top = (y_max > 0) ? y_max * 1.15 : 1.0; 
        
        hists[0]->SetTitle(title.c_str());
        hists[0]->SetMinimum(0);
        hists[0]->SetMaximum(y_top);
        hists[0]->Draw("HIST"); 

        // 1. Rysowanie czerwonych linii (granice plików)
        for (double b_time : boundaries) {
            TLine *line = new TLine(b_time, 0, b_time, y_top);
            line->SetLineColor(kRed);
            line->SetLineWidth(2);
            line->SetLineStyle(2); 
            line->Draw("SAME");
        }

        // 2. Rysowanie zielonych przerywanych linii (wewnętrzne restarty w pliku)
        for (double r_time : internal_resets) {
            TLine *line = new TLine(r_time, 0, r_time, y_top);
            // kGreen+2 daje nieco ciemniejszy zielony, który lepiej widać na białym tle niż jaskrawy kGreen
            line->SetLineColor(kGreen + 2); 
            line->SetLineWidth(2);
            line->SetLineStyle(2); // Linia przerywana
            line->Draw("SAME");
        }

        // Dodanie numerów serii na górze
        for (size_t i = 0; i < file_centers.size(); ++i) {
            TLatex *tex = new TLatex();
            tex->SetTextSize(0.035);
            tex->SetTextAlign(23); 
            tex->SetTextColor(kBlack);
            tex->DrawLatex(file_centers[i], y_top * 0.98, std::to_string(file_numbers[i]).c_str());
        }

        for (size_t i = 1; i < hists.size(); ++i) {
            hists[i]->Draw("HIST SAME");
        }

        TLegend *leg = new TLegend(0.10, 0.01, 0.95, 0.08);
        leg->SetNColumns(hists.size()); 
        leg->SetBorderSize(0);          
        leg->SetFillStyle(0);           
        leg->SetTextSize(0.04);
        
        for (size_t i = 0; i < hists.size(); ++i) {
            leg->AddEntry(hists[i], names[i].c_str(), "l");
        }
        leg->Draw("SAME");

        c1->SaveAs(filename.c_str());
        std::cout << "Saved: " << filename << std::endl;
    };

    // ================= GENERATE PLOTS =================
    std::cout << "\nGenerating plots..." << std::endl;

    save_plot({h_all}, {"All data 2-3-4"}, path + "rate_only_ALL.pdf", "All recorded 2,3,4 coincidences over time");
    save_plot({h_34f}, {"3-4f"}, path + "rate_34f_per_hour.pdf", "3-4 (f) coincidences over time");
    save_plot({h_34},  {"3-4"},  path + "rate_34_per_hour.pdf",  "3-4 coincidences over time");
    save_plot({h_3f},  {"3f"},  path + "rate_3f_per_hour.pdf",  "3 (f) coincidences over time");
    save_plot({h_3},   {"3"},   path + "rate_3_per_hour.pdf",   "3 coincidences over time");
    save_plot({h_4},   {"4"},   path + "rate_4_per_hour.pdf",   "4 coincidences over time");

    save_plot(
        {h_all, h_34f, h_34, h_3f, h_3, h_4}, 
        {"All", "34f", "34", "3f", "3", "4"}, 
        path + "rate_all_combined.pdf", 
        "All data and coincidences"
    );

    save_plot(
        {h_34f, h_34, h_3f, h_3, h_4}, 
        {"34f", "34", "3f", "3", "4"}, 
        path + "rate_only_coincidences.pdf", 
        "Coincidences over time (background excluded)"
    );

    std::cout << "\nFinished! Plots are located in: " << path << std::endl;
}