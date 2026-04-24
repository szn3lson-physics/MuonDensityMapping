#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <fstream>
#include <sstream> 
#include <vector>
#include <string>
#include <iostream>

void plot_best_fit() {
    // ================= PARAMETERS AND PATHS =================
    // 1. Real data file (Z czasem)
    std::string real_data_path = "/home/kacper/MuonDensityMapping/output/all/zussamen_10/coin_time.txt"; 
    
    // 2. Directory containing Monte Carlo simulation files
    std::string sim_dir = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC/MC_one_by_one/";
    
    // 3. Array of 7 radii with specific requested values
    std::vector<double> radii = {10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0};
    double best_radius = 30.0; // This one will be drawn thick and black
    
    // Custom distinct colors for each radius (ROOT color codes)
    std::vector<int> line_colors = {
        kBlue,       // for 10.0
        kGreen+2,    // for 20.0
        kBlack,      // for 30.0 (Best fit)
        kMagenta,    // for 40.0
        kOrange+7,   // for 60.0
        kCyan+2,     // for 80.0
        kViolet-3    // for 100.0
    };

    // Histogram size
    int n_bins = 12;
    double x_min = 0;
    double x_max = 180;

    // ================= 1. LOADING REAL DATA =================
    TH1F *h_real = new TH1F("h_real", "Comparison of distance simulation with real data;Angle [degrees];Normalized counts (per 1000 min)", n_bins, x_min, x_max);
    
    std::ifstream f_real(real_data_path);
    if (f_real.is_open()) {
        std::string line;
        std::getline(f_real, line); // Pomiń nagłówek
        
        while (std::getline(f_real, line)) {
            if (line.empty()) continue;
            
            std::stringstream ss(line);
            std::string s_angle, s_time, s_counts;
            
            if (std::getline(ss, s_angle, ';') && 
                std::getline(ss, s_time, ';') && 
                std::getline(ss, s_counts, ';')) {
                
                try {
                    double angle = std::stod(s_angle);
                    double time_min = std::stod(s_time);
                    double counts = std::stod(s_counts);
                    
                    if (time_min > 0) {
                        // NORMALIZACJA CZASOWA do 1000 minut
                        double normalized_counts = counts * (1000.0 / time_min);
                        h_real->Fill(angle, normalized_counts); 
                    }
                } catch (const std::exception& e) { }
            }
        }
        f_real.close();
    } else {
        std::cout << "[ERROR] Cannot open real data file: " << real_data_path << std::endl;
    }

    // ================= 2. LOADING MC SIMULATION DATA =================
    std::vector<TH1F*> h_sims;
    TH1F* h_best = nullptr;

    double global_max = h_real->GetMaximum();
    double global_min = h_real->GetMinimum(0);

    for (size_t i = 0; i < radii.size(); ++i) {
        double r = radii[i];
        char buf[256];
        snprintf(buf, sizeof(buf), "%.2f_metrow.txt", r);
        std::string path = sim_dir + std::string(buf);
        
        TH1F *h_sim = new TH1F(Form("h_sim_%.0f", r), "", n_bins, x_min, x_max);
        std::fstream f_sim(path, std::ios::in);
        
        double val;
        if (f_sim.is_open()) {
            // WRACAMY DO STAREJ, DZIAŁAJĄCEJ PĘTLI DLA ZDARZEŃ MC
            while (f_sim >> val) {
                // Wycięcie ślepego kąta 0 stopni, żeby był fair wobec detektora
                if (val < 1.0) continue; 
                
                h_sim->Fill(val); 
            }
            f_sim.close();
        } else {
            std::cout << "[WARNING] Cannot open simulation file: " << path << std::endl;
        }
        
        if (h_sim->GetMaximum() > global_max) {
            global_max = h_sim->GetMaximum();
        }
        double current_min = h_sim->GetMinimum(0);
        if (current_min > 0 && current_min < global_min) {
            global_min = current_min;
        }
        
        if (r == best_radius) {
            h_best = h_sim;
        }
        
        h_sims.push_back(h_sim);
    }

    // ================= 3. MEAN/SIGMA CALCULATIONS (Real data) =================
    double manual_sigma = 2.0;
    double sum_of_content = h_real->Integral(); 
    double mean_y = sum_of_content / h_real->GetNbinsX(); 
    double sigma_y = manual_sigma * TMath::Sqrt(mean_y);
    double sigma = sigma_y;

    // ================= 4. GLOBAL STYLING =================
    gStyle->SetOptStat(0);       
    gStyle->SetPadGridX(kTRUE);  
    gStyle->SetPadGridY(kTRUE);  
    gStyle->SetPadTickX(1);      
    gStyle->SetPadTickY(1);

    h_real->SetLineColor(kRed);
    h_real->SetLineWidth(2);
    h_real->SetFillColorAlpha(kRed, 0.3);

    // ================= 5. DRAWING IN CORRECT ORDER =================
    TCanvas *c1 = new TCanvas("c1", "Fit Comparison", 1200, 700);
    c1->SetLeftMargin(0.08);
    c1->SetRightMargin(0.05);
    c1->SetBottomMargin(0.25); 

    double y_max_needed = TMath::Max(global_max * 1.15, mean_y + sigma * 1.5);
    double y_min_needed = global_min * 0.85; 
    
    h_real->SetMaximum(y_max_needed);
    h_real->SetMinimum(y_min_needed);

    h_real->Draw("HIST");

    int dash_counter = 0; 
    for (size_t i = 0; i < h_sims.size(); ++i) {
        if (radii[i] != best_radius) {
            h_sims[i]->SetLineColor(line_colors[i]); 
            h_sims[i]->SetLineWidth(2);              

            if (dash_counter % 2 == 1) {
                h_sims[i]->SetLineStyle(7); 
            } else {
                h_sims[i]->SetLineStyle(1); 
            }
            dash_counter++; 

            h_sims[i]->SetFillStyle(0);              
            h_sims[i]->Draw("HIST SAME");
        }
    }

    if (h_best) {
        h_best->SetLineColor(kBlack);
        h_best->SetLineWidth(4); 
        h_best->SetLineStyle(1); 
        h_best->SetFillStyle(0);
        h_best->Draw("HIST SAME");
    }

    TLine *mean_y_line = new TLine(x_min, mean_y, x_max, mean_y);
    mean_y_line->SetLineColor(kBlue);
    mean_y_line->SetLineWidth(2);
    mean_y_line->SetLineStyle(2);
    mean_y_line->Draw("SAME");
    
    TLatex *mean_y_text = new TLatex();
    mean_y_text->SetTextSize(0.035);
    mean_y_text->SetTextColor(kBlue);
    mean_y_text->SetTextAlign(31); 
    mean_y_text->DrawLatex(x_max * 0.98, mean_y + (y_max_needed * 0.02), Form("Mean (Y): %.2f", mean_y));

    double y_high = mean_y + sigma;
    TLine *sig_high = new TLine(x_min, y_high, x_max, y_high);
    sig_high->SetLineColor(kBlack);
    sig_high->SetLineWidth(1);
    sig_high->SetLineStyle(3); 
    sig_high->Draw("SAME");

    double y_low = mean_y - sigma;
    if (y_low < y_min_needed) y_low = y_min_needed; 
    TLine *sig_low = new TLine(x_min, y_low, x_max, y_low);
    sig_low->SetLineColor(kBlack);
    sig_low->SetLineWidth(1);
    sig_low->SetLineStyle(3); 
    sig_low->Draw("SAME");

    TLatex *sigma_text = new TLatex();
    sigma_text->SetTextSize(0.03);
    sigma_text->SetTextColor(kBlack);
    sigma_text->SetTextAlign(12); 
    
    sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01, y_high + ((y_max_needed - y_min_needed) * 0.02), Form("+%.1f#sigma", manual_sigma));
    if (y_low > y_min_needed) {
        sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01, y_low - ((y_max_needed - y_min_needed) * 0.02), Form("-%.1f#sigma", manual_sigma));
    }

    // ================= 6. LEGEND =================
    TLegend *leg = new TLegend(0.08, 0.02, 0.95, 0.15); 
    leg->SetNColumns(4); 
    leg->SetBorderSize(0); 
    
    leg->AddEntry(h_real, "Real data (normalized to 1000 min)", "f");
    
    if (h_best) {
        leg->AddEntry(h_best, Form("Best fit (%.0f m)", best_radius), "l");
    }
    
    for (size_t i = 0; i < h_sims.size(); ++i) {
        if (radii[i] != best_radius) {
            leg->AddEntry(h_sims[i], Form("Simulation %.0f m", radii[i]), "l");
        }
    }
    
    leg->Draw("SAME");

    // ================= 7. SAVE FILE =================
    c1->SaveAs("/home/kacper/MuonDensityMapping/figures/best_fit.pdf");
    
    std::cout << "\n[DONE] Saved file: /home/kacper/MuonDensityMapping/figures/best_fit.pdf" << std::endl;
}