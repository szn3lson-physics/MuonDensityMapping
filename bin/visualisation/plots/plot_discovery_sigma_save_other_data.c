//====================================
// Kacper Dorszewski                 ||
// plot_muography_analysis.C         ||
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
#include <TPad.h>
#include <vector>

    const char* INPUT_FILE = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_3_with_hole.txt";
    const char* INPUT_BACKGROUND = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_3_without_hole.txt";

    const char* OUTPUT_PLOT = "/home/kacper/MuonDensityMapping/figures/residua_detektor_3.pdf";
    const char* OUTPUT_TEXT = "/home/kacper/MuonDensityMapping/bin/visualisation/output/residua_detektor_3.txt";

void plot_discovery_sigma_save_other_data() {
    // --- Global Style Settings ---
    gStyle->SetOptStat(0);          
    gStyle->SetPadGridX(kTRUE);     
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);         
    gStyle->SetPadTickY(1);
    gStyle->SetTextFont(132);

    int bin_width_deg = 5;
    int n_bins = 180 / bin_width_deg; // Poprawiono na 180, aby uwzględnić kąt 180 stopni
    double x_min = 0.0;               // Zmieniono na 0, aby (180-0)/180 dawało równy 1 stopień szerokości
    double x_max = 180.0;

    // --- Main Histograms ---
    TH1F *h_dat = new TH1F("h_dat", "Muography Detector-Background Data Analysis;Zenith angle [#circ];Muon Rate [Counts / Hour]", n_bins, x_min, x_max);
    TH1F *h_bkg = new TH1F("h_bkg", "Background", n_bins, x_min, x_max);
    TH1F *h_sig = new TH1F("h_sig", "", n_bins, x_min, x_max);
    
    // --- Wektory akumulacji danych ---
    std::vector<double> sum_counts_dat(n_bins, 0.0);
    std::vector<double> sum_counts_bkg(n_bins, 0.0);
    std::vector<double> sum_time_min(n_bins, 0.0);

    std::string line;

    // --- 1. LOAD DETECTOR DATA (Z dziurą) ---
    std::ifstream f_dat_txt(INPUT_FILE);
    if (f_dat_txt.is_open()) {
        while (std::getline(f_dat_txt, line)) {
            // Pomijamy puste linie i nagłówek
            if (line.empty() || line.find("Angle") != std::string::npos) continue;
            
            std::stringstream ss(line);
            std::string item;
            int angle = 0; double time_min = 0.0; double counts = 0.0;

            if (std::getline(ss, item, ';')) angle = std::stoi(item);
            if (std::getline(ss, item, ';')) time_min = std::stod(item);
            if (std::getline(ss, item, ';')) counts = std::stod(item);

            int bin_idx = (angle - 1) / bin_width_deg;
            // Warunek upewnia się, że pomijamy kąt 0 (bin_idx będzie < 0)
            if (bin_idx >= 0 && bin_idx < n_bins) {
                sum_counts_dat[bin_idx] += counts;
                sum_time_min[bin_idx] += time_min;
            }
        }
        f_dat_txt.close();
    } else {
        std::cerr << "[ERROR] Cannot open data file: " << INPUT_FILE << std::endl;
        return;
    }

    // --- 2. LOAD BACKGROUND MODEL (Bez dziury) ---
    std::ifstream f_bkg_txt(INPUT_BACKGROUND);
    if (f_bkg_txt.is_open()) {
        while (std::getline(f_bkg_txt, line)) {
            // Pomijamy puste linie i nagłówek
            if (line.empty() || line.find("Angle") != std::string::npos) continue;
            
            std::stringstream ss(line);
            std::string item;
            int ang = 0; double time_min_bkg = 0.0; double theo_counts = 0.0;

            // Zastosowano to samo czytanie średnikami co w przypadku Data
            if (std::getline(ss, item, ';')) ang = std::stoi(item);
            if (std::getline(ss, item, ';')) time_min_bkg = std::stod(item);
            if (std::getline(ss, item, ';')) theo_counts = std::stod(item);

            int bin_idx = (ang - 1) / bin_width_deg;
            if (bin_idx >= 0 && bin_idx < n_bins) {
                sum_counts_bkg[bin_idx] += theo_counts;
            }
        }
        f_bkg_txt.close();
    } else {
        std::cerr << "[ERROR] Cannot open background model file: " << INPUT_BACKGROUND << std::endl;
        return;
    }

    // --- 3. CALCULATE GLOBAL SIGMA ---
    double global_counts_dat = 0.0;
    double global_time_hours = 0.0;

    for (int i = 0; i < n_bins; ++i) {
        global_counts_dat += sum_counts_dat[i];
        if (sum_time_min[i] > 0) global_time_hours += (sum_time_min[i] / 60.0);
    }

    double avg_counts = global_counts_dat / n_bins;
    double avg_time_h = global_time_hours / n_bins;
    double global_sigma = 0.0;
    if (avg_time_h > 0) {
        global_sigma = TMath::Sqrt(avg_counts) / avg_time_h;
    }

    // --- 4. FILL HISTOGRAMS & WRITE TEXT FILE ---
    std::ofstream f_out_txt(OUTPUT_TEXT);
    if (f_out_txt.is_open()) {
        f_out_txt << "Angle;Time_h;Substraction;Significance\n";
    } else {
        std::cerr << "[ERROR] Cannot open output text file: " << OUTPUT_TEXT << std::endl;
    }

    for(int i = 0; i < n_bins; ++i) {
        double current_angle = (i * bin_width_deg) + 1.0; 
        
        if (sum_time_min[i] > 0) {
            double time_hours = sum_time_min[i] / 60.0;
            double rate_dat = sum_counts_dat[i] / time_hours;
            double rate_bkg = sum_counts_bkg[i] / time_hours;

            h_dat->SetBinContent(i + 1, rate_dat);
            h_bkg->SetBinContent(i + 1, rate_bkg);

            double diff = rate_dat - rate_bkg;
            double significance = 0.0;
            
            if (global_sigma > 0) {
                significance = diff / global_sigma;
                h_sig->SetBinContent(i + 1, significance);
            }
            
            if (f_out_txt.is_open()) {
                f_out_txt << current_angle << ";" << time_hours << ";" << diff << ";" << significance << "\n";
            }
            
        } else {
            if (f_out_txt.is_open()) {
                f_out_txt << current_angle << ";0;0;0\n";
            }
        }
    }
    
    if (f_out_txt.is_open()) {
        f_out_txt.close();
        std::cout << "[SUCCESS] Residua text file saved as " << OUTPUT_TEXT << std::endl;
    }

    // --- 5. GRAPHICAL FORMATTING ---
    h_dat->SetLineColor(kBlue + 1);
    h_dat->SetLineWidth(3);
    h_dat->SetFillColorAlpha(kBlue, 0.3);

    h_bkg->SetLineColor(kRed + 1);
    h_bkg->SetLineWidth(3);
    h_bkg->SetFillColorAlpha(kRed, 0.3);

    h_sig->SetLineColor(kBlack);
    h_sig->SetLineWidth(2);
    h_sig->SetFillStyle(0);

    // --- 6. CANVAS AND PADS ---
    TCanvas *c1 = new TCanvas("c1", "Muography Detector-Background Data Analysis", 1000, 900);

    // Top Pad
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0.02); 
    pad1->SetLeftMargin(0.12);
    pad1->Draw();
    pad1->cd();

    h_dat->SetMaximum(h_dat->GetMaximum() * 1.5);
    h_dat->SetMinimum(0.0);
    h_dat->GetXaxis()->SetLabelSize(0); 
    h_dat->GetYaxis()->SetTitle("Muon Rate [Counts / Hour]");
    h_dat->GetYaxis()->SetTitleSize(0.06);
    h_dat->GetYaxis()->SetTitleOffset(0.8);

    h_dat->Draw("HIST");        
    h_bkg->Draw("HIST SAME");   
    h_dat->Draw("HIST SAME");   

    TLegend *leg = new TLegend(0.6, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_dat, "Measured Data (Detector)", "f");
    leg->AddEntry(h_bkg, "Background Model", "f");
    leg->Draw();

    // Bottom Pad
    c1->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.12);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    h_sig->SetTitle("");
    h_sig->GetXaxis()->SetTitle("Zenith angle [#circ]");
    h_sig->GetXaxis()->SetTitleSize(0.12);
    h_sig->GetXaxis()->SetLabelSize(0.10);
    h_sig->GetXaxis()->SetTitleOffset(1.0);
    
    h_sig->GetYaxis()->SetTitle("Data - Bkg [#sigma]");
    h_sig->GetYaxis()->SetTitleSize(0.09);
    h_sig->GetYaxis()->SetTitleOffset(0.5);
    h_sig->GetYaxis()->SetLabelSize(0.08);
    h_sig->GetYaxis()->SetNdivisions(505);

    double limit = TMath::Max(TMath::Abs(h_sig->GetMaximum()), TMath::Abs(h_sig->GetMinimum())) * 1.3;
    if (limit < 3.5) limit = 3.5; 
    
    h_sig->SetMaximum(limit);
    h_sig->SetMinimum(-limit);
    h_sig->Draw("HIST");

    TLine *line0 = new TLine(x_min, 0, x_max, 0);
    line0->SetLineColor(kGray+2);
    line0->SetLineStyle(2);
    line0->Draw("SAME");

    TLine *line_p3 = new TLine(x_min, 3, x_max, 3);
    line_p3->SetLineColor(kRed);
    line_p3->SetLineStyle(2);
    line_p3->Draw("SAME");

    TLine *line_m3 = new TLine(x_min, -3, x_max, -3);
    line_m3->SetLineColor(kBlue);
    line_m3->SetLineStyle(2);
    line_m3->Draw("SAME");

    c1->SaveAs(OUTPUT_PLOT);
    std::cout << "[SUCCESS] Plot saved as " << OUTPUT_PLOT << std::endl;
}