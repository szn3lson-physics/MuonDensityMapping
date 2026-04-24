#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TMath.h>
#include <TPad.h>

// -------------------------------------------------------------------------
// Funkcje pomocnicze
// -------------------------------------------------------------------------
std::string VecToStr(const std::vector<int>& v) {
    std::string s;
    for (size_t i = 0; i < v.size(); ++i) {
        s += std::to_string(v[i]);
        if (i != v.size() - 1) s += "_";
    }
    return s;
}

std::string VecToTitle(const std::vector<int>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) {
        s += std::to_string(v[i]);
        if (i != v.size() - 1) s += ",";
    }
    s += "}";
    return s;
}

// -------------------------------------------------------------------------
// Ładowanie danych
// -------------------------------------------------------------------------
void LoadData(const std::vector<int>& series_list, int bin_width_deg, 
              std::vector<double>& sum_counts, std::vector<double>& sum_time_h) {
    
    int n_bins = 180 / bin_width_deg;
    sum_counts.assign(n_bins, 0.0);
    sum_time_h.assign(n_bins, 0.0);

    for (int series_id : series_list) {
        std::string input_file = "/home/kacper/MuonDensityMapping/output/single/dane_" 
                                 + std::to_string(series_id) + "/processed/time_per_angle_" 
                                 + std::to_string(series_id) + ".txt";
        
        std::fstream file(input_file, std::ios::in);
        
        if (!file.is_open()) {
            std::cerr << "[OSTRZEZENIE] Nie znaleziono pliku: " << input_file << std::endl;
            continue; 
        }

        std::string line;
        std::getline(file, line); 

        while (std::getline(file, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string item;
            int angle = 0; double time_min = 0.0, counts = 0.0;
            
            if (std::getline(ss, item, ';')) angle = std::stoi(item);
            if (std::getline(ss, item, ';')) time_min = std::stod(item);
            if (std::getline(ss, item, ';')) counts = std::stod(item);

            int bin_idx = angle / bin_width_deg;
            if (bin_idx >= 0 && bin_idx < n_bins) {
                sum_counts[bin_idx] += counts;
                sum_time_h[bin_idx] += (time_min / 60.0); 
            }
        }
        file.close();
    }
}

// -------------------------------------------------------------------------
// GŁÓWNA FUNKCJA MAKRA
// -------------------------------------------------------------------------
void his_compare_rotations() {
    
    // =========================================================================
    // KONFIGURACJA ZESTAWÓW DANYCH
    // =========================================================================
    std::vector<int> data_1 = {8};
    std::vector<int> data_2 = {4};
    int bin_width = 5;
    // =========================================================================

    // --- Ustawienia Globalne ---
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(kTRUE);     
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);         
    gStyle->SetPadTickY(1);
    gStyle->SetTextFont(132);

    int n_bins = 180 / bin_width;

    std::vector<double> counts_1, time_1, counts_2, time_2;
    LoadData(data_1, bin_width, counts_1, time_1);
    LoadData(data_2, bin_width, counts_2, time_2);

    TH1F *h1 = new TH1F("h1", "Comparison Data 1 vs Data 2;Zenith angle [#circ];Muon Rate [Counts / Hour]", n_bins, 0, 180);
    TH1F *h2 = new TH1F("h2", "", n_bins, 0, 180);
    
    h1->Sumw2();
    h2->Sumw2();

    for (int i = 0; i < n_bins; ++i) {
        if (time_1[i] > 0) {
            h1->SetBinContent(i + 1, counts_1[i] / time_1[i]);
            h1->SetBinError(i + 1, TMath::Sqrt(counts_1[i]) / time_1[i]);
        }
        if (time_2[i] > 0) {
            h2->SetBinContent(i + 1, counts_2[i] / time_2[i]);
            h2->SetBinError(i + 1, TMath::Sqrt(counts_2[i]) / time_2[i]);
        }
    }

    // --- STYLIZACJA GŁÓWNA ZGODNA Z PLOT_DISCOVERY ---
    
    // Data 1 (Sygnał/Czerwony)
    h1->SetLineColor(kRed + 1);
    h1->SetLineWidth(3);
    h1->SetFillColorAlpha(kRed, 0.3);
    
    // Data 2 (Referencja/Niebieski)
    h2->SetLineColor(kBlue + 1);
    h2->SetLineWidth(3);
    h2->SetFillColorAlpha(kBlue, 0.3);

    TCanvas *c1 = new TCanvas("c1", "Porownanie Roznicy", 1000, 900);
    
    // Panel górny
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0.02); 
    pad1->SetLeftMargin(0.12);
    pad1->Draw();
    pad1->cd();

    double max_val = TMath::Max(h1->GetMaximum(), h2->GetMaximum());
    h1->SetMaximum(max_val * 1.5); 
    h1->SetMinimum(0.0);

    h1->GetXaxis()->SetLabelSize(0); // Ukrywamy podpisy osi X na górnym panelu
    h1->GetYaxis()->SetTitle("Muon Rate [Counts / Hour]");
    h1->GetYaxis()->SetTitleSize(0.06);
    h1->GetYaxis()->SetTitleOffset(0.8);
    
    // Rysowanie z wypełnieniem ("HIST")
    h1->Draw("HIST");       
    h2->Draw("HIST SAME"); 
    h1->Draw("HIST SAME"); // Rysujemy h1 jeszcze raz na wierzch, by zachować odpowiednią warstwę

    TLegend *leg = new TLegend(0.6, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1, ("Data 1: Serie " + VecToTitle(data_1)).c_str(), "f");
    leg->AddEntry(h2, ("Data 2: Serie " + VecToTitle(data_2)).c_str(), "f");
    leg->Draw();

    // --- PANEL DOLNY (RÓŻNICA: DATA 1 - DATA 2) ---
    c1->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.12);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1F *h_diff = (TH1F*)h1->Clone("h_diff");
    // Obliczenie różnicy (Data 1 minus Data 2)
    h_diff->Add(h2, -1.0); 

    // USUWANIE WYPEŁNIENIA (TO NAPRAWIA PROBLEM!)
    h_diff->SetFillStyle(0);
    h_diff->SetFillColor(0);
    
    // Czysta czarna linia
    h_diff->SetLineColor(kBlack);
    h_diff->SetLineWidth(2);
    h_diff->SetTitle(""); 

    h_diff->GetXaxis()->SetTitle("Zenith angle [#circ]");
    h_diff->GetXaxis()->SetTitleSize(0.12);
    h_diff->GetXaxis()->SetLabelSize(0.10);
    h_diff->GetXaxis()->SetTitleOffset(1.0);
    
    h_diff->GetYaxis()->SetTitle("Data 1 - Data 2");
    h_diff->GetYaxis()->SetTitleSize(0.09);
    h_diff->GetYaxis()->SetLabelSize(0.08);
    h_diff->GetYaxis()->SetTitleOffset(0.5);
    h_diff->GetYaxis()->SetNdivisions(505);
    
    // --- Symetryczne Skalowanie wokół Zera ---
    double diff_max = TMath::Abs(h_diff->GetMaximum());
    double diff_min = TMath::Abs(h_diff->GetMinimum());
    double limit = TMath::Max(diff_max, diff_min) * 1.3;
    
    // Zabezpieczenie przed totalnie płaskim wykresem, jeśli serie są identyczne
    if (limit == 0) limit = 1.0; 

    h_diff->SetMaximum(limit);
    h_diff->SetMinimum(-limit);
    
    // Rysowanie samej linii ("HIST")
    h_diff->Draw("HIST");

    // Linia referencyjna dokładnie w punkcie zero
    TLine *line0 = new TLine(0, 0, 180, 0);
    line0->SetLineColor(kGray+2);
    line0->SetLineStyle(2);
    line0->Draw("SAME");

    // --- ZAPIS PLIKU ---
    std::string out_dir = "/home/kacper/MuonDensityMapping/results/relative/";
    gSystem->mkdir(out_dir.c_str(), kTRUE);
    
    std::string out_name = out_dir + "difference_" + VecToStr(data_1) + "_vs_" + VecToStr(data_2) + "_" + std::to_string(bin_width) + "deg.pdf";
    c1->SaveAs(out_name.c_str());

    std::cout << "[SUCCESS] Wykres zostal zapisany jako: " << out_name << std::endl;
}