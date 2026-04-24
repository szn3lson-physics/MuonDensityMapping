#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

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
// Funkcja ładująca zliczenia z plików coin_X_34f.txt (sama kolumna kątów)
// -------------------------------------------------------------------------
void LoadRawData(const std::vector<int>& series_list, TH1F* hist) {
    for (int series_id : series_list) {
        // Zaktualizowana ścieżka do Twoich plików 34f
        std::string input_file = "/home/kacper/MuonDensityMapping/output/single/dane_" 
                                 + std::to_string(series_id) + "/processed/coin_" 
                                 + std::to_string(series_id) + "_34f.txt";
        
        std::fstream file(input_file, std::ios::in);
        
        if (!file.is_open()) {
            std::cerr << "[OSTRZEZENIE] Nie znaleziono pliku: " << input_file << std::endl;
            continue; 
        }

        double value;
        int licznik_zdarzen = 0;
        
        // Zliczanie na żywo
        while (file >> value) {
            hist->Fill(value);
            licznik_zdarzen++;
        }
        file.close();
        
        // Diagnostyka
        std::cout << "  -> Plik " << series_id << " zaladowal " << licznik_zdarzen << " zliczen." << std::endl;
    }
}

// -------------------------------------------------------------------------
// GŁÓWNA FUNKCJA MAKRA
// -------------------------------------------------------------------------
void his_compare_rotations() {
    
    // =========================================================================
    // KONFIGURACJA ZESTAWÓW DANYCH
    // =========================================================================
    std::vector<int> data_1 = {3};      // Pierwszy zestaw ("Sygnał")
    std::vector<int> data_2 = {8};   // Drugi zestaw ("Tło")
    int bin_width = 5;                         // Szerokość binu (w stopniach)
    // =========================================================================

    gStyle->SetOptStat(0);          
    gStyle->SetPadGridX(kTRUE);     
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);         
    gStyle->SetPadTickY(1);
    gStyle->SetTextFont(132); 

    int n_bins = 180 / bin_width;

    TH1F *h1 = new TH1F("h1", "Muography Detector-Background Analysis;Zenith angle [#circ];Muon Counts", n_bins, 0, 180);
    TH1F *h2 = new TH1F("h2", "", n_bins, 0, 180);
    TH1F *h_sig = new TH1F("h_sig", "", n_bins, 0, 180);

    h1->Sumw2();
    h2->Sumw2();

    std::cout << "\n=== Ladowanie Data 1 " << VecToTitle(data_1) << " ===" << std::endl;
    LoadRawData(data_1, h1);
    
    std::cout << "\n=== Ladowanie Data 2 " << VecToTitle(data_2) << " ===" << std::endl;
    LoadRawData(data_2, h2);

    double total_dat = h1->Integral();
    double total_bkg = h2->Integral();
    
    std::cout << "\n[INFO] Calkowita liczba w h1 (Data 1): " << total_dat << std::endl;
    std::cout << "[INFO] Calkowita liczba w h2 (Data 2): " << total_bkg << std::endl;

    if (total_dat == 0) {
        std::cerr << "\n[BLAD KRYTYCZNY] Histogram Data 1 jest calkowicie pusty!" << std::endl;
    }

    // --- WYSKALOWANIE TŁA (SCALING) ---
    if (total_bkg > 0 && total_dat > 0) {
        double scale_factor = total_dat / total_bkg;
        std::cout << "[INFO] Skalowanie Data 2 do Data 1. Wyliczony mnoznik: " << scale_factor << std::endl;
        h2->Scale(scale_factor);
    }

    // --- OBLICZANIE RÓŻNICY W SIGMACH ---
    for(int i = 1; i <= n_bins; i++) {
        double N1 = h1->GetBinContent(i);
        double N2 = h2->GetBinContent(i); 

        double diff = N1 - N2;
        double error = TMath::Sqrt(N1); 

        if(error > 0) {
            h_sig->SetBinContent(i, diff / error);
        } else {
            h_sig->SetBinContent(i, 0);
        }
    }

    // --- STYLIZACJA ---
    h1->SetLineColor(kRed + 1);
    h1->SetLineWidth(3);
    h1->SetFillColorAlpha(kRed, 0.3);
    
    h2->SetLineColor(kBlue + 1);
    h2->SetLineWidth(3);
    h2->SetFillColorAlpha(kBlue, 0.3);

    h_sig->SetLineColor(kBlack);
    h_sig->SetLineWidth(2);
    h_sig->SetFillStyle(0); // BRAK WYPEŁNIENIA NA DOLNYM PLOCIE

    TCanvas *c1 = new TCanvas("c1", "Porownanie Roznicy", 1000, 900);
    
    // --- TOP PAD ---
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0.02); 
    pad1->SetLeftMargin(0.12);
    pad1->Draw();
    pad1->cd();

    double max_val = TMath::Max(h1->GetMaximum(), h2->GetMaximum());
    if (max_val == 0) max_val = 1.0; 
    
    h1->SetMaximum(max_val * 1.5); 
    h1->SetMinimum(0.0);

    h1->GetXaxis()->SetLabelSize(0); 
    h1->GetYaxis()->SetTitle("Muon Counts");
    h1->GetYaxis()->SetTitleSize(0.06);
    h1->GetYaxis()->SetTitleOffset(0.8);
    
    h1->Draw("HIST");       
    h2->Draw("HIST SAME"); 
    h1->Draw("HIST SAME"); 

    TLegend *leg = new TLegend(0.6, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1, ("Data 1: Serie " + VecToTitle(data_1)).c_str(), "f");
    leg->AddEntry(h2, ("Data 2: (Scaled) Serie " + VecToTitle(data_2)).c_str(), "f");
    leg->Draw();

    // --- BOTTOM PAD ---
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
    
    h_sig->GetYaxis()->SetTitle("Data 1 - Data 2 [#sigma]");
    h_sig->GetYaxis()->SetTitleSize(0.09);
    h_sig->GetYaxis()->SetTitleOffset(0.5);
    h_sig->GetYaxis()->SetLabelSize(0.08);
    h_sig->GetYaxis()->SetNdivisions(505);

    double limit = TMath::Max(TMath::Abs(h_sig->GetMaximum()), TMath::Abs(h_sig->GetMinimum())) * 1.3;
    if (limit < 3.5) limit = 3.5; 
    
    h_sig->SetMaximum(limit);
    h_sig->SetMinimum(-limit);

    h_sig->Draw("HIST");

    TLine *line0 = new TLine(0, 0, 180, 0);
    line0->SetLineColor(kGray+2);
    line0->SetLineStyle(2);
    line0->Draw("SAME");

    TLine *line_p3 = new TLine(0, 3, 180, 3);
    line_p3->SetLineColor(kRed);
    line_p3->SetLineStyle(2);
    line_p3->Draw("SAME");

    TLine *line_m3 = new TLine(0, -3, 180, -3);
    line_m3->SetLineColor(kBlue);
    line_m3->SetLineStyle(2);
    line_m3->Draw("SAME");

    // --- ZAPIS PLIKU ---
    std::string out_dir = "/home/kacper/MuonDensityMapping/results/comparison_raw/";
    gSystem->mkdir(out_dir.c_str(), kTRUE);
    
    std::string out_name = out_dir + "raw_diff_" + VecToStr(data_1) + "_vs_" + VecToStr(data_2) + ".pdf";
    c1->SaveAs(out_name.c_str());

    std::cout << "\n[SUCCESS] Wykres zapisano jako: " << out_name << std::endl;
}