#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

void plot_best_fit() {
    // ================= PARAMETRY I ŚCIEŻKI =================
    // 1. Plik z danymi rzeczywistymi (upewnij się, że ścieżka jest poprawna)
    std::string real_data_path = "D:\\MuonDensityMapping\\output\\all\\zussamen_10\\coin_34f.txt"; 
    
    // 2. Folder, w którym C++ wygenerował pliki Monte Carlo
    std::string sim_dir = "D:\\MuonDensityMapping\\bin\\simulation\\simulation_output\\MC_one_by_one\\";
    
    // 3. Tablica 10 promieni, które chcemy narysować (1 czarny + 9 szarych)
    std::vector<double> radii = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
    double best_radius = 30.0; // Ten zostanie narysowany na czarno
    
    // Rozmiar histogramu
    int n_bins = 20;
    double x_min = 0;
    double x_max = 180;

    // ================= 1. WCZYTYWANIE DANYCH RZECZYWISTYCH =================
    TH1F *h_real = new TH1F("h_real", "Porownanie symulacji odleglosci z danymi rzeczywistymi;Kąt [stopnie];Liczba zliczen", n_bins, x_min, x_max);
    
    std::fstream f_real(real_data_path, std::ios::in);
    double val;
    while (f_real >> val) {
        h_real->Fill(val);
    }
    f_real.close();

    // ================= 2. WCZYTYWANIE DANYCH Z SYMULACJI MC =================
    std::vector<TH1F*> h_sims;
    TH1F* h_best = nullptr;

    // Zmienna trzymająca globalne maksimum, żeby oś Y pomieściła wszystkie wykresy
    double global_max = h_real->GetMaximum();

    for (double r : radii) {
        char buf[256];
        snprintf(buf, sizeof(buf), "%.2f_metrow.txt", r);
        std::string path = sim_dir + std::string(buf);
        
        TH1F *h_sim = new TH1F(Form("h_sim_%.0f", r), "", n_bins, x_min, x_max);
        std::fstream f_sim(path, std::ios::in);
        
        if (f_sim.is_open()) {
            while (f_sim >> val) {
                h_sim->Fill(val);
            }
            f_sim.close();
        } else {
            std::cout << "[OSTRZEZENIE] Nie mozna otworzyc pliku symulacji: " << path << std::endl;
        }
        
        // Aktualizacja maksimum osi Y
        if (h_sim->GetMaximum() > global_max) {
            global_max = h_sim->GetMaximum();
        }
        
        // Zapisanie wskaźnika do najlepszego fita, żeby go wyciągnąć potem na wierzch
        if (r == best_radius) {
            h_best = h_sim;
        }
        
        h_sims.push_back(h_sim);
    }

    // ================= 3. OBLICZENIA MEAN/SIGMA (Dla danych rzeczywistych) =================
    double manual_sigma = 2.0;
    double sum_of_content = h_real->Integral(); 
    double mean_y = sum_of_content / h_real->GetNbinsX(); 
    double sigma_y = manual_sigma * TMath::Sqrt(mean_y);
    double sigma = sigma_y;

    // ================= 4. STYLIZACJA GLOBALNA =================
    gStyle->SetOptStat(0);       
    gStyle->SetPadGridX(kTRUE);  
    gStyle->SetPadGridY(kTRUE);  
    gStyle->SetPadTickX(1);      
    gStyle->SetPadTickY(1);

    // Styl danych rzeczywistych
    h_real->SetLineColor(kRed);
    h_real->SetLineWidth(2);
    h_real->SetFillColorAlpha(kRed, 0.3); // Półprzezroczysty czerwony

    // ================= 5. RYSOWANIE W ODPOWIEDNIEJ KOLEJNOŚCI =================
    TCanvas *c1 = new TCanvas("c1", "Porownanie Fitow", 1200, 600);
    c1->SetLeftMargin(0.08);
    c1->SetRightMargin(0.05);
    c1->SetBottomMargin(0.12);

    double y_max_needed = TMath::Max(global_max * 1.15, mean_y + sigma * 1.5);
    h_real->SetMaximum(y_max_needed);
    h_real->SetMinimum(0);

    // 1. Najpierw dane rzeczywiste (tło/wypełnienie)
    h_real->Draw("HIST");

    // 2. Następnie wszystkie SZARE linie dla odrzuconych promieni
    for (size_t i = 0; i < h_sims.size(); ++i) {
        if (radii[i] != best_radius) {
            h_sims[i]->SetLineColor(kGray);
            h_sims[i]->SetLineWidth(1);
            h_sims[i]->SetFillStyle(0); // Brak wypełnienia
            h_sims[i]->Draw("HIST SAME");
        }
    }

    // 3. Na samym końcu gruby CZARNY wykres na wierzchu
    if (h_best) {
        h_best->SetLineColor(kBlack);
        h_best->SetLineWidth(3);
        h_best->SetFillStyle(0); // Brak wypełnienia
        h_best->Draw("HIST SAME");
    }

    // 4. Linie znaczącości narysowane na samym wierzchu
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
    if (y_low < 0) y_low = 0; 
    TLine *sig_low = new TLine(x_min, y_low, x_max, y_low);
    sig_low->SetLineColor(kBlack);
    sig_low->SetLineWidth(1);
    sig_low->SetLineStyle(3); 
    sig_low->Draw("SAME");

    TLatex *sigma_text = new TLatex();
    sigma_text->SetTextSize(0.03);
    sigma_text->SetTextColor(kBlack);
    sigma_text->SetTextAlign(12); 
    
    sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01, y_high + (y_max_needed * 0.02), Form("+%.1f#sigma", manual_sigma));
    if (y_low > 0) {
        sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01, y_low - (y_max_needed * 0.02), Form("-%.1f#sigma", manual_sigma));
    }

    // ================= 6. LEGENDA =================
    TLegend *leg = new TLegend(0.70, 0.72, 0.93, 0.88); // Prawy górny róg
    leg->AddEntry(h_real, "Dane rzeczywiste", "f");
    
    if (h_best) {
        leg->AddEntry(h_best, Form("Najlepszy fit (%.0f m)", best_radius), "l");
    }
    
    // Szukamy jednego szarego, żeby służył za wzór w legendzie
    for (size_t i = 0; i < h_sims.size(); ++i) {
        if (radii[i] != best_radius) {
            leg->AddEntry(h_sims[i], "Inne badane promienie", "l");
            break;
        }
    }
    
    leg->Draw("SAME");

    // Zapisz do pliku PDF / PNG
    c1->SaveAs("D:\\MuonDensityMapping\\porownanie_promieni.pdf");
    
    std::cout << "\nGotowe! Zapisano plik: D:\\MuonDensityMapping\\porownanie_promieni.pdf" << std::endl;
}