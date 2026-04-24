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

// Funkcja generująca siatki wykresów oraz zbiorcze podsumowanie dla JEDNEJ serii pomiarowej
void generate_grid_for_series(int series_id, int bin_width_deg = 3, double sigma_input = 3.0) {
    
    // Ustawienia globalne stylu
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    std::string base_dir = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/single_rotation/" + std::to_string(series_id) + "_series/";
    std::string out_dir = "/home/kacper/MuonDensityMapping/results/rotation/" + std::to_string(series_id) + "_series/";
    // Tworzenie folderu docelowego
    gSystem->mkdir(out_dir.c_str(), kTRUE);

    int n_bins = 180 / bin_width_deg;
    int current_rot_start = 1;

    // --- TWORZENIE GŁÓWNEGO PŁÓTNA NA PODSUMOWANIA ---
    // Będzie to jedna duża siatka 4x2 przechowująca do 8 podsumowań.
    TCanvas *c_all_summaries = new TCanvas(Form("c_all_sum_s%d", series_id), 
                                           Form("Wszystkie podsumowania Serii %d", series_id), 
                                           2000, 1000);
    c_all_summaries->Divide(4, 2);
    int summary_pad_idx = 1; // Licznik okienek na podsumowania

    // Pętla paczkująca (co 8 obrotów)
    while (true) {
        int start_rot = current_rot_start;
        int end_rot = start_rot + 7; // Paczka 8 obrotów (np. 1-8, 9-16)

        // Sprawdzanie, czy w paczce są dane
        bool batch_has_data = false;
        for (int rot = start_rot; rot <= end_rot; ++rot) {
            std::string check_file = base_dir + std::to_string(rot) + "_rotation.txt";
            std::fstream file(check_file, std::ios::in);
            if (file.is_open()) {
                batch_has_data = true;
                file.close();
                break;
            }
        }

        if (!batch_has_data) break; // Brak danych = koniec serii pomiarowej

        // Zmienne do podsumowania całej paczki (8 obrotów)
        std::vector<double> batch_sum_counts(n_bins, 0.0);
        std::vector<double> batch_sum_time(n_bins, 0.0);

        // --- Tworzenie płótna 4x2 dla bieżącej paczki POJEDYNCZYCH obrotów ---
        TCanvas *c1 = new TCanvas(Form("c_grid_s%d_%d_%d", series_id, start_rot, end_rot), 
                                  Form("Series %d Rotations %d-%d", series_id, start_rot, end_rot), 
                                  2000, 1000); 
        c1->Divide(4, 2); 

        for (int rot = start_rot; rot <= end_rot; ++rot) {
            int pad_idx = rot - start_rot + 1; // 1 do 8
            c1->cd(pad_idx);

            std::string input_file = base_dir + std::to_string(rot) + "_rotation.txt";
            std::fstream file(input_file, std::ios::in);

            if (!file.is_open()) {
                TLatex *empty_text = new TLatex();
                empty_text->SetTextSize(0.08);
                empty_text->SetTextAlign(22);
                empty_text->DrawLatexNDC(0.5, 0.5, Form("Brak danych (Obr %d)", rot));
                continue;
            }

            std::vector<double> sum_counts(n_bins, 0.0);
            std::vector<double> sum_time_min(n_bins, 0.0);
            std::string line;
            std::getline(file, line); // Pomijanie nagłówka

            // Wczytywanie danych
            while (std::getline(file, line)) {
                if (line.empty()) continue;
                std::stringstream ss(line);
                std::string item;
                int angle = 0;
                double time_min = 0.0, counts = 0.0;

                if (std::getline(ss, item, ';')) angle = std::stoi(item);
                if (std::getline(ss, item, ';')) time_min = std::stod(item);
                if (std::getline(ss, item, ';')) counts = std::stod(item);

                int bin_idx = angle / bin_width_deg; 
                if (bin_idx >= 0 && bin_idx < n_bins) {
                    sum_counts[bin_idx] += counts;
                    sum_time_min[bin_idx] += time_min;

                    // Zbieranie do globalnego podsumowania całej paczki
                    batch_sum_counts[bin_idx] += counts;
                    batch_sum_time[bin_idx] += time_min;
                }
            }
            file.close();

            // Rysowanie histogramu dla pojedynczego obrotu
            TH1F *hist = new TH1F(Form("h_s%d_r%d", series_id, rot), Form("Obrot %d;Katy [#circ];Zlicz/h", rot), n_bins, 0, 180);
            
            double total_rate = 0.0, global_counts = 0.0, global_time_hours = 0.0;
            for (int i = 0; i < n_bins; ++i) {
                double rate = 0.0;
                if (sum_time_min[i] > 0) {
                    double time_hours = sum_time_min[i] / 60.0;
                    rate = sum_counts[i] / time_hours;
                    global_counts += sum_counts[i];
                    global_time_hours += time_hours;
                }
                hist->SetBinContent(i + 1, rate);
                total_rate += rate;
            }

            double mean_y = total_rate / n_bins;
            double avg_counts = global_counts / n_bins;
            double avg_time_h = global_time_hours / n_bins;
            double sigma_1 = (avg_time_h > 0) ? TMath::Sqrt(avg_counts) / avg_time_h : 0.0;
            double sigma = sigma_input * sigma_1;

            hist->SetLineColor(kRed);
            hist->SetLineWidth(2);
            hist->SetFillColorAlpha(kRed, 0.3);
            hist->SetMaximum(TMath::Max(hist->GetMaximum() * 1.15, mean_y + sigma * 1.2));
            hist->SetMinimum(0);
            hist->Draw("HIST");

            double x_min = hist->GetXaxis()->GetXmin();
            double x_max = hist->GetXaxis()->GetXmax();

            TLine *mean_y_line = new TLine(x_min, mean_y, x_max, mean_y);
            mean_y_line->SetLineColor(kBlue);
            mean_y_line->SetLineWidth(2);
            mean_y_line->SetLineStyle(2);
            mean_y_line->Draw("SAME");

            double y_high = mean_y + sigma;
            double y_low  = TMath::Max(0.0, mean_y - sigma);

            TLine *sig_high = new TLine(x_min, y_high, x_max, y_high);
            sig_high->SetLineColor(kBlack); sig_high->SetLineWidth(1); sig_high->SetLineStyle(3); sig_high->Draw("SAME");

            TLine *sig_low = new TLine(x_min, y_low, x_max, y_low);
            sig_low->SetLineColor(kBlack); sig_low->SetLineWidth(1); sig_low->SetLineStyle(3); sig_low->Draw("SAME");

            TLatex *text = new TLatex();
            text->SetTextSize(0.04);
            text->SetTextColor(kBlue);
            text->SetTextAlign(31); 
            text->DrawLatexNDC(0.85, 0.85, Form("#mu: %.1f/h", mean_y));
        }

        // 1. Zapis pliku z SIATKĄ POJEDYNCZYCH OBROTÓW 4x2
        std::string output_grid_pdf = out_dir + std::to_string(start_rot) + "_" + std::to_string(end_rot) + "_grid.pdf";
        c1->SaveAs(output_grid_pdf.c_str());
        delete c1;

        // --- 2. TWORZENIE ZIELONEGO PODSUMOWANIA NA WSPÓLNYM PŁÓTNIE ---
        // Zabezpieczenie przed przekroczeniem 8 okienek na siatce 4x2
        if (summary_pad_idx <= 8) {
            c_all_summaries->cd(summary_pad_idx);

            TH1F *hist_sum = new TH1F(Form("h_sum_s%d_%d_%d", series_id, start_rot, end_rot),
                                      Form("Podsumowanie %d obrotow (%d - %d);Katy zenitalne [#circ];Zliczenia 34f / h", 
                                           (end_rot - start_rot + 1), start_rot, end_rot),
                                      n_bins, 0, 180);

            double total_rate_sum = 0.0, global_counts_sum = 0.0, global_time_hours_sum = 0.0;
            
            for (int i = 0; i < n_bins; ++i) {
                double rate = 0.0;
                if (batch_sum_time[i] > 0) {
                    double time_hours = batch_sum_time[i] / 60.0;
                    rate = batch_sum_counts[i] / time_hours; 
                    global_counts_sum += batch_sum_counts[i];
                    global_time_hours_sum += time_hours;
                }
                hist_sum->SetBinContent(i + 1, rate);
                total_rate_sum += rate;
            }

            double mean_y_sum = total_rate_sum / n_bins;
            double avg_counts_sum = global_counts_sum / n_bins;
            double avg_time_h_sum = global_time_hours_sum / n_bins;
            double sigma_1_sum = (avg_time_h_sum > 0) ? TMath::Sqrt(avg_counts_sum) / avg_time_h_sum : 0.0;
            double sigma_sum = sigma_input * sigma_1_sum;

            hist_sum->SetLineColor(kGreen+2);
            hist_sum->SetLineWidth(3);
            hist_sum->SetFillColorAlpha(kGreen+2, 0.3);
            hist_sum->SetMaximum(TMath::Max(hist_sum->GetMaximum() * 1.15, mean_y_sum + sigma_sum * 1.2));
            hist_sum->SetMinimum(0);
            hist_sum->Draw("HIST");

            double x_min_sum = hist_sum->GetXaxis()->GetXmin();
            double x_max_sum = hist_sum->GetXaxis()->GetXmax();

            TLine *mean_line_sum = new TLine(x_min_sum, mean_y_sum, x_max_sum, mean_y_sum);
            mean_line_sum->SetLineColor(kBlue); mean_line_sum->SetLineWidth(3); mean_line_sum->SetLineStyle(2);
            mean_line_sum->Draw("SAME");

            double y_high_sum = mean_y_sum + sigma_sum;
            double y_low_sum  = TMath::Max(0.0, mean_y_sum - sigma_sum);

            TLine *sig_high_sum = new TLine(x_min_sum, y_high_sum, x_max_sum, y_high_sum);
            sig_high_sum->SetLineColor(kBlack); sig_high_sum->SetLineWidth(2); sig_high_sum->SetLineStyle(3); sig_high_sum->Draw("SAME");

            TLine *sig_low_sum = new TLine(x_min_sum, y_low_sum, x_max_sum, y_low_sum);
            sig_low_sum->SetLineColor(kBlack); sig_low_sum->SetLineWidth(2); sig_low_sum->SetLineStyle(3); sig_low_sum->Draw("SAME");

            TLatex *text_sum = new TLatex();
            text_sum->SetTextSize(0.04); text_sum->SetTextColor(kBlue); text_sum->SetTextAlign(31);
            text_sum->DrawLatexNDC(0.85, 0.85, Form("Mean: %.1f / h", mean_y_sum));
            text_sum->SetTextColor(kBlack); text_sum->SetTextAlign(11);
            text_sum->DrawLatexNDC(0.15, 0.85, Form("#pm%.1f#sigma", sigma_input));

            summary_pad_idx++; // Przejście do kolejnego okienka dla następnej paczki
        }

        // Przejście do kolejnych 8 obrotów w bieżącej serii
        current_rot_start += 8;
    }

    // --- ZAPISUJEMY WSPÓLNE PŁÓTNO PODSUMOWAŃ PO ZAKOŃCZENIU SERII ---
    if (summary_pad_idx > 1) { // Upewniamy się, że narysowaliśmy chociaż jedno podsumowanie
        std::string output_all_sum_pdf = out_dir + "series_" + std::to_string(series_id) + "_all_summaries.pdf";
        c_all_summaries->SaveAs(output_all_sum_pdf.c_str());
    }
    
    // Zwalniamy z pamięci
    delete c_all_summaries;
}

// Funkcja główna uruchamiana w konsoli ROOT
void his_rotation() {
    std::vector<int> series_list = {11};
    int bin_width = 3;  
    double sigma_level = 3.0; 
    
    std::cout << "Rozpoczynam automatyczne generowanie wszystkich obrotow i podsumowan..." << std::endl;

    for (int series_id : series_list) {
        std::cout << "--> Przetwarzanie Serii: " << series_id << "..." << std::endl;
        generate_grid_for_series(series_id, bin_width, sigma_level);
    }

    std::cout << "\nZakonczono sukcesem!" << std::endl;
    std::cout << "Pliki (siatki *_grid.pdf i zbiorcze *_all_summaries.pdf) zapisano w: " << std::endl;
    std::cout << "/home/kacper/MuonDensityMapping/results/rotation/" << std::endl;
}