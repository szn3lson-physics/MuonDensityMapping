#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>

void plot_temperature() {
    // ================= PARAMETERS =================
    const int START = 1;
    const int END   = 10; 
    
    const std::string INPUT_FILE_PREFIX = "/home/kacper/MuonDensityMapping/bin/data_analysis/data/dane_";
    const std::string OUTPUT_PATH = "/home/kacper/MuonDensityMapping/figures/rate_per_hour/";

    // ==============================================

    std::vector<double> time_hours;
    std::vector<double> temperatures;
    std::vector<double> pressures;
    
    std::vector<double> file_boundaries;  
    std::vector<double> internal_resets;  
    std::vector<double> file_centers;     
    std::vector<int> file_numbers;    

    double global_time_offset = 0.0;

    std::cout << "Rozpoczynam ekstrakcję danych o ciśnieniu i temperaturze..." << std::endl;

    for (int i = START; i <= END; ++i) {
        std::string filename = INPUT_FILE_PREFIX + std::to_string(i) + ".log";
        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "[WARNING] Nie można otworzyć pliku: " << filename << ". Pomijam..." << std::endl;
            continue;
        }

        std::string line;
        bool is_first_line = true;
        
        double local_t0 = 0.0;
        double prev_t_raw = 0.0;
        double local_segment_offset = 0.0;
        
        double max_t_this_file = 0.0;
        int points_this_file = 0;
        
        double start_of_file = global_time_offset;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            if (line.find("Adafruit") != std::string::npos || line.find("Pod") != std::string::npos) {
                continue;
            }

            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            
            while (std::getline(ss, token, ';')) {
                tokens.push_back(token);
            }
            
            if (tokens.size() >= 13) {
                double t_raw = -1.0;
                double pressure = -1.0;
                double temperature = -1000.0; 

                try {
                    t_raw = std::stod(tokens[1]);       
                    pressure = std::stod(tokens[11]);   
                    temperature = std::stod(tokens[12]);
                } catch (...) {
                    continue; 
                }

                if (t_raw >= 0) {
                    if (is_first_line) {
                        local_t0 = t_raw; 
                        prev_t_raw = t_raw;
                        local_segment_offset = 0.0;
                        is_first_line = false;
                    }
                    
                    if (t_raw < prev_t_raw) {
                        local_segment_offset += (prev_t_raw - local_t0);
                        local_t0 = t_raw;
                        
                        double reset_global_time = global_time_offset + (local_segment_offset / 3600000.0);
                        internal_resets.push_back(reset_global_time);
                    }
                    
                    double t_hours = (local_segment_offset + (t_raw - local_t0)) / 3600000.0;
                    double t_global = global_time_offset + t_hours;
                    
                    prev_t_raw = t_raw;
                    
                    time_hours.push_back(t_global);
                    pressures.push_back(pressure);
                    temperatures.push_back(temperature);

                    points_this_file++;
                    if (t_hours > max_t_this_file) max_t_this_file = t_hours;
                }
            }
        }
        file.close();

        std::cout << "Plik " << i << ": Zapisano " << points_this_file 
                  << " punktów | Czas trwania: " << max_t_this_file << " godzin." << std::endl;

        if (points_this_file > 0) {
            global_time_offset += max_t_this_file;
            double end_of_file = global_time_offset;
            
            file_centers.push_back(start_of_file + (end_of_file - start_of_file) / 2.0);
            file_numbers.push_back(i);
            
            if (i < END) file_boundaries.push_back(end_of_file);
        }
    }

    if (time_hours.empty()) {
        std::cerr << "[ERROR] Nie znaleziono żadnych poprawnych danych do narysowania!" << std::endl;
        return;
    }

    // ================= STYLING W ROOT =================
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    
    auto draw_markers = [&](double y_min, double y_max) {
        for (double b_time : file_boundaries) {
            TLine *line = new TLine(b_time, y_min, b_time, y_max);
            line->SetLineColor(kRed);
            line->SetLineWidth(2);
            line->SetLineStyle(2); 
            line->Draw("SAME");
        }

        for (size_t i = 0; i < file_centers.size(); ++i) {
            TLatex *tex = new TLatex();
            tex->SetTextSize(0.04);
            tex->SetTextAlign(23); 
            tex->SetTextColor(kRed);
            tex->DrawLatex(file_centers[i], y_max - (y_max - y_min) * 0.02, std::to_string(file_numbers[i]).c_str());
        }

        for (double r_time : internal_resets) {
            TLine *line = new TLine(r_time, y_min, r_time, y_max);
            line->SetLineColor(kGreen + 2);
            line->SetLineWidth(2);
            line->SetLineStyle(2); 
            line->Draw("SAME");
        }
    };

    // ================= WYKRES TEMPERATURY =================
    TCanvas *c_temp = new TCanvas("c_temp", "Temperatura", 1200, 600);
    c_temp->SetLeftMargin(0.10);
    c_temp->SetRightMargin(0.05);

    TGraph *g_temp = new TGraph(time_hours.size(), time_hours.data(), temperatures.data());
    g_temp->SetTitle("Zaleznosc temperatury od czasu;Czas [h];Temperatura [#circC]");
    g_temp->SetLineColor(kOrange);
    g_temp->SetLineWidth(2);
    
    // <--- NOWE: Ustawienie limitów osi dla temperatury
    g_temp->GetXaxis()->SetLimits(0, 3500); 
    g_temp->SetMinimum(19.5);               
    g_temp->SetMaximum(22.5);               
    
    g_temp->Draw("AL"); 
    draw_markers(19.5, 22.5); // <--- NOWE: Sztywne granice dla linii pomocniczych
    
    c_temp->SaveAs((OUTPUT_PATH + "temperatura_vs_czas.pdf").c_str());
    std::cout << "Zapisano wykres: " << OUTPUT_PATH << "temperatura_vs_czas.pdf" << std::endl;

    // ================= WYKRES CIŚNIENIA =================
    TCanvas *c_press = new TCanvas("c_press", "Cisnienie", 1200, 600);
    c_press->SetLeftMargin(0.10);
    c_press->SetRightMargin(0.05);

    TGraph *g_press = new TGraph(time_hours.size(), time_hours.data(), pressures.data());
    g_press->SetTitle("Zaleznosc cisnienia od czasu;Czas [h];Cisnienie [hPa]");
    g_press->SetLineColor(kBlue+1);
    g_press->SetLineWidth(2);
    
    // <--- NOWE: Ustawienie limitu osi X dla ciśnienia, aby wykresy były równe
    g_press->GetXaxis()->SetLimits(0, 3500);
    
    g_press->Draw("AL"); 
    
    // ROOT potrzebuje odświeżenia płótna, żeby poprawnie wyciągnąć automatyczne Minimum/Maximum z osi Y po jej narysowaniu
    c_press->Update(); 
    draw_markers(g_press->GetHistogram()->GetMinimum(), g_press->GetHistogram()->GetMaximum());

    c_press->SaveAs((OUTPUT_PATH + "cisnienie_vs_czas.pdf").c_str());
    std::cout << "Zapisano wykres: " << OUTPUT_PATH << "cisnienie_vs_czas.pdf" << std::endl;

    std::cout << "\nGotowe! Wykresy zapisano." << std::endl;
}