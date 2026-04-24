#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <filesystem>

#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>

namespace fs = std::filesystem;

struct CategoryData {
    std::vector<double> sum_d1, sum_d2, sum_d3, sum_d4;
    std::vector<int> count;

    void add(int bin_idx, double d1, double d2, double d3, double d4) {
        if (bin_idx >= sum_d1.size()) {
            int new_size = bin_idx + 1;
            sum_d1.resize(new_size, 0.0);
            sum_d2.resize(new_size, 0.0);
            sum_d3.resize(new_size, 0.0);
            sum_d4.resize(new_size, 0.0);
            count.resize(new_size, 0);
        }
        sum_d1[bin_idx] += d1;
        sum_d2[bin_idx] += d2;
        sum_d3[bin_idx] += d3;
        sum_d4[bin_idx] += d4;
        count[bin_idx]++;
    }
};

// Struktura do przechowywania średnich dla serii plików (teraz dla każdej kategorii)
struct FileMean {
    double start_t, end_t;
    double mean_all_d1, mean_all_d2, mean_all_d3, mean_all_d4;
    double mean_4_d1, mean_4_d2, mean_4_d3, mean_4_d4;
    double mean_3_d1, mean_3_d2, mean_3_d3, mean_3_d4;
    double mean_3f_d1, mean_3f_d2, mean_3f_d3, mean_3f_d4;
    double mean_34_d1, mean_34_d2, mean_34_d3, mean_34_d4;
    double mean_34f_d1, mean_34f_d2, mean_34f_d3, mean_34f_d4;
};

void plot_detector_rate() {
    const int START = 1;
    const int END   = 12; 
    const char* INPUT_FILE ="/home/kacper/MuonDensityMapping/bin/data_analysis/data/dane_";
    std::string path = "/home/kacper/MuonDensityMapping/figures/detector_frequency/";

    const double BIN_SIZE_MINUTES = 60.0; 
    const double BIN_SIZE_HOURS = BIN_SIZE_MINUTES / 60.0;

    if (!fs::exists(path)) fs::create_directories(path);

    CategoryData cat_all, cat_4, cat_3, cat_3f, cat_34, cat_34f;
    std::vector<FileMean> file_means;
    
    std::vector<double> boundaries;      
    std::vector<double> internal_resets; 
    std::vector<double> file_centers; 
    std::vector<int> file_numbers;    
    
    double global_time_ms = 0.0;

    for (int i = START; i <= END; ++i) {
        std::string filename = INPUT_FILE + std::to_string(i) + ".log";
        std::ifstream file(filename);
        if (!file.is_open()) continue;

        double file_start_hours = global_time_ms / 3600000.0;
        
        // Lokalne sumy i zliczenia dla każdego pliku, by obliczyć żółte linie na kategoriach
        double l_all_d1=0, l_all_d2=0, l_all_d3=0, l_all_d4=0; int c_all=0;
        double l_4_d1=0, l_4_d2=0, l_4_d3=0, l_4_d4=0; int c_4=0;
        double l_3_d1=0, l_3_d2=0, l_3_d3=0, l_3_d4=0; int c_3=0;
        double l_3f_d1=0, l_3f_d2=0, l_3f_d3=0, l_3f_d4=0; int c_3f=0;
        double l_34_d1=0, l_34_d2=0, l_34_d3=0, l_34_d4=0; int c_34=0;
        double l_34f_d1=0, l_34f_d2=0, l_34f_d3=0, l_34f_d4=0; int c_34f=0;

        std::string line;
        double prev_t_raw = -1.0;

        while (std::getline(file, line)) {
            if (line.empty() || line.rfind("Adafruit", 0) == 0 || line.rfind("Pod", 0) == 0) continue;

            std::stringstream ss(line);
            std::string token;
            int col = 0;
            double t_raw = -1, d1 = -1, d2 = -1, d3 = -1, d4 = -1;
            std::string coin_str = "";
            bool valid = false;
            
            while (std::getline(ss, token, ';')) {
                if (col == 1) { try { t_raw = std::stod(token); } catch(...) {} }
                else if (col == 2) { coin_str = token; }
                else if (col == 3) { try { d1 = std::stod(token); } catch(...) {} }
                else if (col == 4) { try { d2 = std::stod(token); } catch(...) {} }
                else if (col == 5) { try { d3 = std::stod(token); } catch(...) {} }
                else if (col == 6) { try { d4 = std::stod(token); valid = true; } catch(...) {} break; }
                col++;
            }

            if (valid && t_raw >= 0) {
                if (prev_t_raw < 0) prev_t_raw = t_raw; 
                double delta_t_raw = t_raw - prev_t_raw;
                if (delta_t_raw < 0) {
                    delta_t_raw = 0; prev_t_raw = t_raw;
                    internal_resets.push_back((global_time_ms + delta_t_raw) / 3600000.0);
                }
                global_time_ms += delta_t_raw;
                prev_t_raw = t_raw;

                double current_hours = global_time_ms / 3600000.0;
                int bin_idx = std::floor(current_hours / BIN_SIZE_HOURS);

                cat_all.add(bin_idx, d1, d2, d3, d4);
                l_all_d1 += d1; l_all_d2 += d2; l_all_d3 += d3; l_all_d4 += d4; c_all++;

                if (coin_str.length() == 6) {
                    bool is_4=(coin_str=="000000"), is_3=(coin_str=="001011"||coin_str=="111000");
                    bool is_3f=(is_3||coin_str=="010101"||coin_str=="100110");
                    bool is_34=(is_3||is_4), is_34f=(is_3f||is_4);
                    
                    if (is_4) {
                        cat_4.add(bin_idx, d1, d2, d3, d4);
                        l_4_d1+=d1; l_4_d2+=d2; l_4_d3+=d3; l_4_d4+=d4; c_4++;
                    }
                    if (is_3) {
                        cat_3.add(bin_idx, d1, d2, d3, d4);
                        l_3_d1+=d1; l_3_d2+=d2; l_3_d3+=d3; l_3_d4+=d4; c_3++;
                    }
                    if (is_3f) {
                        cat_3f.add(bin_idx, d1, d2, d3, d4);
                        l_3f_d1+=d1; l_3f_d2+=d2; l_3f_d3+=d3; l_3f_d4+=d4; c_3f++;
                    }
                    if (is_34) {
                        cat_34.add(bin_idx, d1, d2, d3, d4);
                        l_34_d1+=d1; l_34_d2+=d2; l_34_d3+=d3; l_34_d4+=d4; c_34++;
                    }
                    if (is_34f) {
                        cat_34f.add(bin_idx, d1, d2, d3, d4);
                        l_34f_d1+=d1; l_34f_d2+=d2; l_34f_d3+=d3; l_34f_d4+=d4; c_34f++;
                    }
                }
            }
        }
        file.close();

        if (c_all > 0) {
            double end_hours = global_time_ms / 3600000.0;
            
            FileMean fm;
            fm.start_t = file_start_hours; fm.end_t = end_hours;
            
            auto calc_mean = [](double sum, int c) { return (c > 0) ? sum/c : 0.0; };
            
            fm.mean_all_d1 = calc_mean(l_all_d1, c_all); fm.mean_all_d2 = calc_mean(l_all_d2, c_all);
            fm.mean_all_d3 = calc_mean(l_all_d3, c_all); fm.mean_all_d4 = calc_mean(l_all_d4, c_all);

            fm.mean_4_d1 = calc_mean(l_4_d1, c_4); fm.mean_4_d2 = calc_mean(l_4_d2, c_4);
            fm.mean_4_d3 = calc_mean(l_4_d3, c_4); fm.mean_4_d4 = calc_mean(l_4_d4, c_4);

            fm.mean_3_d1 = calc_mean(l_3_d1, c_3); fm.mean_3_d2 = calc_mean(l_3_d2, c_3);
            fm.mean_3_d3 = calc_mean(l_3_d3, c_3); fm.mean_3_d4 = calc_mean(l_3_d4, c_3);

            fm.mean_3f_d1 = calc_mean(l_3f_d1, c_3f); fm.mean_3f_d2 = calc_mean(l_3f_d2, c_3f);
            fm.mean_3f_d3 = calc_mean(l_3f_d3, c_3f); fm.mean_3f_d4 = calc_mean(l_3f_d4, c_3f);

            fm.mean_34_d1 = calc_mean(l_34_d1, c_34); fm.mean_34_d2 = calc_mean(l_34_d2, c_34);
            fm.mean_34_d3 = calc_mean(l_34_d3, c_34); fm.mean_34_d4 = calc_mean(l_34_d4, c_34);

            fm.mean_34f_d1 = calc_mean(l_34f_d1, c_34f); fm.mean_34f_d2 = calc_mean(l_34f_d2, c_34f);
            fm.mean_34f_d3 = calc_mean(l_34f_d3, c_34f); fm.mean_34f_d4 = calc_mean(l_34f_d4, c_34f);

            file_means.push_back(fm);
            
            file_centers.push_back(file_start_hours + (end_hours - file_start_hours) / 2.0);
            file_numbers.push_back(i);
            if (i < END) boundaries.push_back(end_hours);
        }
    }

    auto create_graphs = [&](CategoryData& cat, TGraph*& g1, TGraph*& g2, TGraph*& g3, TGraph*& g4) {
        g1 = new TGraph(); g2 = new TGraph(); g3 = new TGraph(); g4 = new TGraph();
        int pt_idx = 0;
        for (size_t i = 0; i < cat.sum_d1.size(); ++i) {
            if (cat.count[i] > 0) {
                double t = (i + 0.5) * BIN_SIZE_HOURS;
                g1->SetPoint(pt_idx, t, cat.sum_d1[i]/cat.count[i]);
                g2->SetPoint(pt_idx, t, cat.sum_d2[i]/cat.count[i]);
                g3->SetPoint(pt_idx, t, cat.sum_d3[i]/cat.count[i]);
                g4->SetPoint(pt_idx, t, cat.sum_d4[i]/cat.count[i]);
                pt_idx++;
            }
        }
    };

    TGraph *g_all_1, *g_all_2, *g_all_3, *g_all_4; create_graphs(cat_all, g_all_1, g_all_2, g_all_3, g_all_4);
    TGraph *g_4_1, *g_4_2, *g_4_3, *g_4_4;         create_graphs(cat_4, g_4_1, g_4_2, g_4_3, g_4_4);
    TGraph *g_3_1, *g_3_2, *g_3_3, *g_3_4;         create_graphs(cat_3, g_3_1, g_3_2, g_3_3, g_3_4);
    TGraph *g_3f_1, *g_3f_2, *g_3f_3, *g_3f_4;     create_graphs(cat_3f, g_3f_1, g_3f_2, g_3f_3, g_3f_4);
    TGraph *g_34_1, *g_34_2, *g_34_3, *g_34_4;     create_graphs(cat_34, g_34_1, g_34_2, g_34_3, g_34_4);
    TGraph *g_34f_1, *g_34f_2, *g_34f_3, *g_34f_4; create_graphs(cat_34f, g_34f_1, g_34f_2, g_34f_3, g_34f_4);

    gStyle->SetOptStat(0);
    auto style_graph = [](TGraph* g, int color) {
        g->SetLineColor(color); g->SetLineWidth(2);
        g->SetMarkerStyle(20); g->SetMarkerSize(0.6); g->SetMarkerColor(color);
        g->GetXaxis()->SetTitle("Experiment Time [hours]");
        g->GetYaxis()->SetTitle("Average Frequency [Hz]");
    };

    TCanvas *c1 = new TCanvas("c1", "Plots", 1400, 650);
    c1->SetBottomMargin(0.18); 

    // det_idx określa który konkretny detektor rysujemy (1, 2, 3 lub 4). 0 = wykres zbiorczy.
    // cat_type: 0=ALL, 1=4, 2=3, 3=3f, 4=34, 5=34f. Pozwala wyciągnąć odpowiednią średnią dla serii z FileMean.
    auto save_plot = [&](std::vector<TGraph*> graphs, std::vector<std::string> names, std::string filename, std::string title, int det_idx, int cat_type) {
        c1->Clear(); 
        double y_max = 0;
        for(auto g : graphs) if(g->GetN()>0) y_max = std::max(y_max, TMath::MaxElement(g->GetN(), g->GetY()));
        double y_top = y_max * 1.15; 
        if (y_top == 0) y_top = 100.0;
        
        graphs[0]->SetTitle(title.c_str());
        graphs[0]->SetMinimum(0); graphs[0]->SetMaximum(y_top);
        graphs[0]->Draw("ALP"); 

        // Rysowanie żółtej kreski średniej dla każdej serii i odpowiedniego detektora i kategorii
        if (det_idx != 0) {
            for (auto &m : file_means) {
                double val = 0;
                if (cat_type == 0) {
                    if (det_idx == 1) val = m.mean_all_d1; else if (det_idx == 2) val = m.mean_all_d2;
                    else if (det_idx == 3) val = m.mean_all_d3; else if (det_idx == 4) val = m.mean_all_d4;
                } else if (cat_type == 1) {
                    if (det_idx == 1) val = m.mean_4_d1; else if (det_idx == 2) val = m.mean_4_d2;
                    else if (det_idx == 3) val = m.mean_4_d3; else if (det_idx == 4) val = m.mean_4_d4;
                } else if (cat_type == 2) {
                    if (det_idx == 1) val = m.mean_3_d1; else if (det_idx == 2) val = m.mean_3_d2;
                    else if (det_idx == 3) val = m.mean_3_d3; else if (det_idx == 4) val = m.mean_3_d4;
                } else if (cat_type == 3) {
                    if (det_idx == 1) val = m.mean_3f_d1; else if (det_idx == 2) val = m.mean_3f_d2;
                    else if (det_idx == 3) val = m.mean_3f_d3; else if (det_idx == 4) val = m.mean_3f_d4;
                } else if (cat_type == 4) {
                    if (det_idx == 1) val = m.mean_34_d1; else if (det_idx == 2) val = m.mean_34_d2;
                    else if (det_idx == 3) val = m.mean_34_d3; else if (det_idx == 4) val = m.mean_34_d4;
                } else if (cat_type == 5) {
                    if (det_idx == 1) val = m.mean_34f_d1; else if (det_idx == 2) val = m.mean_34f_d2;
                    else if (det_idx == 3) val = m.mean_34f_d3; else if (det_idx == 4) val = m.mean_34f_d4;
                }
                
                if (val > 0) {
                    TLine *ml = new TLine(m.start_t, val, m.end_t, val);
                    ml->SetLineColor(kYellow); ml->SetLineWidth(3); ml->Draw("SAME");
                }
            }
        }

        for (double b : boundaries) { TLine *l=new TLine(b,0,b,y_top); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw("SAME"); }
        for (double r : internal_resets) { TLine *l=new TLine(r,0,r,y_top); l->SetLineColor(kGreen+2); l->SetLineStyle(2); l->Draw("SAME"); }
        for (size_t i=0; i<file_centers.size(); ++i) {
            TLatex *tex = new TLatex(); tex->SetTextSize(0.035); tex->SetTextAlign(23);
            tex->DrawLatex(file_centers[i], y_top*0.98, std::to_string(file_numbers[i]).c_str());
        }

        for (size_t i=1; i<graphs.size(); ++i) if(graphs[i]->GetN()>0) graphs[i]->Draw("LP SAME");

        TLegend *leg = new TLegend(0.10, 0.01, 0.95, 0.08); leg->SetNColumns(graphs.size()); leg->SetBorderSize(0);
        for (size_t i=0; i<graphs.size(); ++i) leg->AddEntry(graphs[i], names[i].c_str(), "lp");
        if (det_idx != 0) leg->AddEntry((TObject*)0, "--- Yellow Line: Series Mean ---", "");
        leg->Draw("SAME");

        c1->SaveAs(filename.c_str());
    };

    // Generowanie Stylów
    std::vector<std::string> d_names = {"Det 1", "Det 2", "Det 3", "Det 4"};
    auto style_group = [&](TGraph* g1, TGraph* g2, TGraph* g3, TGraph* g4) {
        style_graph(g1, kBlack); style_graph(g2, kRed); style_graph(g3, kBlue); style_graph(g4, kGreen+2);
    };

    style_group(g_all_1, g_all_2, g_all_3, g_all_4);
    style_group(g_4_1, g_4_2, g_4_3, g_4_4);
    style_group(g_3_1, g_3_2, g_3_3, g_3_4);
    style_group(g_3f_1, g_3f_2, g_3f_3, g_3f_4);
    style_group(g_34_1, g_34_2, g_34_3, g_34_4);
    style_group(g_34f_1, g_34f_2, g_34f_3, g_34f_4);

    std::cout << "\nGenerating individual detector plots with yellow means..." << std::endl;
    // INDYWIDUALNE: ALL (cat_type = 0)
    save_plot({g_all_1}, {"Det 1"}, path+"avg_freq_ALL_det1.pdf", "Avg Freq (All): Det 1", 1, 0);
    save_plot({g_all_2}, {"Det 2"}, path+"avg_freq_ALL_det2.pdf", "Avg Freq (All): Det 2", 2, 0);
    save_plot({g_all_3}, {"Det 3"}, path+"avg_freq_ALL_det3.pdf", "Avg Freq (All): Det 3", 3, 0);
    save_plot({g_all_4}, {"Det 4"}, path+"avg_freq_ALL_det4.pdf", "Avg Freq (All): Det 4", 4, 0);

    // INDYWIDUALNE: Koincydencje 4 (cat_type = 1)
    save_plot({g_4_1}, {"Det 1"}, path+"avg_freq_4_det1.pdf", "Avg Freq (Coinc 4): Det 1", 1, 1);
    save_plot({g_4_2}, {"Det 2"}, path+"avg_freq_4_det2.pdf", "Avg Freq (Coinc 4): Det 2", 2, 1);
    save_plot({g_4_3}, {"Det 3"}, path+"avg_freq_4_det3.pdf", "Avg Freq (Coinc 4): Det 3", 3, 1);
    save_plot({g_4_4}, {"Det 4"}, path+"avg_freq_4_det4.pdf", "Avg Freq (Coinc 4): Det 4", 4, 1);

    // INDYWIDUALNE: Koincydencje 3 (cat_type = 2)
    save_plot({g_3_1}, {"Det 1"}, path+"avg_freq_3_det1.pdf", "Avg Freq (Coinc 3): Det 1", 1, 2);
    save_plot({g_3_2}, {"Det 2"}, path+"avg_freq_3_det2.pdf", "Avg Freq (Coinc 3): Det 2", 2, 2);
    save_plot({g_3_3}, {"Det 3"}, path+"avg_freq_3_det3.pdf", "Avg Freq (Coinc 3): Det 3", 3, 2);
    save_plot({g_3_4}, {"Det 4"}, path+"avg_freq_3_det4.pdf", "Avg Freq (Coinc 3): Det 4", 4, 2);

    // INDYWIDUALNE: Koincydencje 3f (cat_type = 3)
    save_plot({g_3f_1}, {"Det 1"}, path+"avg_freq_3f_det1.pdf", "Avg Freq (Coinc 3f): Det 1", 1, 3);
    save_plot({g_3f_2}, {"Det 2"}, path+"avg_freq_3f_det2.pdf", "Avg Freq (Coinc 3f): Det 2", 2, 3);
    save_plot({g_3f_3}, {"Det 3"}, path+"avg_freq_3f_det3.pdf", "Avg Freq (Coinc 3f): Det 3", 3, 3);
    save_plot({g_3f_4}, {"Det 4"}, path+"avg_freq_3f_det4.pdf", "Avg Freq (Coinc 3f): Det 4", 4, 3);

    // INDYWIDUALNE: Koincydencje 34 (cat_type = 4)
    save_plot({g_34_1}, {"Det 1"}, path+"avg_freq_34_det1.pdf", "Avg Freq (Coinc 34): Det 1", 1, 4);
    save_plot({g_34_2}, {"Det 2"}, path+"avg_freq_34_det2.pdf", "Avg Freq (Coinc 34): Det 2", 2, 4);
    save_plot({g_34_3}, {"Det 3"}, path+"avg_freq_34_det3.pdf", "Avg Freq (Coinc 34): Det 3", 3, 4);
    save_plot({g_34_4}, {"Det 4"}, path+"avg_freq_34_det4.pdf", "Avg Freq (Coinc 34): Det 4", 4, 4);

    // INDYWIDUALNE: Koincydencje 34f (cat_type = 5)
    save_plot({g_34f_1}, {"Det 1"}, path+"avg_freq_34f_det1.pdf", "Avg Freq (Coinc 34f): Det 1", 1, 5);
    save_plot({g_34f_2}, {"Det 2"}, path+"avg_freq_34f_det2.pdf", "Avg Freq (Coinc 34f): Det 2", 2, 5);
    save_plot({g_34f_3}, {"Det 3"}, path+"avg_freq_34f_det3.pdf", "Avg Freq (Coinc 34f): Det 3", 3, 5);
    save_plot({g_34f_4}, {"Det 4"}, path+"avg_freq_34f_det4.pdf", "Avg Freq (Coinc 34f): Det 4", 4, 5);


    std::cout << "Generating combined detector plots..." << std::endl;
    // ZBIORCZE (det_idx = 0 -> brak żółtych linii, 4 detektory naraz)
    save_plot({g_all_1, g_all_2, g_all_3, g_all_4}, d_names, path+"avg_freq_ALL_combined.pdf", "Avg Freq (All): All Detectors", 0, 0);
    save_plot({g_4_1, g_4_2, g_4_3, g_4_4},         d_names, path+"avg_freq_COINC_4.pdf",       "Hz during Coinc 4",           0, 0);
    save_plot({g_3_1, g_3_2, g_3_3, g_3_4},         d_names, path+"avg_freq_COINC_3.pdf",       "Hz during Coinc 3",           0, 0);
    save_plot({g_3f_1, g_3f_2, g_3f_3, g_3f_4},     d_names, path+"avg_freq_COINC_3f.pdf",      "Hz during Coinc 3f",          0, 0);
    save_plot({g_34_1, g_34_2, g_34_3, g_34_4},     d_names, path+"avg_freq_COINC_34.pdf",      "Hz during Coinc 34",          0, 0);
    save_plot({g_34f_1, g_34f_2, g_34f_3, g_34f_4}, d_names, path+"avg_freq_COINC_34f.pdf",     "Hz during Coinc 34f",         0, 0);

    std::cout << "\nDone! All files are in: " << path << std::endl;
}