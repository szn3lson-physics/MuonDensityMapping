#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

struct ModelData {
    std::string name;
    std::string filepath;
    std::vector<double> v_radius;
    std::vector<double> v_red_chi2;
    double min_chi2 = 1e18;
    double best_radius = 0;
    int color;
};

void load_model_data(ModelData& md) {
    std::ifstream file(md.filepath);
    if (!file.is_open()) {
        std::cerr << "[ERROR] File not found: " << md.filepath << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); 

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double radius, chi2_raw, df;

        if (iss >> radius >> chi2_raw >> df) {
            if (df <= 0) df = 1; 
            double red_chi2 = chi2_raw / df;
            
            md.v_radius.push_back(radius);
            md.v_red_chi2.push_back(red_chi2);

            if (red_chi2 < md.min_chi2) {
                md.min_chi2 = red_chi2;
                md.best_radius = radius;
            }
        }
    }
    file.close();
}

void plot_distance_unorm() {
    const char* OUTPUT_PLOT = "/home/kacper/MuonDensityMapping/figures/fit_distance_900.pdf";
    std::string base_dir = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/900/distance_analysis/";

    std::vector<ModelData> models = {
        {"Instant",     base_dir + "instant/distance_one_by_one_analysis.txt", {}, {}, 1e18, 0, kBlack},
        {"Exponential", base_dir + "exp/distance_one_by_one_analysis.txt",     {}, {}, 1e18, 0, kRed},
        {"Gaussian",    base_dir + "gauss/distance_one_by_one_analysis.txt",   {}, {}, 1e18, 0, TColor::GetColor("#0072B2")}
    };

    for (auto& md : models) {
        load_model_data(md);
        std::cout << "Loaded " << md.name << " | Minimum: R=" << md.best_radius << ", Chi2/DF=" << md.min_chi2 << std::endl;
    }

    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132, "t");
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetLabelFont(132, "xyz");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);

    // Wyższe płótno (960px) + wyliczone marginesy = fizyczny obszar wykresu 560x560 px (idealny kwadrat)
    TCanvas* c1 = new TCanvas("c1", "Chi2 Optimization Combined", 800, 960);
    c1->SetGrid(1, 1);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.366); // Potężny margines na dole specjalnie pod legendę

    TMultiGraph* mg = new TMultiGraph();
    //mg->SetTitle("Detection Range Optimization Using #chi^{2} Goodness-of-Fit Test");

    // Legenda umieszczona w strefie marginesu dolnego
    TLegend* leg = new TLegend(0.15, 0.05, 0.85, 0.25); 
    leg->SetTextFont(132);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0); // Bez ramki, by ładnie leżała pod spodem
    leg->SetFillStyle(0);

    for (auto& md : models) {
        if (md.v_radius.empty()) continue;
        
        int nPoints = md.v_radius.size();
        TGraph* gr = new TGraph(nPoints);
        for (int i = 0; i < nPoints; ++i) {
            gr->SetPoint(i, md.v_radius[i], md.v_red_chi2[i]);
        }
        
        gr->Sort();
        gr->SetLineColor(md.color);
        gr->SetLineWidth(3); 

        mg->Add(gr, "L");

        std::string leg_text = Form("%s (Best: R=%.1f m, #chi^{2}/df=%.3f)", md.name.c_str(), md.best_radius, md.min_chi2);
        leg->AddEntry(gr, leg_text.c_str(), "l"); 
    }

    mg->Draw("AL");
    c1->Update(); 

    // ================= CUSTOM AXIS RANGE ALGORITHM =================
    double min_x = 0.0, max_x = 100.0;
    double min_y = 1e18, max_y = -1e18;

    for (const auto& md : models) {
        for (size_t i = 0; i < md.v_radius.size(); ++i) {
            double x = md.v_radius[i];
            double y = md.v_red_chi2[i];
            if (x >= min_x && x <= max_x) {
                if (y < min_y) min_y = y;
                if (y > max_y) max_y = y;
            }
        }
    }
    
    double tick_offset = (max_y - min_y) * 0.10; 
    double y_bottom = std::max(0.0, min_y - tick_offset);

    mg->GetXaxis()->SetLimits(min_x, max_x);
    mg->SetMinimum(y_bottom);
    mg->SetMaximum(max_y); 
    // ===============================================================

    mg->GetXaxis()->SetTitle("Scanning Radius / Attenuation parameter R_{0} [m]");
    mg->GetYaxis()->SetTitle("Reduced #chi^{2} (#chi^{2} / df)");
    mg->GetXaxis()->SetTitleSize(0.04);
    mg->GetYaxis()->SetTitleSize(0.04);
    mg->GetXaxis()->SetLabelSize(0.03);
    mg->GetYaxis()->SetLabelSize(0.03);
    mg->GetYaxis()->SetTitleOffset(1.4);
    mg->GetXaxis()->SetTitleOffset(1.2);

    for (const auto& md : models) {
        if (md.v_radius.empty() || md.best_radius > max_x) continue;
        TLine *lx = new TLine(md.best_radius, y_bottom, md.best_radius, md.min_chi2);
        lx->SetLineColor(md.color);
        lx->SetLineStyle(2); 
        lx->SetLineWidth(2);
        lx->Draw("SAME");
    }

    leg->Draw("SAME");
    c1->SaveAs(OUTPUT_PLOT);
    std::cout << "Plot generated: " << OUTPUT_PLOT << std::endl;
}