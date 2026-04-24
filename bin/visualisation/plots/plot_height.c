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

struct ModelData {
    std::string name;
    std::string filepath;
    std::vector<double> v_height;
    std::vector<double> v_red_chi2;
    double min_chi2 = 1e18;
    double best_height = 0;
    int color;
};

void load_model_data(ModelData& md) {
    std::ifstream file(md.filepath);
    if (!file.is_open()) {
        std::cerr << "[ERROR] Nie znaleziono pliku: " << md.filepath << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Pomijamy nagłówek

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double h, chi2_raw, df;

        if (iss >> h >> chi2_raw >> df) {
            if (df <= 0) df = 1; 
            double red_chi2 = chi2_raw / df;
            
            md.v_height.push_back(h);
            md.v_red_chi2.push_back(red_chi2);

            if (red_chi2 < md.min_chi2) {
                md.min_chi2 = red_chi2;
                md.best_height = h;
            }
        }
    }
    file.close();
}

void plot_height() {
    const char* OUTPUT_PLOT = "/home/kacper/MuonDensityMapping/figures/fit_height_combined.pdf";
    std::string base_dir = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/hight_analysis/";

    // SPRAWDŹ TE ŚCIEŻKI: Upewnij się, że foldery instant, exp i gauss istnieją w height_analysis!
    std::vector<ModelData> models = {
        {"Hard Cutoff", base_dir + "instant/height_analysis.txt", {}, {}, 1e18, 0, kBlack},
        {"Exponential", base_dir + "exp/height_analysis.txt",     {}, {}, 1e18, 0, kRed},
        {"Gaussian",    base_dir + "gauss/height_analysis.txt",   {}, {}, 1e18, 0, TColor::GetColor("#0072B2")}
    };

    int loaded_models = 0;
    for (auto& md : models) {
        load_model_data(md);
        if (!md.v_height.empty()) {
            std::cout << "Loaded " << md.name << " | Minimum: H=" << md.best_height << "m, Chi2/DF=" << md.min_chi2 << std::endl;
            loaded_models++;
        }
    }

    // Zabezpieczenie przed rysowaniem pustego płótna
    if (loaded_models == 0) {
        std::cerr << "\n[KRYTYCZNY BLAD] Zadnen plik z danymi sie nie wczytal. Rysowanie przerwane.\n";
        return;
    }

    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132, "t");
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetLabelFont(132, "xyz");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "Height Optimization Combined", 1400, 800);
    c1->SetGrid(1, 1);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);

    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle("Cavity Height Optimization Using #chi^{2} Goodness-of-Fit Test");

    TLegend* leg = new TLegend(0.45, 0.70, 0.88, 0.88); 
    leg->SetTextFont(132);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(1);

    for (auto& md : models) {
        if (md.v_height.empty()) continue;
        
        int nPoints = md.v_height.size();
        TGraph* gr = new TGraph(nPoints);
        for (int i = 0; i < nPoints; ++i) {
            gr->SetPoint(i, md.v_height[i], md.v_red_chi2[i]);
        }
        
        gr->Sort();
        gr->SetLineColor(md.color);
        gr->SetLineWidth(3); // Grubsza linia dla lepszej widoczności

        mg->Add(gr, "L"); // Używamy "L" dla dodania jako linia

        std::string leg_text = Form("%s (Best: H=%.1f m, #chi^{2}/df=%.3f)", md.name.c_str(), md.best_height, md.min_chi2);
        leg->AddEntry(gr, leg_text.c_str(), "l"); 
    }

    // POPRAWKA: Wymuszenie linii ("AL")
    mg->Draw("AL");
    
    // POPRAWKA: Obowiązkowy update canvasu przed pytaniem o granice osi!
    c1->Update(); 

    mg->GetXaxis()->SetTitle("Cavity Height H [meters]");
    mg->GetYaxis()->SetTitle("Reduced #chi^{2} (#chi^{2} / df)");
    mg->GetXaxis()->SetTitleSize(0.045);
    mg->GetYaxis()->SetTitleSize(0.045);
    mg->GetXaxis()->SetLabelSize(0.04);
    mg->GetYaxis()->SetLabelSize(0.04);
    mg->GetYaxis()->SetTitleOffset(1.1);
    mg->GetXaxis()->SetTitleOffset(1.1);

    for (const auto& md : models) {
        if (md.v_height.empty()) continue;
        
        double y_bottom = mg->GetYaxis()->GetXmin(); // Teraz po Update() zadziała prawidłowo
        
        TLine *lx = new TLine(md.best_height, y_bottom, md.best_height, md.min_chi2);
        lx->SetLineColor(md.color);
        lx->SetLineStyle(2); 
        lx->SetLineWidth(2);
        lx->Draw("SAME");
    }

    leg->Draw("SAME");
    c1->SaveAs(OUTPUT_PLOT);
    std::cout << "Plot generated: " << OUTPUT_PLOT << std::endl;
}