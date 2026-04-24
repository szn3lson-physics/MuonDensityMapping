#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLine.h>
#include <TLatex.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

void plot_distance_chi() {

    // ================= INPUT / OUTPUT SETTINGS =================
    const char* FILE_PATH = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/distance_analysis_prof/distance_one_by_one_analysis.txt";
    const char* OUTPUT_PLOT = "/home/kacper/MuonDensityMapping/figures/fit_exp_distance_chi2.pdf";

    // ================= DATA LOADING =================
    std::ifstream file(FILE_PATH);
    if (!file.is_open()) {
        std::cerr << "[ERROR] File not found: " << FILE_PATH << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Pominięcie nagłówka

    std::vector<double> v_radius, v_chi2;
    double min_chi2 = 1e18;
    double best_radius = 0;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        
        double radius, chi2_val, df, red_chi2, p_val;

        // Wczytujemy zgodnie z nowym formatem kolumn: Radius, Chi2, DF...
        if (iss >> radius >> chi2_val) {
            v_radius.push_back(radius);
            v_chi2.push_back(chi2_val);

            // Szukanie globalnego minimum
            if (chi2_val < min_chi2) {
                min_chi2 = chi2_val;
                best_radius = radius;
            }
        }
    }
    file.close();

    int nPoints = v_radius.size();
    if (nPoints == 0) return;

    // ================= STYLE SETTINGS =================
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132, "t");
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetLabelFont(132, "xyz");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(2);
    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "Chi2 Optimization", 1400, 800);
    c1->SetGrid(1, 1);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);

    // Tworzenie wykresu
    TGraph* gr = new TGraph(nPoints, &v_radius[0], &v_chi2[0]);
    gr->Sort();

    Int_t colorBlue = TColor::GetColor("#0072B2");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->SetMarkerColor(colorBlue);
    gr->SetLineColor(colorBlue);
    gr->SetLineWidth(2);

    // ================= AXIS FORMATTING =================
    gr->SetTitle("Detection Range Optimization Using #chi^{2} Goodness-of-Fit Test");
    gr->GetXaxis()->SetTitle("Scanning Radius R [meters]");
    gr->GetYaxis()->SetTitle("#chi^{2} Value"); // Zmienione z Reduced Chi2
    
    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetTitleOffset(1.1);

    // Ustawienie zakresu osi Y (od 0 do 110% wartości max dla przejrzystości)
    double yMax = *std::max_element(v_chi2.begin(), v_chi2.end());
    gr->SetMinimum(0.0);
    gr->SetMaximum(yMax * 1.1);

    gr->Draw("APL");

    // ================= MINIMUM MARKERS (Lines & Red Labels) =================
    
    // Czerwona linia do osi X
    TLine *lx = new TLine(best_radius, 0, best_radius, min_chi2);
    lx->SetLineColor(kRed);
    lx->SetLineStyle(2);
    lx->SetLineWidth(1);
    lx->Draw("SAME");

    // Czerwona linia do osi Y
    TLine *ly = new TLine(gr->GetXaxis()->GetXmin(), min_chi2, best_radius, min_chi2);
    ly->SetLineColor(kRed);
    ly->SetLineStyle(2);
    ly->SetLineWidth(1);
    ly->Draw("SAME");

    // Wartość na osi X (Radius)
    TLatex *texX = new TLatex();
    texX->SetTextColor(kRed);
    texX->SetTextFont(132);
    texX->SetTextSize(0.04);
    texX->SetTextAlign(23); // Centrowanie pod punktem
    // Przesunięcie lekko pod oś X
    double yOffset = yMax * 0.03;
    texX->DrawLatex(best_radius, -yOffset, Form("#bf{%.1f}", best_radius));

    // Wartość na osi Y (Chi2)
    TLatex *texY = new TLatex();
    texY->SetTextColor(kRed);
    texY->SetTextFont(132);
    texY->SetTextSize(0.04);
    texY->SetTextAlign(32); // Do prawej, wyrównane w pionie
    // Przesunięcie lekko na lewo od osi Y
    double xRange = gr->GetXaxis()->GetXmax() - gr->GetXaxis()->GetXmin();
    double xOffset = xRange * 0.02;
    texY->DrawLatex(gr->GetXaxis()->GetXmin() - xOffset, min_chi2, Form("#bf{%.2f}", min_chi2));

    // ================= SAVE =================
    c1->SaveAs(OUTPUT_PLOT);
    
    std::cout << "--- Analysis Results ---" << std::endl;
    std::cout << "Minimum Chi2: " << min_chi2 << " at R = " << best_radius << " m" << std::endl;
}