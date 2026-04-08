#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

void plot_height() {

    // ================= INPUT / OUTPUT SETTINGS =================
    const char* FILE_PATH =
        "D:\\MuonDensityMapping\\bin\\simulation\\simulation_output\\height_analysis.txt";

    const char* OUTPUT_PLOT =
        "D:\\MuonDensityMapping\\figures\\fit_height.pdf";

    // ================= DATA LOADING =================
    std::ifstream file(FILE_PATH);
    if (!file.is_open()) {
        std::cerr << "[ERROR] File not found: "
                  << FILE_PATH << std::endl;
        return;
    }

    std::string line;

    // Skip header
    std::getline(file, line);

    std::vector<double> v_height, v_chi2;

    // Read file line by line
    while (std::getline(file, line)) {

        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);

        double h, chi2;

        // Read first two columns (ignore remaining columns if present)
        if (iss >> h >> chi2) {
            v_height.push_back(h);
            v_chi2.push_back(chi2);
        }
    }

    file.close();

    int nPoints = v_height.size();

    if (nPoints == 0) {
        std::cerr << "[ERROR] File is empty or contains no valid data."
                  << std::endl;
        return;
    }

    // ================= CREATE GRAPH =================
    TGraph* gr = new TGraph(nPoints,
                            v_height.data(),
                            v_chi2.data());

    gr->Sort();  // Ensure sorting by height

    // ================= STYLE SETTINGS =================
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132, "t");
    gStyle->SetTitleFont(132, "x");
    gStyle->SetTitleFont(132, "y");
    gStyle->SetLabelFont(132, "x");
    gStyle->SetLabelFont(132, "y");

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(2);
    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "Height Chi2 Standard", 1400, 800);
    c1->SetGrid(1, 1);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);

    Int_t color1 = TColor::GetColor("#0072B2");

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.4);
    gr->SetMarkerColor(color1);
    gr->SetLineColor(color1);
    gr->SetLineWidth(2);

    // ================= AXIS FORMATTING =================
    gr->SetTitle("#chi^{2} Goodness-of-Fit Test for Different Cavity Height Assumptions");

    gr->GetXaxis()->SetTitle("Cavity Height H [meters] (D = 37 m, Bin = 15#circ)");
    gr->GetYaxis()->SetTitle("Standard #chi^{2} (Fit to Detector Data)");

    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);

    gr->GetYaxis()->SetTitleOffset(1.15);
    gr->GetXaxis()->SetTitleOffset(1.1);

    // ================= X-AXIS PADDING =================
    double min_h = *std::min_element(v_height.begin(), v_height.end());
    double max_h = *std::max_element(v_height.begin(), v_height.end());
    double padding_x = (max_h - min_h) * 0.05;

    gr->Draw("APL");

    gr->GetXaxis()->SetLimits(min_h - padding_x,
                              max_h + padding_x);

    // ================= LEGEND =================
    TLegend* leg = new TLegend(0.65, 0.80, 0.88, 0.88);
    leg->SetTextFont(132);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(1);

    leg->AddEntry(gr, "Tested heights", "pl");
    leg->Draw("SAME");

    // ================= SAVE AS VECTOR PDF =================
    c1->SaveAs(OUTPUT_PLOT);

    std::cout << "Plot successfully generated and saved as: "
              << OUTPUT_PLOT << std::endl;
}
