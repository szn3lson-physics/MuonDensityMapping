#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TColor.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

void plot_distance() {

    // ================= INPUT / OUTPUT SETTINGS =================
    const char* FILE_PATH =
        //"D:\\MuonDensityMapping\\bin\\simulation\\simulation_output\\distance_one_by_one_analysis.txt";
        "D:\\MuonDensityMapping\\bin\\simulation\\simulation_output\\distance_loop_analysis.txt";
    const char* OUTPUT_PLOT =
        //"D:\\MuonDensityMapping\\figures\\fit_distance_one_by_one.pdf";
        "D:\\MuonDensityMapping\\figures\\fit_distance_loop.pdf";

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

    std::vector<double> v_radius, v_red_chi2;

    // Safe line-by-line reading
    while (std::getline(file, line)) {

        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);

        double radius, chi2_raw, df;

        if (iss >> radius >> chi2_raw >> df) {
            v_radius.push_back(radius);
            v_red_chi2.push_back(chi2_raw / df);
        }
    }

    file.close();

    int nPoints = v_radius.size();

    if (nPoints == 0) {
        std::cerr << "[ERROR] File is empty or contains no valid data."
                  << std::endl;
        return;
    }

    // ================= CREATE GRAPH =================
    TGraph* gr = new TGraph(nPoints);

    for (int i = 0; i < nPoints; ++i)
        gr->SetPoint(i, v_radius[i], v_red_chi2[i]);

    gr->Sort();  // Sort by radius

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

    TCanvas* c1 = new TCanvas("c1", "Chi2 Optimization", 1400, 800);
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
    gr->SetTitle("Detection Range Optimization Using #chi^{2} Goodness-of-Fit Test");

    gr->GetXaxis()->SetTitle("Scanning Radius R [meters]");
    gr->GetYaxis()->SetTitle("Reduced #chi^{2} (#chi^{2} / df)");

    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);

    gr->GetYaxis()->SetTitleOffset(1.1);
    gr->GetXaxis()->SetTitleOffset(1.1);

    gr->SetMinimum(0.0);

    gr->Draw("APL");  // Axes + Points + Line

    // ================= LEGEND =================
    TLegend* leg = new TLegend(0.65, 0.78, 0.88, 0.88);
    leg->SetTextFont(132);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(1);

    leg->AddEntry(gr, "Tested radii", "pl");
    leg->Draw("SAME");

    // ================= SAVE AS VECTOR PDF =================
    c1->SaveAs(OUTPUT_PLOT);

    std::cout << "Plot successfully generated and saved as: "
              << OUTPUT_PLOT << std::endl;
}
