#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TColor.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

void plot_bin() {

    // ================= USER SETTINGS =================
    const char* FILE_PATH  = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/bin_analysis/bin_analysis_small.txt";
    const char* OUTPUT_PLOT = "/home/kacper/MuonDensityMapping/figures/fit_bin_small.pdf";

    // ================= DATA LOADING =================
    std::ifstream file(FILE_PATH);
    if (!file.is_open()) {
        std::cerr << "[ERROR] File not found: " << FILE_PATH << std::endl;
        return;
    }

    // Skip header line
    std::string header;
    std::getline(file, header);

    std::vector<double> v_bin, v_chi2, v_err;
    double bin, col1, df, chi2red, pval;

    // Read five columns from file
    while (file >> bin >> col1 >> df >> chi2red >> pval) {
        v_bin.push_back(bin);
        v_chi2.push_back(chi2red);

        // Statistical uncertainty of reduced chi-square: sqrt(2 / df)
        v_err.push_back(std::sqrt(2.0 / df));
    }
    file.close();

    int nPoints = v_bin.size();
    if (nPoints == 0) {
        std::cerr << "[ERROR] File is empty or has invalid format." << std::endl;
        return;
    }

    // ================= CREATE GRAPH =================
    TGraphErrors* gr = new TGraphErrors(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        gr->SetPoint(i, v_bin[i], v_chi2[i]);
        gr->SetPointError(i, 0.0, v_err[i]);  // No X error, Y error from statistics
    }

    // Sort points along X axis
    gr->Sort();

    // ================= STYLE SETTINGS (PUBLICATION READY) =================
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

    // High-resolution canvas (not critical for PDF, but improves screen preview)
    TCanvas* c1 = new TCanvas("c1", "Bin Optimization", 1600, 1000);
    c1->SetGrid(1, 1);
    c1->SetBottomMargin(0.18);
    c1->SetLeftMargin(0.12);

    // Custom color (#0072B2)
    Int_t color1 = TColor::GetColor("#0072B2");

    gr->SetMarkerStyle(21);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerColor(color1);
    gr->SetLineColor(kGray+2);
    gr->SetLineWidth(2);

    gr->SetTitle("Bin Width Optimization Based on Goodness-of-Fit Test");
    gr->GetXaxis()->SetTitle("Bin Width #Delta#theta [degrees]");
    gr->GetYaxis()->SetTitle("Reduced #chi^{2} (#chi^{2} / df)");

    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetTitleOffset(1.2);
    gr->GetXaxis()->SetTitleOffset(1.2);

    gr->Draw("AP");

    // ================= IDEAL MODEL LINE =================
    double xmin = gr->GetXaxis()->GetXmin();
    double xmax = gr->GetXaxis()->GetXmax();

    TLine* ideal_line = new TLine(xmin, 1.0, xmax, 1.0);
    ideal_line->SetLineColor(kBlack);
    ideal_line->SetLineWidth(2);
    ideal_line->Draw("SAME");

    // ================= LEGEND =================
    TLegend* leg = new TLegend(0.2, 0.02, 0.8, 0.08);
    leg->SetNColumns(2);
    leg->SetBorderSize(0);
    leg->SetTextFont(132);
    leg->SetTextSize(0.04);

    leg->AddEntry(gr, "Reduced #chi^{2} (#pm 1#sigma)", "pe");
    leg->AddEntry(ideal_line, "Ideal model (#chi^{2}_{red} = 1)", "l");

    leg->Draw("SAME");

    // ================= SAVE AS VECTOR PDF =================
    c1->SaveAs(OUTPUT_PLOT);

    std::cout << "Plot successfully generated and saved as: "
              << OUTPUT_PLOT << std::endl;
}
