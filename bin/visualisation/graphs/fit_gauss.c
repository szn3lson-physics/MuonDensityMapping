#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TColor.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

void fit_gauss() {

    // ================= INPUT / OUTPUT SETTINGS =================
    const char* filename =
        "D:\\MuonDensityMapping\\output\\all\\10_zussamen\\coin__34f.txt";

    const char* OUTPUT_PLOT =
        "D:\\MuonDensityMapping\\figures\\gauss_fit.pdf";

    // ================= DATA LOADING =================
    std::vector<double> data;
    std::ifstream fin(filename);

    if (fin.is_open()) {
        double val;
        while (fin >> val) {
            data.push_back(fmod(val, 180.0));
        }
        fin.close();
        std::cout << "Loaded " << data.size()
                  << " events from file." << std::endl;
    } else {
        std::cout << "[WARNING] File not found. Generating test bimodal data..."
                  << std::endl;

        TRandom3 rand(0);

        // Background-like component (center ~40 deg)
        for (int i = 0; i < 20000; ++i)
            data.push_back(fmod(rand.Gaus(40, 10), 180.0));

        // Signal-like component (center ~120 deg)
        for (int i = 0; i < 30000; ++i)
            data.push_back(fmod(rand.Gaus(120, 20), 180.0));
    }

    // ================= HISTOGRAM =================
    // 12 bins between 0 and 180 degrees
    TH1F* h = new TH1F(
        "h",
        "Angular Distribution (Bimodal Gaussian Fit);Angle [degrees];Probability Density",
        12, 0, 180);

    for (double v : data)
        h->Fill(v);

    // Normalize to probability density (area = 1)
    if (h->Integral() > 0)
        h->Scale(1.0 / h->Integral("width"));

    // ================= BASIC STATISTICS =================
    double mu     = h->GetMean();
    double sigma  = h->GetStdDev();
    double skew   = h->GetSkewness();
    double kurt   = h->GetKurtosis();  // Fisher definition (Gaussian = 0)

    std::cout << "\nHistogram statistics:" << std::endl;
    std::cout << "Mean: " << mu
              << " | Std Dev: " << sigma << std::endl;
    std::cout << "Skewness: " << skew << std::endl;
    std::cout << "Kurtosis: " << kurt << std::endl;

    // ================= FIT: SUM OF TWO GAUSSIANS =================
    TF1* bimodal = new TF1("bimodal", "gaus(0) + gaus(3)", 0, 180);

    double expected_A = h->GetMaximum() / 2.0;

    bimodal->SetParameters(expected_A, mu - sigma, sigma/2.0,
                           expected_A, mu + sigma, sigma/2.0);

    bimodal->SetParLimits(1, 0, 180);
    bimodal->SetParLimits(4, 0, 180);

    bimodal->SetLineColor(kRed);
    bimodal->SetLineWidth(3);

    h->Fit(bimodal, "S 0");

    double A1     = bimodal->GetParameter(0);
    double mu1    = bimodal->GetParameter(1);
    double sigma1 = bimodal->GetParameter(2);

    double A2     = bimodal->GetParameter(3);
    double mu2    = bimodal->GetParameter(4);
    double sigma2 = bimodal->GetParameter(5);

    std::cout << "\nFitted Gaussian 1: mu=" << mu1
              << ", sigma=" << sigma1 << std::endl;
    std::cout << "Fitted Gaussian 2: mu=" << mu2
              << ", sigma=" << sigma2 << std::endl;

    TF1* g1 = new TF1("g1", "gaus", 0, 180);
    g1->SetParameters(A1, mu1, sigma1);
    g1->SetLineColor(kGreen+2);
    g1->SetLineStyle(2);
    g1->SetLineWidth(2);

    TF1* g2 = new TF1("g2", "gaus", 0, 180);
    g2->SetParameters(A2, mu2, sigma2);
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
    g2->SetLineWidth(2);

    // ================= VISUALIZATION =================
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(132);

    TCanvas* c1 = new TCanvas("c1", "Bimodal Fit", 1400, 800);
    c1->SetGrid(1, 1);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);

    h->SetFillColorAlpha(kMagenta+2, 0.4);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);

    double max_y = std::max(h->GetMaximum(),
                            bimodal->GetMaximum(0, 180));
    h->SetMaximum(max_y * 1.2);

    h->Draw("HIST");
    g1->Draw("SAME");
    g2->Draw("SAME");
    bimodal->Draw("SAME");

    // Legend
    TLegend* leg = new TLegend(0.12, 0.65, 0.45, 0.88);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(1);

    leg->AddEntry(h, "Angular distribution", "f");
    leg->AddEntry(g1,
        Form("Gaussian 1 (#mu=%.1f, #sigma=%.1f)", mu1, sigma1), "l");
    leg->AddEntry(g2,
        Form("Gaussian 2 (#mu=%.1f, #sigma=%.1f)", mu2, sigma2), "l");
    leg->AddEntry(bimodal, "Sum of two Gaussians", "l");

    leg->Draw("SAME");

    // Statistics box
    TPaveText* pt = new TPaveText(0.72, 0.75, 0.88, 0.88, "NDC");
    pt->SetFillColor(kWhite);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);

    pt->AddText(Form("Skewness: %.2f", skew));
    pt->AddText(Form("Kurtosis: %.2f", kurt));
    pt->Draw("SAME");

    // ================= SAVE AS PDF (VECTOR) =================
    c1->SaveAs(OUTPUT_PLOT);

    std::cout << "\nPlot successfully saved as: "
              << OUTPUT_PLOT << std::endl;
}
