#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TLatex.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

// =============================================================================
//  PARAMETERS
// =============================================================================
//const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_test_det_1.txt";
const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_1_without_hole.txt";
//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_test_det_2.txt";
const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_2_without_hole.txt";

const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/jala.pdf";

const double D1X_px = 1000,  D1Y_px = 850;
const double D2X_px = 925,   D2Y_px = 1000;

const double D1X =  0.0,  D1Y =  0.0;
const double D2X =  -(D1X_px-D2X_px)*0.13333, D2Y = (D1Y_px-D2Y_px)*0.13333;

const double SPREAD_DEG = 11.0;

const double MAX_DIST   = 50.0;
const double MAP_DIST   = 75.0;
const double MIN_DIST   = 0.2;
const int RES           = 300;

// 0 = Hard Cutoff
// 1 = Exponential
// 2 = Gaussian
const int ATTENUATION_MODE = 0;

// =============================================================================

double GetAngularDiff(double a, double b) {
    double d = std::abs(a - b);
    if (d > 180.0) d = 360.0 - d;
    return d;
}

std::vector<double> get_angular_density(const char* filename, double spread_deg) {
    std::vector<double> raw_counts(181, 0.0); 
    std::vector<bool>   has_data(181, false);

    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "[ERROR] Cannot open: " << filename << std::endl;
        return std::vector<double>(360, 0.0);
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line.find("Angle") != std::string::npos) continue;
        if (line.find(';') == std::string::npos) continue;

        std::stringstream ss(line);
        std::string item;
        std::getline(ss, item, ';'); int ang = std::stoi(item);
        std::getline(ss, item, ';');  
        std::getline(ss, item, ';'); double cnt = std::stod(item);

        if (ang >= 0 && ang <= 180 && cnt > 0) {
            raw_counts[ang] = cnt;
            has_data[ang] = true;
        }
    }
    fin.close();

    // Wyciąganie "Baseline" - odcinanie tła litych skał
    double baseline = 1e18;
    for (int i = 1; i <= 180; i++)
        if (has_data[i]) baseline = std::min(baseline, raw_counts[i]);
    if (baseline > 1e17) baseline = 0.0;

    std::vector<double> delta(181, 0.0);
    for (int i = 1; i <= 180; i++) {
        if (has_data[i]) {
            delta[i] = std::max(0.0, raw_counts[i] - baseline);
        }
    }

    std::vector<double> density(360, 0.0);
    double sigma = spread_deg / 2.355;  

    for (int ang_half = 1; ang_half <= 180; ang_half++) {
        if (!has_data[ang_half] || delta[ang_half] < 1e-12) continue;

        double dirs[2] = { (double)ang_half, fmod(ang_half + 180.0, 360.0) };

        for (int k = 0; k < 2; k++) {
            double center = dirs[k];

            for (int idx = 0; idx < 360; idx++) {
                double diff = GetAngularDiff(idx, center);
                if (diff > 3.0 * spread_deg) continue;  

                double gauss = std::exp(-0.5 * (diff / sigma) * (diff / sigma));
                density[idx] += delta[ang_half] * gauss;
            }
        }
    }

    double sum_d = 0.0;
    for (double v : density) sum_d += v;
    if (sum_d > 0)
        for (double& v : density) v /= sum_d;

    return density;
}

// =============================================================================
// KULOODPORNE SKALOWANIE: TRIK Z 1e-9 ZABEZPIECZAJĄCY PRZED BIAŁYM TŁEM
// =============================================================================
void apply_clean_scale(TH2F* h, int mode) {
    double max_val = h->GetMaximum();
    if (max_val <= 1e-12) return;

    if (mode == 0) {
        // Twarde odcięcie: wszystko puste dostaje 1e-9 (zamiast 0.0), żeby ROOT
        // wymusił rysowanie na czarno i zlikwidował białą przezroczystość.
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            for (int j = 1; j <= h->GetNbinsY(); ++j) {
                if (h->GetBinContent(i, j) <= 1e-12) {
                    h->SetBinContent(i, j, 1e-9); 
                }
            }
        }
        h->SetMinimum(0.0); 
        return;
    }

    // --- Dla trybu Exp i Gauss ---
    // Odcięcie szumu tła na poziomie 2% maximum
    double cut = max_val * 0.02; 
    double range = max_val - cut;
    
    if (range <= 1e-12) return;

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);

            if (v <= cut) {
                // Trik: mikroskopijna wartość 1e-9 zamiast 0.0. ROOT to narysuje!
                h->SetBinContent(i, j, 1e-9); 
            } else {
                double norm = (v - cut) / range;
                if (norm > 1.0) norm = 1.0; 
                if (norm < 0.0) norm = 0.0; 

                // Podbijamy sygnał pierwiastkiem
                h->SetBinContent(i, j, std::pow(norm, 0.5));
            }
        }
    }

    // Skala ustawiona na twardo, a najniższy piksel ma wartość niezerową
    h->SetMinimum(0.0);
    h->SetMaximum(1.0);
}

// =============================================================================

void draw_decorations(bool d1, bool d2) {
    auto draw_one = [](double x, double y) {
        TEllipse* e = new TEllipse(x, y, MAX_DIST, MAX_DIST);
        e->SetLineColor(kYellow);
        e->SetLineWidth(1);
        e->SetLineStyle(2);
        e->SetFillStyle(0);
        e->Draw("SAME");

        TGraph* g = new TGraph();
        g->SetPoint(0, x, y);
        g->SetMarkerStyle(29);
        g->SetMarkerSize(2.0);
        g->SetMarkerColor(kRed);
        g->Draw("P SAME");
    };

    if (d1) draw_one(D1X, D1Y);
    if (d2) draw_one(D2X, D2Y);
}

// =============================================================================

void PSF2_final()
{
    gStyle->SetOptStat(0);

    const Int_t NRGBs = 2;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 1.00 };

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 999);
    gStyle->SetNumberContours(999);

    auto dens1 = get_angular_density(FILE_DET1, SPREAD_DEG);
    auto dens2 = get_angular_density(FILE_DET2, SPREAD_DEG);

    double Xmin = std::min(D1X, D2X) - MAP_DIST - 2.0;
    double Xmax = std::max(D1X, D2X) + MAP_DIST + 2.0;
    double Ymin = std::min(D1Y, D2Y) - MAP_DIST - 2.0;
    double Ymax = std::max(D1Y, D2Y) + MAP_DIST + 2.0;

    TH2F* h_det1     = new TH2F("h_det1", "Teleskop 1;X [m];Y [m]", RES, Xmin, Xmax, RES, Ymin, Ymax);
    TH2F* h_det2     = new TH2F("h_det2", "Teleskop 2;X [m];Y [m]", RES, Xmin, Xmax, RES, Ymin, Ymax);
    TH2F* h_combined = new TH2F("h_combined", "Zrekonstruowana Gestosc (T1 * T2);X [m];Y [m]", RES, Xmin, Xmax, RES, Ymin, Ymax);

    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_det1->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y = h_det1->GetYaxis()->GetBinCenter(iy);

            double dx1 = x - D1X, dy1 = y - D1Y;
            double dx2 = x - D2X, dy2 = y - D2Y;

            double dist1 = std::sqrt(dx1*dx1 + dy1*dy1);
            double dist2 = std::sqrt(dx2*dx2 + dy2*dy2);

            auto get_n = [&](double dx, double dy, double dist, const std::vector<double>& dens) {

                if (dist < MIN_DIST) return 0.0;

                double phi = TMath::ATan2(dx, dy) * TMath::RadToDeg();
                if (phi < 0) phi += 360.0;

                double angular_flux = dens[(int)std::round(phi) % 360];

                double spatial_weight = 0.0;

                if (ATTENUATION_MODE == 0)
                    spatial_weight = (dist <= MAX_DIST) ? 1.0 : 0.0;
                else if (ATTENUATION_MODE == 1)
                    spatial_weight = std::exp(-dist / MAX_DIST);
                else
                    spatial_weight = std::exp(-(dist * dist) / (MAX_DIST * MAX_DIST));

                return angular_flux * spatial_weight;
            };

            double n1 = get_n(dx1, dy1, dist1, dens1);
            double n2 = get_n(dx2, dy2, dist2, dens2);

            h_det1->SetBinContent(ix, iy, n1);
            h_det2->SetBinContent(ix, iy, n2);

            if (dist1 >= MIN_DIST && dist2 >= MIN_DIST)
                h_combined->SetBinContent(ix, iy, n1 * n2);
            else
                h_combined->SetBinContent(ix, iy, 0.0);
        }
    }

    apply_clean_scale(h_det1, ATTENUATION_MODE);
    apply_clean_scale(h_det2, ATTENUATION_MODE);
    apply_clean_scale(h_combined, ATTENUATION_MODE);

    TCanvas* c = new TCanvas("c", "Muon Tomography", 1800, 600);
    c->Divide(3, 1, 0.005, 0.005);

    c->cd(1); gPad->SetRightMargin(0.15); h_det1->Draw("COLZ"); draw_decorations(true, false);
    c->cd(2); gPad->SetRightMargin(0.15); h_det2->Draw("COLZ"); draw_decorations(false, true);
    c->cd(3); gPad->SetRightMargin(0.15); h_combined->Draw("COLZ"); draw_decorations(true, true);

    c->SaveAs(OUTPUT_NAME);
}