//====================================
// Kacper Dorszewski                 ||
// plot_muography_analysis.C         ||
//====================================

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
//  PARAMETERS (4 INPUT FILES)
// =============================================================================

const char* FILE_DET1_NOHOLE = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_1_without_hole.txt";
const char* FILE_DET2_NOHOLE = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_4_without_hole.txt";

const char* FILE_DET1_HOLE   = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_1_with_hole.txt";
const char* FILE_DET2_HOLE   = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_4_with_hole.txt";

const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_1_4_Comparison_BW_Inverted.pdf";

const double D1X_px = 915,   D1Y_px = 1116; //det_1
//const double D2X_px = 916,   D2Y_px = 1163; //det_2
//const double D2X_px = 948,   D2Y_px = 1030; //det_3
const double D2X_px = 963,   D2Y_px = 1247; //det_4

const double D1X =  0.0,  D1Y =  0.0;
const double D2X =  -(D1X_px - D2X_px) * 0.13333;
const double D2Y =   (D1Y_px - D2Y_px) * 0.13333;

const double SPREAD_DEG = 11.0;
const double MAX_DIST   = 48.0;
const double MAP_DIST   = 75.0;
const double MIN_DIST   = 0.2;
const int    RES        = 300;

// 0 = Hard Cutoff, 1 = Exponential, 2 = Gaussian
const int ATTENUATION_MODE = 0;

const double GAMMA_SINGLE   = 1.0;   
const double GAMMA_COMBINED = 1.0;   
const double BASELINE_PERCENTILE_SINGLE   = 0.0;
const double BASELINE_PERCENTILE_COMBINED = 0.0;

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

    double baseline = 1e18;
    for (int i = 1; i <= 180; i++)
        if (has_data[i]) baseline = std::min(baseline, raw_counts[i]);
    if (baseline > 1e17) baseline = 0.0;

    std::vector<double> delta(181, 0.0);
    for (int i = 1; i <= 180; i++)
        if (has_data[i]) delta[i] = std::max(0.0, raw_counts[i] - baseline);

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

void apply_gamma_scale(TH2F* h, double gamma, double baseline_frac) {
    double max_val = h->GetMaximum();
    if (max_val <= 1e-12) return;

    double cut = max_val * baseline_frac;
    double range = max_val - cut;
    if (range <= 1e-12) range = max_val;

    const double EMPTY_FILL = 1e-9; 

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);

            // Wartości ujemne lub poniżej cut ustawiamy na tło
            if (v <= cut || v <= 1e-12) {
                h->SetBinContent(i, j, EMPTY_FILL);
                continue;
            }

            double norm = (v - cut) / range;
            if (norm > 1.0) norm = 1.0;
            if (norm < 0.0) norm = 0.0;

            double v_out = std::pow(norm, gamma);
            h->SetBinContent(i, j, v_out);
        }
    }
}

void draw_single_decoration(double x, double y, double radius) {
    TEllipse* e = new TEllipse(x, y, radius, radius);
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
}

// =============================================================================
void PSF2_difference() {
    gStyle->SetOptStat(0);

    // --- ZMIANA: Odwrócona paleta Czarno-Biała (Czarne tło 0 -> Biały sygnał 1) ---
    const Int_t NRGBs = 2;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 1.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 999);
    gStyle->SetNumberContours(999);

    // 1. Pobranie danych wejściowych
    auto dens1_nh = get_angular_density(FILE_DET1_NOHOLE, SPREAD_DEG);
    auto dens2_nh = get_angular_density(FILE_DET2_NOHOLE, SPREAD_DEG);
    auto dens1_h  = get_angular_density(FILE_DET1_HOLE, SPREAD_DEG);
    auto dens2_h  = get_angular_density(FILE_DET2_HOLE, SPREAD_DEG);

    // Wyliczenie "sygnału różnicy" z poziomu kątów
    std::vector<double> dens1_diff(360, 0.0);
    std::vector<double> dens2_diff(360, 0.0);
    for (int i = 0; i < 360; ++i) {
        dens1_diff[i] = dens1_h[i] - dens1_nh[i];
        dens2_diff[i] = dens2_h[i] - dens2_nh[i];
    }

    double Xmin_s = -MAP_DIST - 2.0, Xmax_s = MAP_DIST + 2.0;
    double Ymin_s = -MAP_DIST - 2.0, Ymax_s = MAP_DIST + 2.0;

    double Xmin_c = std::min(D1X, D2X) - MAX_DIST - 2.0;
    double Xmax_c = std::max(D1X, D2X) + MAX_DIST + 2.0;
    double Ymin_c = std::min(D1Y, D2Y) - MAX_DIST - 2.0;
    double Ymax_c = std::max(D1Y, D2Y) + MAX_DIST + 2.0;

    // Utworzenie 9 histogramów
    TH2F* h1_nh = new TH2F("h1_nh", "Det 1 (Bez Dziury);X [m];Y [m]", RES, Xmin_s, Xmax_s, RES, Ymin_s, Ymax_s);
    TH2F* h2_nh = new TH2F("h2_nh", "Det 2 (Bez Dziury);X [m];Y [m]", RES, Xmin_s, Xmax_s, RES, Ymin_s, Ymax_s);
    TH2F* hc_nh = new TH2F("hc_nh", "Kombinacja (Bez Dziury);X [m];Y [m]", RES, Xmin_c, Xmax_c, RES, Ymin_c, Ymax_c);

    TH2F* h1_h = new TH2F("h1_h", "Det 1 (Z Dziura);X [m];Y [m]", RES, Xmin_s, Xmax_s, RES, Ymin_s, Ymax_s);
    TH2F* h2_h = new TH2F("h2_h", "Det 2 (Z Dziura);X [m];Y [m]", RES, Xmin_s, Xmax_s, RES, Ymin_s, Ymax_s);
    TH2F* hc_h = new TH2F("hc_h", "Kombinacja (Z Dziura);X [m];Y [m]", RES, Xmin_c, Xmax_c, RES, Ymin_c, Ymax_c);

    TH2F* h1_d = new TH2F("h1_d", "Roznica Det 1 (Z dziura - Bez dziury);X [m];Y [m]", RES, Xmin_s, Xmax_s, RES, Ymin_s, Ymax_s);
    TH2F* h2_d = new TH2F("h2_d", "Roznica Det 2 (Z dziura - Bez dziury);X [m];Y [m]", RES, Xmin_s, Xmax_s, RES, Ymin_s, Ymax_s);
    TH2F* hc_d = new TH2F("hc_d", "Rekonstrukcja z roznic;X [m];Y [m]", RES, Xmin_c, Xmax_c, RES, Ymin_c, Ymax_c);

    // --- Lambdy pomocnicze ---
    auto get_n_raw = [&](double dx, double dy, double dist, const std::vector<double>& dens) -> double {
        if (dist < MIN_DIST) return 0.0;
        double phi = TMath::ATan2(dx, dy) * TMath::RadToDeg();
        if (phi < 0) phi += 360.0;
        return dens[(int)std::round(phi) % 360];
    };

    auto get_n_attenuated = [&](double dx, double dy, double dist, const std::vector<double>& dens) -> double {
        if (dist < MIN_DIST) return 0.0;
        double phi = TMath::ATan2(dx, dy) * TMath::RadToDeg();
        if (phi < 0) phi += 360.0;
        double angular_flux = dens[(int)std::round(phi) % 360];
        double w = 0.0;
        if      (ATTENUATION_MODE == 0) w = (dist <= MAX_DIST) ? 1.0 : 0.0;
        else if (ATTENUATION_MODE == 1) w = std::exp(-dist / MAX_DIST);
        else                            w = std::exp(-(dist * dist) / (MAX_DIST * MAX_DIST));
        return angular_flux * w;
    };

    auto fill_single = [&](TH2F* h, const std::vector<double>& dens) {
        for (int ix = 1; ix <= RES; ++ix) {
            double x = h->GetXaxis()->GetBinCenter(ix);
            for (int iy = 1; iy <= RES; ++iy) {
                double y = h->GetYaxis()->GetBinCenter(iy);
                double dist = std::sqrt(x*x + y*y);
                h->SetBinContent(ix, iy, get_n_raw(x, y, dist, dens));
            }
        }
    };

    auto fill_combined = [&](TH2F* h, const std::vector<double>& dens1, const std::vector<double>& dens2) {
        for (int ix = 1; ix <= RES; ++ix) {
            double x = h->GetXaxis()->GetBinCenter(ix);
            for (int iy = 1; iy <= RES; ++iy) {
                double y   = h->GetYaxis()->GetBinCenter(iy);
                double dx1 = x - D1X, dy1 = y - D1Y, dist1 = std::sqrt(dx1*dx1 + dy1*dy1);
                double dx2 = x - D2X, dy2 = y - D2Y, dist2 = std::sqrt(dx2*dx2 + dy2*dy2);

                if (dist1 > MAX_DIST && dist2 > MAX_DIST) {
                    h->SetBinContent(ix, iy, 0.0);
                    continue;
                }
                if (dist1 >= MIN_DIST && dist2 >= MIN_DIST)
                    h->SetBinContent(ix, iy, get_n_attenuated(dx1, dy1, dist1, dens1) * get_n_attenuated(dx2, dy2, dist2, dens2));
                else
                    h->SetBinContent(ix, iy, 0.0);
            }
        }
    };

    // 2. Wypełnianie wszystkich histogramów
    fill_single(h1_nh, dens1_nh);
    fill_single(h2_nh, dens2_nh);
    fill_combined(hc_nh, dens1_nh, dens2_nh);

    fill_single(h1_h, dens1_h);
    fill_single(h2_h, dens2_h);
    fill_combined(hc_h, dens1_h, dens2_h);

    fill_single(h1_d, dens1_diff);
    fill_single(h2_d, dens2_diff);
    fill_combined(hc_d, dens1_diff, dens2_diff);

    // 3. Aplikacja kontrastu (Gamma)
    apply_gamma_scale(h1_nh, GAMMA_SINGLE,   BASELINE_PERCENTILE_SINGLE);
    apply_gamma_scale(h2_nh, GAMMA_SINGLE,   BASELINE_PERCENTILE_SINGLE);
    apply_gamma_scale(hc_nh, GAMMA_COMBINED, BASELINE_PERCENTILE_COMBINED);
    
    apply_gamma_scale(h1_h,  GAMMA_SINGLE,   BASELINE_PERCENTILE_SINGLE);
    apply_gamma_scale(h2_h,  GAMMA_SINGLE,   BASELINE_PERCENTILE_SINGLE);
    apply_gamma_scale(hc_h,  GAMMA_COMBINED, BASELINE_PERCENTILE_COMBINED);
    
    apply_gamma_scale(h1_d,  GAMMA_SINGLE,   BASELINE_PERCENTILE_SINGLE);
    apply_gamma_scale(h2_d,  GAMMA_SINGLE,   BASELINE_PERCENTILE_SINGLE);
    apply_gamma_scale(hc_d,  GAMMA_COMBINED, BASELINE_PERCENTILE_COMBINED);

    // Wymuszenie skali 0-1 dla kolorów
    TH2F* all_hists[9] = {h1_nh, h2_nh, hc_nh, h1_h, h2_h, hc_h, h1_d, h2_d, hc_d};
    for (int i=0; i<9; ++i) {
        all_hists[i]->SetMinimum(0.0);
        all_hists[i]->SetMaximum(1.0);
    }

    // --- Rysowanie siatki 3x3 ---
    TCanvas* c = new TCanvas("c", "Muon Tomography - PSF 3x3", 1800, 1500);
    c->Divide(3, 3, 0.005, 0.005);

    // Rząd 1: Bez Dziury
    c->cd(1); gPad->SetRightMargin(0.15); h1_nh->Draw("COLZ"); draw_single_decoration(0.0, 0.0, MAX_DIST);
    c->cd(2); gPad->SetRightMargin(0.15); h2_nh->Draw("COLZ"); draw_single_decoration(0.0, 0.0, MAX_DIST);
    c->cd(3); gPad->SetRightMargin(0.15); hc_nh->Draw("COLZ"); draw_single_decoration(D1X, D1Y, MAX_DIST); draw_single_decoration(D2X, D2Y, MAX_DIST);

    // Rząd 2: Z Dziurą
    c->cd(4); gPad->SetRightMargin(0.15); h1_h->Draw("COLZ"); draw_single_decoration(0.0, 0.0, MAX_DIST);
    c->cd(5); gPad->SetRightMargin(0.15); h2_h->Draw("COLZ"); draw_single_decoration(0.0, 0.0, MAX_DIST);
    c->cd(6); gPad->SetRightMargin(0.15); hc_h->Draw("COLZ"); draw_single_decoration(D1X, D1Y, MAX_DIST); draw_single_decoration(D2X, D2Y, MAX_DIST);

    // Rząd 3: Rekonstrukcja różnicowa
    c->cd(7); gPad->SetRightMargin(0.15); h1_d->Draw("COLZ"); draw_single_decoration(0.0, 0.0, MAX_DIST);
    c->cd(8); gPad->SetRightMargin(0.15); h2_d->Draw("COLZ"); draw_single_decoration(0.0, 0.0, MAX_DIST);
    c->cd(9); gPad->SetRightMargin(0.15); hc_d->Draw("COLZ"); draw_single_decoration(D1X, D1Y, MAX_DIST); draw_single_decoration(D2X, D2Y, MAX_DIST);

    c->SaveAs(OUTPUT_NAME);
    std::cout << "[SUCCESS] Saved 3x3 BW (Inverted) PSF comparison to: " << OUTPUT_NAME << std::endl;
}