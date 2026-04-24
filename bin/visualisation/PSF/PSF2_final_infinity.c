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
//const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_rot_series/series_400.txt";
//const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_rot_series/series_900.txt";
//const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/glory_hole_pos_1.txt";

//const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_1_without_hole.txt";
//const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/MC_det_1_without_hole.txt";

const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_1_with_hole.txt";

//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_2_without_hole.txt";
//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_3_without_hole.txt";
//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_4_without_hole.txt";

//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_2_with_hole.txt";
//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_3_with_hole.txt";
const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_4_with_hole.txt";

//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/MC_det_2_without_hole.txt";
//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/MC_det_3_without_hole.txt";
//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/MC_det_4_without_hole.txt";

//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/400/MC/pos_2_400_bkg_instant.txt";
//const char* FILE_DET2 ="/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/900/MC/pos_2_900_bkg_instant.txt";
//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/glory_hole_pos_2.txt";

//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/real_400_PSF_2_Detectors_Gauss_Combined.pdf";
//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/real_900_PSF_2_Detectors_Gauss_Combined.pdf";
//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/sim_PSF_2_Detectors_Gauss_Combined.pdf";

//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_1_2.pdf";
//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_1_3.pdf";
//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_1_4.pdf";

//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_MC_1_2.pdf";
//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_MC_1_3.pdf";
//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_MC_1_4.pdf";

//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_1_2_hole.pdf";
//const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_1_3_hole.pdf";
const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_1_4_hole.pdf";

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
const double ZOOM_FACTOR = 2.0;

// 0 = Hard Cutoff, 1 = Exponential, 2 = Gaussian
const int ATTENUATION_MODE = 0;

// =============================================================================
//  KONTRAST – GŁÓWNY PARAMETR DO STROJENIA
// =============================================================================
//
//  GAMMA_SINGLE  – wykładnik dla map pojedynczych detektorów (Wykresy 1 i 2).
//  GAMMA_COMBINED – wykładnik dla mapy zrekonstruowanej (Wykres 3).
//
//  Jak działa:
//    - Każda wartość piksela jest normalizowana do zakresu [0, 1],
//      a następnie podniesiona do potęgi GAMMA: v_out = v_norm^GAMMA
//
//    GAMMA = 1.0  ->  brak modyfikacji (liniowa skala, wyjście jak poprzednio)
//    GAMMA < 1.0  ->  rozciąga ciemne obszary (uwypukla słabe sygnały)
//                     np. 0.3 daje bardzo mocny efekt "rozjaśnienia tła"
//    GAMMA > 1.0  ->  ściska ciemne obszary (uwypukla TYLKO silny sygnał)
//                     np. 3.0 "wypala" słabe kierunki i zostają tylko dominujące
//
//  ZALECANE WARTOŚCI DO PRÓB:
//    Wykres 1/2 (detektory): zacznij od GAMMA = 0.3 lub 0.4
//    Wykres 3 (combined):    zacznij od GAMMA = 0.25 lub 0.5
//    Jeśli mapa nadal płaska: zmniejsz GAMMA (np. 0.15)
//    Jeśli mapa zbyt przepalona: zwiększ GAMMA (np. 0.6)
//
//  BASELINE_PERCENTILE – procent minimalnej wartości do odjęcia przed gamma.
//    0.0  = brak odejmowania tła (wartości absolutne)
//    0.02 = odejmij 2% maksimum jako baseline (zalecane dla combined)
//    0.05 = agresywne odjęcie tła
//
// =============================================================================
const double GAMMA_SINGLE   = 1;   // dla Wykresu 1 i 2
const double GAMMA_COMBINED = 1;   // dla Wykresu 3
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

    // Wyznacz linię bazową (minimum po wszystkich kątach) i odejmij
    double baseline = 1e18;
    for (int i = 1; i <= 180; i++)
        if (has_data[i]) baseline = std::min(baseline, raw_counts[i]);
    if (baseline > 1e17) baseline = 0.0;

    std::vector<double> delta(181, 0.0);
    for (int i = 1; i <= 180; i++)
        if (has_data[i]) delta[i] = std::max(0.0, raw_counts[i] - baseline);

    // Nałożenie Gaussowskiej PSF kątowej
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

    // Normalizacja do sumy = 1
    double sum_d = 0.0;
    for (double v : density) sum_d += v;
    if (sum_d > 0)
        for (double& v : density) v /= sum_d;

    return density;
}

// =============================================================================
//  apply_gamma_scale:
//    Skaluje histogram do [0,1] z uwzględnieniem odcięcia tła (baseline_frac),
//    a następnie aplikuje korekcję gamma: v_out = v_norm ^ gamma.
//    Wartość 0 po odcięciu jest ustawiana na minimalną niezerową wartość
//    żeby ROOT nie renderował jej jako "brak danych".
// =============================================================================
void apply_gamma_scale(TH2F* h, double gamma, double baseline_frac) {
    double max_val = h->GetMaximum();
    if (max_val <= 1e-12) return;

    double cut = max_val * baseline_frac;
    double range = max_val - cut;
    if (range <= 1e-12) range = max_val;

    const double EMPTY_FILL = 1e-9; // wartość dla "czarnych" pikseli

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);

            if (v <= cut || v <= 1e-12) {
                h->SetBinContent(i, j, EMPTY_FILL);
                continue;
            }

            // Normalizuj do [0,1] po odjęciu baseline
            double norm = (v - cut) / range;
            if (norm > 1.0) norm = 1.0;
            if (norm < 0.0) norm = 0.0;

            // Aplikuj gamma: v^gamma
            // gamma < 1 rozciąga ciemne obszary (więcej kontrastu w słabym sygnale)
            // gamma > 1 ściska ciemne obszary (wyostrza dominujące kierunki)
            double v_out = std::pow(norm, gamma);

            h->SetBinContent(i, j, v_out);
        }
    }

    h->SetMinimum(0.0);
    h->SetMaximum(1.0);
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
void PSF2_final_infinity() {
    gStyle->SetOptStat(0);

    // Skala szarości
    const Int_t NRGBs = 2;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 1.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 999);
    gStyle->SetNumberContours(999);

    auto dens1 = get_angular_density(FILE_DET1, SPREAD_DEG);
    auto dens2 = get_angular_density(FILE_DET2, SPREAD_DEG);

    // Przestrzeń dla wykresów 1 i 2 (wycentrowana na 0,0)
    double Xmin_s = -MAP_DIST - 2.0, Xmax_s = MAP_DIST + 2.0;
    double Ymin_s = -MAP_DIST - 2.0, Ymax_s = MAP_DIST + 2.0;

    // Przestrzeń dla wykresu 3 — granice dopasowane do kółek o promieniu MAX_DIST
    double Xmin_c = std::min(D1X, D2X) - MAX_DIST - 2.0;
    double Xmax_c = std::max(D1X, D2X) + MAX_DIST + 2.0;
    double Ymin_c = std::min(D1Y, D2Y) - MAX_DIST - 2.0;
    double Ymax_c = std::max(D1Y, D2Y) + MAX_DIST + 2.0;

    TH2F* h_det1     = new TH2F("h_det1",     "Teleskop 1;X [m];Y [m]",                      RES, Xmin_s, Xmax_s, RES, Ymin_s, Ymax_s);
    TH2F* h_det2     = new TH2F("h_det2",     "Teleskop 2;X [m];Y [m]",                      RES, Xmin_s, Xmax_s, RES, Ymin_s, Ymax_s);
    TH2F* h_combined = new TH2F("h_combined", "Zrekonstruowana Gestosc (T1 * T2);X [m];Y [m]", RES, Xmin_c, Xmax_c, RES, Ymin_c, Ymax_c);

    // --- Lambdy do obliczania wartości piksela ---

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

    // --- Wypełnianie histogramów ---

    // Wykres 1: detektor 1, środek (0,0)
    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_det1->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y = h_det1->GetYaxis()->GetBinCenter(iy);
            double dist = std::sqrt(x*x + y*y);
            h_det1->SetBinContent(ix, iy, get_n_raw(x, y, dist, dens1));
        }
    }

    // Wykres 2: detektor 2, środek (0,0)
    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_det2->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y = h_det2->GetYaxis()->GetBinCenter(iy);
            double dist = std::sqrt(x*x + y*y);
            h_det2->SetBinContent(ix, iy, get_n_raw(x, y, dist, dens2));
        }
    }

    // Wykres 3: combined, prawdziwe pozycje D1 i D2
    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_combined->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y   = h_combined->GetYaxis()->GetBinCenter(iy);
            double dx1 = x - D1X, dy1 = y - D1Y, dist1 = std::sqrt(dx1*dx1 + dy1*dy1);
            double dx2 = x - D2X, dy2 = y - D2Y, dist2 = std::sqrt(dx2*dx2 + dy2*dy2);

            // Zeruj piksele poza zasięgiem obu kółek
            if (dist1 > MAX_DIST && dist2 > MAX_DIST) {
                h_combined->SetBinContent(ix, iy, 0.0);
                continue;
            }

            if (dist1 >= MIN_DIST && dist2 >= MIN_DIST)
                h_combined->SetBinContent(ix, iy, get_n_attenuated(dx1, dy1, dist1, dens1)
                                                * get_n_attenuated(dx2, dy2, dist2, dens2));
            else
                h_combined->SetBinContent(ix, iy, 0.0);
        }
    }

    // --- Aplikacja korekcji gamma ---
    apply_gamma_scale(h_det1,     GAMMA_SINGLE,   BASELINE_PERCENTILE_SINGLE);
    apply_gamma_scale(h_det2,     GAMMA_SINGLE,   BASELINE_PERCENTILE_SINGLE);
    apply_gamma_scale(h_combined, GAMMA_COMBINED, BASELINE_PERCENTILE_COMBINED);

    // --- Rysowanie ---
    TCanvas* c = new TCanvas("c", "Muon Tomography – Contrast Enhanced", 1800, 600);
    c->Divide(3, 1, 0.005, 0.005);

    c->cd(1); gPad->SetRightMargin(0.15);
    h_det1->Draw("COLZ");
    draw_single_decoration(0.0, 0.0, MAX_DIST);

    c->cd(2); gPad->SetRightMargin(0.15);
    h_det2->Draw("COLZ");
    draw_single_decoration(0.0, 0.0, MAX_DIST);

    c->cd(3); gPad->SetRightMargin(0.15);
    h_combined->Draw("COLZ");
    draw_single_decoration(D1X, D1Y, MAX_DIST);
    draw_single_decoration(D2X, D2Y, MAX_DIST);

    c->SaveAs(OUTPUT_NAME);
    std::cout << "Saved: " << OUTPUT_NAME << std::endl;
    std::cout << "Gamma single:   " << GAMMA_SINGLE   << std::endl;
    std::cout << "Gamma combined: " << GAMMA_COMBINED << std::endl;
}