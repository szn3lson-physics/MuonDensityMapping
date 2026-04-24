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
//  PARAMETERS FOR MANUAL ADJUSTMENT
// =============================================================================
//const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_test_det_1.txt";
//const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_test_det_2.txt";

const char* FILE_DET1 = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_series/series_7_8_9_10/coin_34f.txt";
const char* FILE_DET2 = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/det_test_det_2.txt";

const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_2_Detectors_Gauss_Combined.pdf";

const double D1X_px = 1000,  D1Y_px = 850;
const double D2X_px = 925,   D2Y_px = 1000;

const double D1X =  0.0,  D1Y =  0.0;
const double D2X =  -(D1X_px-D2X_px)*0.13333, D2Y = (D1Y_px-D2Y_px)*0.13333;

// Rozmycie kątowe (akceptacja detektora). 11.0 oznacza stożek +/- 11 stopni
const double SPREAD_DEG = 11.0; 

const double MAX_DIST   = 50.0; // Promień R0 dla wygaszania
const double MAP_DIST   = 75.0; // Rozmiar siatki mapy
const double MIN_DIST   = 0.2;
const int RES           = 300; // Zwiększona rozdzielczość

// =======================================================
// --- PRZEŁĄCZNIK MODELU WYGASZANIA PROMIENI NA MAPIE ---
// 0 = Hard Cutoff (Natychmiastowe obcięcie po MAX_DIST)
// 1 = Exponential (Zanik eksponencjalny e^(-r/MAX_DIST))
// 2 = Gaussian (Zanik gaussowski e^(-(r/MAX_DIST)^2))
// =======================================================
const int ATTENUATION_MODE = 1; 

// =============================================================================

double GetAngularDiff(double a, double b) {
    double d = std::abs(a - b);
    if (d > 180.0) d = 360.0 - d;
    return d;
}

std::vector<double> get_angular_density(const char* filename, double spread_deg) {
    std::vector<double> density(360, 0.0);
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "[ERROR] Cannot open file: " << filename << std::endl;
        return density;
    }

    std::string line;
    
    while (std::getline(fin, line)) {
        if (line.empty() || line.find("Angle") != std::string::npos) continue;

        double ang = 0.0;
        double weight = 1.0;

        if (line.find(';') != std::string::npos) {
            std::stringstream ss(line);
            std::string item;
            std::getline(ss, item, ';'); ang = std::stod(item);
            std::getline(ss, item, ';'); 
            std::getline(ss, item, ';'); weight = std::stod(item); 
        } else {
            ang = std::stod(line);
        }

        ang = fmod(ang, 180.0);
        if (ang < 0) ang += 180.0;

        double angles_to_fill[2] = { ang, fmod(ang + 180.0, 360.0) };

        for (int k = 0; k < 2; ++k) {
            double center = angles_to_fill[k];
            int half = (int)std::ceil(spread_deg);
            
            for (int i = -half; i <= half; ++i) {
                int idx = ((int)std::round(center) + i + 360) % 360;
                double diff = GetAngularDiff(idx, center);
                
                if (diff <= spread_deg) {
                    double triangle_weight = 1.0 - (diff / spread_deg);
                    density[idx] += weight * triangle_weight;
                }
            }
        }
    }
    fin.close();

    double sum_weights = 0.0;
    for (double w : density) sum_weights += w;
    if (sum_weights > 0) {
        for (double& w : density) w /= sum_weights;
    }

    std::cout << "Loaded events/lines from: " << filename << " (Normalized to integral 1.0)" << std::endl;
    return density;
}

// =============================================================================
// ZMODYFIKOWANE SKALOWANIE: SPÓJNE DLA MODE 0 I PŁYNNE DLA MODE 1/2
// =============================================================================
void apply_clean_scale(TH2F* h, int mode) {
    const int nx = h->GetNbinsX();
    const int ny = h->GetNbinsY();
    
    double max_val = h->GetMaximum();
    if (max_val <= 0.0) return;

    if (mode == 0) {
        // --- 1. Klasyczne podejście dla HARD CUTOFF ---
        // Szukamy absolutnego tła skały, pomijając czyste matematyczne zera.
        double min_val = max_val;
        for (int i = 1; i <= nx; ++i) {
            for (int j = 1; j <= ny; ++j) {
                double v = h->GetBinContent(i, j);
                if (v > 1e-12 && v < min_val) min_val = v;
            }
        }
        if (min_val >= max_val) min_val = 0.0;
        
        // Zamazujemy puste kąty i dalekie zera czernią (min_val).
        for (int i = 1; i <= nx; ++i) {
            for (int j = 1; j <= ny; ++j) {
                if (h->GetBinContent(i, j) < min_val) {
                    h->SetBinContent(i, j, min_val);
                }
            }
        }
        h->SetMinimum(min_val);
        h->SetMaximum(max_val);
    } 
    else {
        // --- 2. Odwracanie wygaszania z wykorzystaniem zbalansowanego Logarytmu (dla EXP i GAUSS) ---
        // Obcinamy ultra-zera matematyczne (wszystko poniżej np. e^-5 * max_val). 
        double noise_floor = max_val * 1e-4; 
        
        // Offset dodany przed logarytmowaniem ratuje nas przed nieskończonością.
        // Gwarantuje płynny gradient dla słabych wiązek zamiast ostrego "ucinanego ogona".
        double offset = noise_floor; 
        
        double new_max = std::log10(max_val + offset);
        double new_min = std::log10(noise_floor + offset);
        
        for (int i = 1; i <= nx; ++i) {
            for (int j = 1; j <= ny; ++j) {
                double v = h->GetBinContent(i, j);
                
                // Jeśli sygnał jest mniejszy niż noise_floor, traktujemy go jako czerń tła.
                if (v < noise_floor) {
                    h->SetBinContent(i, j, new_min);
                } else {
                    double log_v = std::log10(v + offset);
                    h->SetBinContent(i, j, log_v);
                }
            }
        }
        
        h->SetMinimum(new_min);
        h->SetMaximum(new_max);
    }
}

void draw_decorations(bool d1, bool d2) {
    auto draw_one = [](double x, double y) {
        TEllipse* e = new TEllipse(x, y, MAX_DIST, MAX_DIST);
        e->SetLineColor(kYellow); e->SetLineWidth(1); e->SetLineStyle(2); e->SetFillStyle(0); e->Draw("SAME");
        TGraph* g = new TGraph(); g->SetPoint(0, x, y);
        g->SetMarkerStyle(29); g->SetMarkerSize(2.0); g->SetMarkerColor(kRed); g->Draw("P SAME");
    };
    if (d1) draw_one(D1X, D1Y);
    if (d2) draw_one(D2X, D2Y);
}

// =============================================================================
void PSF2_prof() {
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

    std::string mode_str = (ATTENUATION_MODE == 0) ? "Hard Cutoff" : (ATTENUATION_MODE == 1) ? "Exponential" : "Gaussian";
    std::cout << "Processing grid with " << mode_str << " attenuation..." << std::endl;

    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_det1->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y = h_det1->GetYaxis()->GetBinCenter(iy);

            double dx1 = x - D1X, dy1 = y - D1Y, dist1 = std::sqrt(dx1*dx1 + dy1*dy1);
            double dx2 = x - D2X, dy2 = y - D2Y, dist2 = std::sqrt(dx2*dx2 + dy2*dy2);

            auto get_n = [&](double dx, double dy, double dist, const std::vector<double>& dens) {
                if (dist < MIN_DIST) return 0.0;
                
                double phi = TMath::ATan2(dx, dy) * TMath::RadToDeg();
                if (phi < 0) phi += 360.0;
                
                double angular_flux = dens[(int)std::round(phi) % 360];

                double spatial_weight = 0.0;
                if (ATTENUATION_MODE == 0) {
                    spatial_weight = (dist <= MAX_DIST) ? 1.0 : 0.0;
                } else if (ATTENUATION_MODE == 1) {
                    spatial_weight = std::exp(-dist / MAX_DIST);
                } else if (ATTENUATION_MODE == 2) {
                    spatial_weight = std::exp(-(dist * dist) / (MAX_DIST * MAX_DIST));
                }

                return angular_flux * spatial_weight;
            };

            double n1 = get_n(dx1, dy1, dist1, dens1);
            double n2 = get_n(dx2, dy2, dist2, dens2);

            h_det1->SetBinContent(ix, iy, n1);
            h_det2->SetBinContent(ix, iy, n2);

            if (dist1 >= MIN_DIST && dist2 >= MIN_DIST) {
                h_combined->SetBinContent(ix, iy, n1 * n2);
            } else {
                h_combined->SetBinContent(ix, iy, 0.0);
            }
        }
    }

    // APLIKACJA NOWEJ FUNKCJI SKALUJĄCEJ
    apply_clean_scale(h_det1, ATTENUATION_MODE);
    apply_clean_scale(h_det2, ATTENUATION_MODE);
    apply_clean_scale(h_combined, ATTENUATION_MODE);

    TCanvas* c = new TCanvas("c", "Muon Tomography PSF 2-Detectors", 1800, 600);
    c->Divide(3, 1, 0.005, 0.005);

    c->cd(1); gPad->SetRightMargin(0.15); h_det1->Draw("COLZ"); draw_decorations(true, false);
    c->cd(2); gPad->SetRightMargin(0.15); h_det2->Draw("COLZ"); draw_decorations(false, true);
    c->cd(3); gPad->SetRightMargin(0.15); h_combined->Draw("COLZ"); draw_decorations(true, true);

    c->SaveAs(OUTPUT_NAME);
    std::cout << "Finished. File " << OUTPUT_NAME << " saved." << std::endl;
}