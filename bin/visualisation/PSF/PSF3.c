#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

// =============================================================================
//  PARAMETERS FOR MANUAL ADJUSTMENT
// =============================================================================
//const char* FILE_DET1   = "/home/kacper/MuonDensityMapping/output/all/zussamen_10/coin_34f.txt";
//const char* FILE_DET2   = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_data/sim_2_root_events.txt";
//const char* FILE_DET3   = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_data/sim_3_root_events.txt";
const char* FILE_DET1   = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/1_hole_sim_root_C.txt";
const char* FILE_DET2   = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/2_hole_sim_root_C.txt";
const char* FILE_DET3   = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/3_hole_sim_root_C.txt";
const char* OUTPUT_FILE = "/home/kacper/MuonDensityMapping/figures/PSF_3_Detectors_Combined.pdf";

const double D1X_px = 915,  D1Y_px = 1116;
const double D2X_px = 948,  D2Y_px = 1030;
const double D3X_px = 916,  D3Y_px = 1163;

const double D1X =  0.0,  D1Y =  0.0;
const double D2X =  -(D1X_px-D2X_px)*0.13333, D2Y = (D1Y_px-D2Y_px)*0.13333;
const double D3X =  -(D1X_px-D3X_px)*0.13333, D3Y = (D1Y_px-D3Y_px)*0.13333;

//const double D1X =  0.00,  D1Y =  0.00;
//const double D2X =  4.40, D2Y = 11.47;
//const double D3X =  0.13, D3Y = -6.27;

const int SMOOTH_RADIUS = 8;
const double SPREAD_DEG = 11.0;
const double MAX_DIST   = 40.0; // Właściwy zasięg detektora (żółty okrąg)
const double MAP_DIST   = 75.0; // Maksymalny promień nanoszony na mapę
const double MIN_DIST   = 0.2;
const int RES           = 500;
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
    double ang;
    int count = 0;
    while (fin >> ang) {
        ang = fmod(ang, 180.0);
        if (ang < 0) ang += 180.0;
        double angles_to_fill[2] = { ang, fmod(ang + 180.0, 360.0) };
        for (int k = 0; k < 2; ++k) {
            double center = angles_to_fill[k];
            int half = (int)std::ceil(spread_deg);
            for (int i = -half; i <= half; ++i) {
                int idx = ((int)std::round(center) + i % 360 + 360) % 360;
                double diff = GetAngularDiff(idx, center);
                if (diff <= spread_deg) {
                    density[idx] += 1.0 - (diff / spread_deg);
                }
            }
        }
        count++;
    }
    fin.close();
    std::cout << "Loaded " << count << " events from: " << filename << std::endl;
    return density;
}

void apply_dynamic_scale_det(TH2F* h) {
    double max_val = h->GetMaximum();
    double min_val = max_val;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);
            if (v > 1e-9 && v < min_val) min_val = v;
        }
    if (min_val >= max_val) { min_val = 0.0; if (max_val == 0.0) max_val = 1e-6; }
    
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            if (h->GetBinContent(i, j) < min_val) h->SetBinContent(i, j, min_val);
        }
    }
    h->SetMinimum(min_val - 1e-6); 
    h->SetMaximum(max_val);
}

void apply_dynamic_scale_map(TH2F* h) {
    double max_val = h->GetMaximum();
    double sum = 0.0; int count = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);
            if (v > 1e-6) { sum += v; count++; }
        }
    double mean_val = (count > 0) ? (sum / count) : 0.0;
    double min_val = mean_val * 0.90;
    if (min_val >= max_val) { min_val = 0.0; if (max_val == 0.0) max_val = 1e-6; }
    
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            if (h->GetBinContent(i, j) < min_val) h->SetBinContent(i, j, min_val);
        }
    }
    h->SetMinimum(min_val - 1e-6);
    h->SetMaximum(max_val);
}

void apply_extreme_contrast(TH2F* h) {
    double max_val = h->GetMaximum();
    double sum = 0.0; int count = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);
            if (v > 1e-6) { sum += v; count++; }
        }
    }
    double mean_val = (count > 0) ? (sum / count) : 0.0;
    double min_val = mean_val * 1.20; 
    
    if (min_val >= max_val) { min_val = 0.0; if (max_val == 0.0) max_val = 1e-6; }
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            if (h->GetBinContent(i, j) < min_val) h->SetBinContent(i, j, min_val);
        }
    }
    h->SetMinimum(min_val - 1e-6);
    h->SetMaximum(max_val);
}

void interpolate_single_detector_zones(TH2F* h_signal, TH2F* h_coverage, int radius) {
    int nx = h_signal->GetNbinsX();
    int ny = h_signal->GetNbinsY();
    std::vector<std::vector<double>> orig(nx+1, std::vector<double>(ny+1, 0.0));
    std::vector<std::vector<double>> cov (nx+1, std::vector<double>(ny+1, 0.0));
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j) {
            orig[i][j] = h_signal->GetBinContent(i, j);
            cov [i][j] = h_coverage->GetBinContent(i, j);
        }
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            if (cov[i][j] != 1) continue;
            double sum_w = 0.0, sum_wv = 0.0;
            for (int di = -radius; di <= radius; ++di) {
                for (int dj = -radius; dj <= radius; ++dj) {
                    int ni = i+di, nj = j+dj;
                    if (ni < 1 || ni > nx || nj < 1 || nj > ny) continue;
                    if (cov[ni][nj] < 1) continue; 
                    double sigma = radius / 2.0;
                    double w = std::exp(-(di*di + dj*dj) / (2.0*sigma*sigma));
                    sum_w  += w;
                    sum_wv += w * orig[ni][nj];
                }
            }
            if (sum_w > 1e-9) h_signal->SetBinContent(i, j, sum_wv / sum_w);
        }
    }
}

void draw_decorations(bool d1, bool d2, bool d3) {
    auto draw_one = [](double x, double y) {
        TEllipse* e = new TEllipse(x, y, MAX_DIST, MAX_DIST);
        e->SetLineColor(kYellow); e->SetLineWidth(1); e->SetFillStyle(0); e->Draw("SAME");
        TGraph* g = new TGraph(); g->SetPoint(0, x, y);
        g->SetMarkerStyle(29); g->SetMarkerSize(2.0); g->SetMarkerColor(kRed); g->Draw("P SAME");
    };
    if (d1) draw_one(D1X, D1Y);
    if (d2) draw_one(D2X, D2Y);
    if (d3) draw_one(D3X, D3Y);
}

// =============================================================================
void PSF3() {
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
    auto dens3 = get_angular_density(FILE_DET3, SPREAD_DEG);

    double max1 = *std::max_element(dens1.begin(), dens1.end()); if (max1 < 1e-9) max1 = 1.0;
    double max2 = *std::max_element(dens2.begin(), dens2.end()); if (max2 < 1e-9) max2 = 1.0;
    double max3 = *std::max_element(dens3.begin(), dens3.end()); if (max3 < 1e-9) max3 = 1.0;

    double Xmin = std::min({D1X, D2X, D3X}) - MAP_DIST - 2.0;
    double Xmax = std::max({D1X, D2X, D3X}) + MAP_DIST + 2.0;
    double Ymin = std::min({D1Y, D2Y, D3Y}) - MAP_DIST - 2.0;
    double Ymax = std::max({D1Y, D2Y, D3Y}) + MAP_DIST + 2.0;

    TH2F* h_det1    = new TH2F("h_det1", "Detector 1;X [m];Y [m]", RES, Xmin, Xmax, RES, Ymin, Ymax);
    TH2F* h_det2    = new TH2F("h_det2", "Detector 2;X [m];Y [m]", RES, Xmin, Xmax, RES, Ymin, Ymax);
    TH2F* h_det3    = new TH2F("h_det3", "Detector 3;X [m];Y [m]", RES, Xmin, Xmax, RES, Ymin, Ymax);
    
    TH2F* h_overlap  = new TH2F("h_overlap", "Superposition;X [m];Y [m]",   RES, Xmin, Xmax, RES, Ymin, Ymax);
    TH2F* h_corr     = new TH2F("h_corr",    "Triangulation;X [m];Y [m]",    RES, Xmin, Xmax, RES, Ymin, Ymax);
    TH2F* h_recon    = new TH2F("h_recon",   "Reconstruction;X [m];Y [m]",   RES, Xmin, Xmax, RES, Ymin, Ymax);
    TH2F* h_coverage = new TH2F("h_cov",     "coverage",                     RES, Xmin, Xmax, RES, Ymin, Ymax);

    std::cout << "Processing grid..." << std::endl;

    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_det1->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y = h_det1->GetYaxis()->GetBinCenter(iy);

            double dx1 = x - D1X, dy1 = y - D1Y, dist1 = std::sqrt(dx1*dx1 + dy1*dy1);
            double dx2 = x - D2X, dy2 = y - D2Y, dist2 = std::sqrt(dx2*dx2 + dy2*dy2);
            double dx3 = x - D3X, dy3 = y - D3Y, dist3 = std::sqrt(dx3*dx3 + dy3*dy3);

            auto get_n = [&](double dx, double dy, double dist, const std::vector<double>& dens, double mmax) {
                if (dist < MIN_DIST) return 0.0;
                double phi = TMath::ATan2(dx, dy) * TMath::RadToDeg();
                if (phi < 0) phi += 360.0;
                return dens[(int)std::round(phi) % 360] / mmax;
            };

            double n1 = get_n(dx1, dy1, dist1, dens1, max1);
            double n2 = get_n(dx2, dy2, dist2, dens2, max2);
            double n3 = get_n(dx3, dy3, dist3, dens3, max3);

            // Zmiana: Rysowanie detektorów na pełnym zasięgu osi, ignorując MAX_DIST
            h_det1->SetBinContent(ix, iy, n1);
            h_det2->SetBinContent(ix, iy, n2);
            h_det3->SetBinContent(ix, iy, n3);

            bool a1 = (dist1 >= MIN_DIST), a2 = (dist2 >= MIN_DIST), a3 = (dist3 >= MIN_DIST);

            // Superposition
            h_overlap->SetBinContent(ix, iy, (a1?n1:0) + (a2?n2:0) + (a3?n3:0));

            // Triangulation
            if (a1 && a2 && a3) h_corr->SetBinContent(ix, iy, n1 * n2 * n3);
            else h_corr->SetBinContent(ix, iy, 0.0);

            // Reconstruction
            double rv = 0.0;
            int n_alive = (int)a1 + (int)a2 + (int)a3;
            if (n_alive == 3)      rv = std::cbrt(n1 * n2 * n3);
            else if (n_alive == 2) {
                if      (a1&&a2) rv = std::sqrt(n1*n2);
                else if (a1&&a3) rv = std::sqrt(n1*n3);
                else             rv = std::sqrt(n2*n3);
            } else if (n_alive == 1) {
                if      (a1) rv = n1;
                else if (a2) rv = n2;
                else         rv = n3;
            }
            h_recon->SetBinContent(ix, iy, rv);
            
            bool in1 = (dist1 >= MIN_DIST && dist1 <= MAX_DIST);
            bool in2 = (dist2 >= MIN_DIST && dist2 <= MAX_DIST);
            bool in3 = (dist3 >= MIN_DIST && dist3 <= MAX_DIST);
            h_coverage->SetBinContent(ix, iy, (double)((int)in1 + (int)in2 + (int)in3));
        }
    }

    std::cout << "Interpolation..." << std::endl;
    interpolate_single_detector_zones(h_recon, h_coverage, SMOOTH_RADIUS);

    apply_dynamic_scale_det(h_det1);
    apply_dynamic_scale_det(h_det2);
    apply_dynamic_scale_det(h_det3);
    
    apply_extreme_contrast(h_overlap); 
    apply_dynamic_scale_map(h_corr);
    apply_dynamic_scale_map(h_recon);

    // ════════════════════════════════════════════════════════════════════════
    //  COMBINED CANVAS (3x2 Grid)
    // ════════════════════════════════════════════════════════════════════════
    TCanvas* c = new TCanvas("c", "Muon Tomography PSF 3-Detectors", 1600, 1000);
    c->Divide(3, 2, 0.005, 0.005);

    c->cd(1); gPad->SetRightMargin(0.15); h_det1->Draw("COLZ"); draw_decorations(true, false, false);
    c->cd(2); gPad->SetRightMargin(0.15); h_det2->Draw("COLZ"); draw_decorations(false, true, false);
    c->cd(3); gPad->SetRightMargin(0.15); h_det3->Draw("COLZ"); draw_decorations(false, false, true);
    
    c->cd(4); gPad->SetRightMargin(0.15); h_overlap->Draw("COLZ"); draw_decorations(true, true, true);
    c->cd(5); gPad->SetRightMargin(0.15); h_corr->Draw("COLZ");    draw_decorations(true, true, true);
    c->cd(6); gPad->SetRightMargin(0.15); h_recon->Draw("COLZ");   draw_decorations(true, true, true);

    c->SaveAs(OUTPUT_FILE);
    std::cout << "Finished. File " << OUTPUT_FILE << " saved." << std::endl;
}