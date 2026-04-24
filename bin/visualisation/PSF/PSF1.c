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
const char* INPUT_FILE = "/home/kacper/MuonDensityMapping/output/all/zussamen_10/coin_34f.txt";
const char* OUTPUT_NAME = "/home/kacper/MuonDensityMapping/figures/PSF_Reconstruction_1.pdf";

const int SMOOTH_RADIUS = 8;
const double SPREAD_DEG = 11.0;
const double MAX_DIST   = 70.0;
const double MAP_DIST   = 75.0;
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

void draw_decorations(double x, double y, double r) {
    TEllipse* e = new TEllipse(x, y, r, r);
    e->SetLineColor(kYellow);
    e->SetLineWidth(1);
    e->SetFillStyle(0);
    e->Draw("SAME");

    TGraph* g = new TGraph();
    g->SetPoint(0, x, y);
    g->SetMarkerStyle(29);
    g->SetMarkerSize(2.0);
    g->SetMarkerColor(kRed);
    g->Draw("P SAME");
}

void PSF1() {
    const double D1X = 0.0, D1Y = 0.0;
    const double span = MAP_DIST + 2.0;

    gStyle->SetOptStat(0);
    const Int_t NRGBs = 2;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 1.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 999);
    gStyle->SetNumberContours(999);

    auto dens1 = get_angular_density(INPUT_FILE, SPREAD_DEG);
    double max1 = *std::max_element(dens1.begin(), dens1.end()); 
    if (max1 < 1e-9) max1 = 1.0;

    TH2F* h_det1     = new TH2F("h_det1", "Raw Detector View;X [m];Y [m]",      RES, -span, span, RES, -span, span);
    TH2F* h_recon    = new TH2F("h_recon", "Matter Reconstruction;X [m];Y [m]", RES, -span, span, RES, -span, span);
    TH2F* h_coverage = new TH2F("h_cov",   "coverage",                          RES, -span, span, RES, -span, span);

    std::cout << "Generating maps..." << std::endl;

    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_det1->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y = h_det1->GetYaxis()->GetBinCenter(iy);
            double d = std::sqrt(x*x + y*y);
            
            if (d < MIN_DIST) continue;
            double phi = TMath::ATan2(x, y) * TMath::RadToDeg();
            if (phi < 0) phi += 360.0;
            int idx = (int)std::round(phi) % 360;

            h_det1->SetBinContent(ix, iy, dens1[idx]);
            h_recon->SetBinContent(ix, iy, dens1[idx] / max1);
            if (d <= MAX_DIST) h_coverage->SetBinContent(ix, iy, 1.0);
        }
    }

    std::cout << "Smoothing..." << std::endl;
    interpolate_single_detector_zones(h_recon, h_coverage, SMOOTH_RADIUS);

    apply_dynamic_scale_det(h_det1);
    apply_dynamic_scale_map(h_recon);

    TCanvas* c = new TCanvas("c", "Muon Tomography PSF", 1200, 600);
    c->Divide(2, 1, 0.01, 0.01);

    c->cd(1);
    gPad->SetRightMargin(0.15);
    h_det1->Draw("COLZ");
    draw_decorations(D1X, D1Y, MAX_DIST);

    c->cd(2);
    gPad->SetRightMargin(0.15);
    h_recon->Draw("COLZ");
    draw_decorations(D1X, D1Y, MAX_DIST);

    c->SaveAs(OUTPUT_NAME);
    std::cout << "Done. Output saved to: " << OUTPUT_NAME << std::endl;
}