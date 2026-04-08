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
const double THRESHOLD = 0.5;
const int SMOOTH_RADIUS = 8;
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
    if (min_val >= max_val) {
        min_val = 0.0;
        if (max_val == 0.0) max_val = 1e-6;
    }
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        for (int j = 1; j <= h->GetNbinsY(); ++j)
            if (h->GetBinContent(i, j) < min_val)
                h->SetBinContent(i, j, min_val);
    h->SetMinimum(min_val);
    h->SetMaximum(max_val);
}

void apply_dynamic_scale_map(TH2F* h) {
    double max_val = h->GetMaximum();
    double sum = 0.0;
    int count = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);
            if (v > 1e-6) { sum += v; count++; }
        }
    double mean_val = (count > 0) ? (sum / count) : 0.0;
    double min_val = mean_val * 0.90;
    if (min_val >= max_val) {
        min_val = 0.0;
        if (max_val == 0.0) max_val = 1e-6;
    }
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        for (int j = 1; j <= h->GetNbinsY(); ++j)
            if (h->GetBinContent(i, j) < min_val)
                h->SetBinContent(i, j, min_val);
    h->SetMinimum(min_val);
    h->SetMaximum(max_val);
}

// --- NEW SPECIFIC FOR SUPERPOSITION PLOT (Sharp contrast) ---
void apply_extreme_contrast(TH2F* h) {
    double max_val = h->GetMaximum();
    double sum = 0.0;
    int count = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double v = h->GetBinContent(i, j);
            if (v > 1e-6) { sum += v; count++; }
        }
    }
    double mean_val = (count > 0) ? (sum / count) : 0.0;
    
    // Raise the cutoff threshold up to 1.2 x Background Mean 
    double min_val = mean_val * 1.20; 
    
    if (min_val >= max_val) {
        min_val = 0.0;
        if (max_val == 0.0) max_val = 1e-6;
    }
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            if (h->GetBinContent(i, j) < min_val) {
                h->SetBinContent(i, j, min_val);
            }
        }
    }
    h->SetMinimum(min_val);
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
            if (sum_w > 1e-9)
                h_signal->SetBinContent(i, j, sum_wv / sum_w);
        }
    }
}

void draw_yellow_circle(double x, double y, double r) {
    TEllipse* e = new TEllipse(x, y, r, r);
    e->SetLineColor(kYellow);
    e->SetLineWidth(2);
    e->SetLineStyle(1);
    e->SetFillStyle(0);
    e->Draw("SAME");
}

void draw_det(double x, double y) {
    TGraph* g = new TGraph();
    g->SetPoint(0, x, y);
    g->SetMarkerStyle(29);
    g->SetMarkerSize(2.5);
    g->SetMarkerColor(kRed);
    g->Draw("P SAME");
}

// =============================================================================
void PSF2() {

    const double D1X =  0.0,  D1Y =  0.0;
    const double D2X =  6.66, D2Y = 11.1;

    const double max_dist   = 70.0;
    const double map_dist   = 75.0;
    const double spread_deg = 11.0;
    const double min_dist   = 0.2;

    const double Xs1 = -map_dist - 2.0, Xs2 = map_dist + 2.0;
    const double Ys1 = -map_dist - 2.0, Ys2 = map_dist + 2.0;

    double cx = (std::min({D1X,D2X}) + std::max({D1X,D2X})) / 2.0;
    double cy = (std::min({D1Y,D2Y}) + std::max({D1Y,D2Y})) / 2.0;
    double half_span = map_dist + std::max({
        std::abs(D1X - cx), std::abs(D2X - cx),
        std::abs(D1Y - cy), std::abs(D2Y - cy)
    }) + 2.0;
    const double Xc1 = cx - half_span, Xc2 = cx + half_span;
    const double Yc1 = cy - half_span, Yc2 = cy + half_span;

    const int RES = 500;

    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.10);

    const Int_t NRGBs = 2;
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 1.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 999);
    gStyle->SetNumberContours(999);

    auto dens1 = get_angular_density("D:\\MuonDensityMapping\\sim_2_root_events.txt", spread_deg);
    auto dens2 = get_angular_density("D:\\MuonDensityMapping\\sim_3_root_events.txt", spread_deg);

    double max1 = *std::max_element(dens1.begin(), dens1.end()); if (max1 < 1e-9) max1 = 1.0;
    double max2 = *std::max_element(dens2.begin(), dens2.end()); if (max2 < 1e-9) max2 = 1.0;

    TH2F* h_det1 = new TH2F("h_det1","Detector 1;X [m];Y [m]",                  RES,Xs1,Xs2, RES,Ys1,Ys2);
    TH2F* h_det2 = new TH2F("h_det2","Detector 2;X [m];Y [m]",                  RES,Xs1,Xs2, RES,Ys1,Ys2);
    
    TH2F* h_overlap  = new TH2F("h_overlap", "Superposition 2 det.;X [m];Y [m]",  RES,Xc1,Xc2, RES,Yc1,Yc2);
    TH2F* h_corr     = new TH2F("h_corr",    "Triangulation (full coverage);X [m];Y [m]", RES,Xc1,Xc2, RES,Yc1,Yc2);
    TH2F* h_recon    = new TH2F("h_recon",   "Matter Reconstruction;X [m];Y [m]", RES,Xc1,Xc2, RES,Yc1,Yc2);
    TH2F* h_binary   = new TH2F("h_binary",  "Binary Map;X [m];Y [m]",          RES,Xc1,Xc2, RES,Yc1,Yc2);
    TH2F* h_coverage = new TH2F("h_cov",     "coverage",                                  RES,Xc1,Xc2, RES,Yc1,Yc2);

    std::cout << "Processing grid..." << std::endl;

    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_det1->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y = h_det1->GetYaxis()->GetBinCenter(iy);
            double dist = std::sqrt(x*x + y*y);
            if (dist < min_dist) continue;
            double phi_deg = TMath::ATan2(x, y) * TMath::RadToDeg();
            if (phi_deg < 0) phi_deg += 360.0;
            int idx = (int)std::round(phi_deg) % 360;
            h_det1->SetBinContent(ix, iy, dens1[idx]);
            h_det2->SetBinContent(ix, iy, dens2[idx]);
        }
    }

    // ── Main calculations ────────────────────────────────────────────────────
    for (int ix = 1; ix <= RES; ++ix) {
        double x = h_corr->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= RES; ++iy) {
            double y = h_corr->GetYaxis()->GetBinCenter(iy);

            double dx1=x-D1X, dy1=y-D1Y, dist1=std::sqrt(dx1*dx1+dy1*dy1);
            double dx2=x-D2X, dy2=y-D2Y, dist2=std::sqrt(dx2*dx2+dy2*dy2);

            auto get_n = [&](double dx, double dy, double dist,
                             const std::vector<double>& dens, double mmax) -> double {
                if (dist < min_dist) return 0.0;
                double phi = TMath::ATan2(dx, dy) * TMath::RadToDeg();
                if (phi < 0) phi += 360.0;
                return dens[(int)std::round(phi) % 360] / mmax;
            };

            double n1 = get_n(dx1,dy1,dist1,dens1,max1);
            double n2 = get_n(dx2,dy2,dist2,dens2,max2);

            bool in1 = (dist1 >= min_dist && dist1 <= max_dist);
            bool in2 = (dist2 >= min_dist && dist2 <= max_dist);

            // SUPERPOSITION (Full map, no clipping by yellow circles)
            double o1_full = (dist1 >= min_dist) ? n1 : 0.0;
            double o2_full = (dist2 >= min_dist) ? n2 : 0.0;
            h_overlap->SetBinContent(ix, iy, o1_full + o2_full);

            // Triangulation (on the whole map)
            if ((dist1 >= min_dist) && (dist2 >= min_dist))
                h_corr->SetBinContent(ix, iy, n1 * n2);
            else
                h_corr->SetBinContent(ix, iy, 0.0);

            // Reconstruction
            double rv = 0.0;
            bool a1 = (dist1 >= min_dist);
            bool a2 = (dist2 >= min_dist);

            int n_alive = (int)a1 + (int)a2;
            if (n_alive == 2)      rv = std::sqrt(n1 * n2);
            else if (n_alive == 1) {
                if      (a1) rv = n1;
                else         rv = n2;
            }
            h_recon->SetBinContent(ix, iy, rv);
            
            // Coverage for binarization (uses max_dist)
            h_coverage->SetBinContent(ix, iy, (double)((int)in1+(int)in2));
        }
    }

    std::cout << "Interpolation..." << std::endl;
    interpolate_single_detector_zones(h_recon, h_coverage, SMOOTH_RADIUS);

    apply_dynamic_scale_map(h_recon);
    double r_min = h_recon->GetMinimum();
    double r_max = h_recon->GetMaximum();
    if (r_max - r_min < 1e-9) r_max = r_min + 1.0;

    for (int i = 1; i <= RES; ++i) {
        for (int j = 1; j <= RES; ++j) {
            double cv = h_coverage->GetBinContent(i, j);
            if (cv < 1e-6) {
                h_binary->SetBinContent(i, j, 0.0);
            } else {
                double v = h_recon->GetBinContent(i, j);
                double norm_v = (v - r_min) / (r_max - r_min);
                h_binary->SetBinContent(i, j, norm_v >= THRESHOLD ? 1.0 : 0.0);
            }
        }
    }
    h_binary->SetMinimum(0.0);
    h_binary->SetMaximum(1.0);

    apply_dynamic_scale_det(h_det1);
    apply_dynamic_scale_det(h_det2);
    apply_extreme_contrast(h_overlap); 
    apply_dynamic_scale_map(h_corr);

    // ════════════════════════════════════════════════════════════════════════
    //  CANVAS 1: 1x3 - Detectors and Superposition
    // ════════════════════════════════════════════════════════════════════════
    TCanvas* c1 = new TCanvas("c1", "PSF 2 detectors", 1680, 560);
    c1->Divide(3, 1, 0.002, 0.002);

    c1->cd(1); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.13); gPad->SetBottomMargin(0.11); gPad->SetTopMargin(0.09);
    h_det1->Draw("COLZ");
    draw_yellow_circle(0.0, 0.0, max_dist);
    draw_det(0.0, 0.0);

    c1->cd(2); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.13); gPad->SetBottomMargin(0.11); gPad->SetTopMargin(0.09);
    h_det2->Draw("COLZ");
    draw_yellow_circle(0.0, 0.0, max_dist);
    draw_det(0.0, 0.0);

    c1->cd(3); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.13); gPad->SetBottomMargin(0.11); gPad->SetTopMargin(0.09);
    h_overlap->Draw("COLZ");
    draw_yellow_circle(D1X, D1Y, max_dist);
    draw_yellow_circle(D2X, D2Y, max_dist);
    draw_det(D1X, D1Y);
    draw_det(D2X, D2Y);

    c1->SaveAs("PSF_2_detectors.pdf");

    // ════════════════════════════════════════════════════════════════════════
    //  CANVAS 2: Triangulation and Binary Map
    // ════════════════════════════════════════════════════════════════════════
    TCanvas* c2 = new TCanvas("c2", "Matter Reconstruction", 1120, 560);
    c2->Divide(2, 1, 0.002, 0.002);

    c2->cd(1); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.13); gPad->SetBottomMargin(0.11); gPad->SetTopMargin(0.09);
    h_corr->Draw("COLZ");
    draw_yellow_circle(D1X, D1Y, max_dist);
    draw_yellow_circle(D2X, D2Y, max_dist);
    draw_det(D1X, D1Y);
    draw_det(D2X, D2Y);

    c2->cd(2); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.13); gPad->SetBottomMargin(0.11); gPad->SetTopMargin(0.09);
    h_binary->Draw("COLZ");
    draw_yellow_circle(D1X, D1Y, max_dist);
    draw_yellow_circle(D2X, D2Y, max_dist);
    draw_det(D1X, D1Y);
    draw_det(D2X, D2Y);

    TLatex* label = new TLatex();
    label->SetNDC();
    label->SetTextSize(0.04);
    label->SetTextColor(kYellow);
    label->DrawLatex(0.15, 0.93, Form("Threshold: %.0f%%  |  white=void  black=rock", THRESHOLD*100.0));

    c2->SaveAs("PSF_2_reconstruction.pdf");

    std::cout << "Finished. Files: PSF_2_detectors.pdf and PSF_2_reconstruction.pdf saved." << std::endl;
}