#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraph.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>

// This ROOT script has been made to visualise detection structure
// Still the programme is in development stage.

void his2() {

    // --- Plot style ---
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetNumberContours(600);

    const int RES = 1000;       // grid resolution
    const double Rmax = 1;      // drawing area radius

    TCanvas *c1 = new TCanvas("c1","PSF Test",1200,1200);

    TH2F *hist = new TH2F("hist",
                          "Coincidence 4;X;Y",
                          RES, -Rmax, Rmax,
                          RES, -Rmax, Rmax);
    hist->SetStats(0);

    // --- Load angles (0..180°) with filtering of invalid values ---
    std::vector<double> angles;
    std::ifstream fin("Output/test.txt");
    double ang;
    while (fin >> ang) {
        if (ang < 0.0 || ang > 180.0) continue;
        angles.push_back(ang);
    }
    fin.close();

    if (angles.empty()) angles.push_back(45.0); // fallback

    // --- Parameter: width of the uniform distribution interval (degrees) ---
    double spread_deg = 5;           // <-- set e.g. 1.0 or 9.0 here

    // Helper lambda: minimum angular difference in degrees (0..360)
    auto minimal_angle_diff_deg = [](double a_deg, double b_deg) {
        double diff = fabs(a_deg - b_deg);
        if (diff > 360.0) diff = fmod(diff, 360.0);
        if (diff > 180.0) diff = 360.0 - diff;
        return diff;
    };

    // --- For each angle, generate two directions: theta and theta + 180 ---
    for (double theta_deg : angles) {

        double thetas_deg[2] = { theta_deg, fmod(theta_deg + 180.0, 360.0) };

        for (int idir = 0; idir < 2; ++idir) {

            double theta_check_deg = thetas_deg[idir];

            // iterate over the grid of points
            for (int ix = 0; ix < RES; ++ix) {
                double x = -Rmax + 2.0 * Rmax * ix / (RES - 1);

                for (int iy = 0; iy < RES; ++iy) {
                    double y = -Rmax + 2.0 * Rmax * iy / (RES - 1);

                    double r = sqrt(x*x + y*y);
                    if (r > Rmax) continue; // out of area

                    // direction vector theta (in radians) for "in front of detector" test
                    double theta_rad = theta_check_deg * TMath::DegToRad();
                    double cosT = TMath::Cos(theta_rad);
                    double sinT = TMath::Sin(theta_rad);

                    // Point must lie on the "source side" for this direction
                    double dot = x * cosT + y * sinT;
                    if (dot <= 0) continue; // point behind the detector for this direction

                    // calculate φ in degrees in range [0,360)
                    double phi_rad = TMath::ATan2(y, x); // -pi..pi
                    double phi_deg = phi_rad * TMath::RadToDeg();
                    if (phi_deg < 0) phi_deg += 360.0;

                    // minimum angular difference (degrees)
                    double ddeg = minimal_angle_diff_deg(phi_deg, theta_check_deg);

                    // uniform distribution: if within distance <= spread_deg
                    if (ddeg <= spread_deg) {
                        hist->Fill(x, y, 1.0);
                    }
                }
            }
        } // end of idir loop
    } // end of angles loop

    // --- Z-AXIS NORMALIZATION (COLOR SCALE) to [0, 1] range ---
    
    // 1. Calculate the maximum value in the histogram
    double max_content = hist->GetMaximum();
    
    if (max_content > 0) {
        // 2. Scale all cells by 1/max_content
        hist->Scale(1.0 / max_content);
    }
    
    // 3. Set the color scale (Z) to [0, 1]
    hist->SetMaximum(1.0);
    // hist->SetMinimum(0.0); // Optional

    // --- Drawing ---
    hist->SetContour(800);
    hist->Draw("COLZ");

    // Detector: black dot
    TGraph *det = new TGraph();
    det->SetPoint(0, 0.0, 0.0);
    det->SetMarkerStyle(20);
    det->SetMarkerSize(3.0);
    det->SetMarkerColor(kBlack);
    det->Draw("P SAME");
}