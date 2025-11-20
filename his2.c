#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>

void his2() {
    // Tworzymy okno
    TCanvas *c1 = new TCanvas();
    gStyle->SetPalette(kRainBow);

    // Histogram 2D dla wizualizacji
    TH2F *hist = new TH2F("hist", "Histogram", 100, -1, 1, 100, -1, 1);
    hist->SetStats(0);

    // Parametry prostej
    double theta = 45; // kąt w stopniach
    double a = TMath::Tan(theta * TMath::DegToRad()); // współczynnik kierunkowy
    int n = 100;     // liczba punktów
    double x_min = -0.9, x_max = 0.9;

    double dx = (x_max - x_min) / (n-1);
    double x = x_min;

    // Wypełnianie histogramu punktami na prostej
    for(int i = 0; i < n; i++) {
        double y = a * x;
        hist->Fill(x, y);
        x += dx; // liniowy wzrost X
    }

    // Opisy osi
    hist->GetXaxis()->SetTitle("x");
    hist->GetYaxis()->SetTitle("y");
    hist->GetZaxis()->SetTitle("z");

    hist->SetContour(1000);
    hist->Draw("colz");
}
