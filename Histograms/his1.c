void his1() {

    // -----------------------------
    //   1. Create histogram
    // -----------------------------
    TH1F *hist = new TH1F("hist", "Histogram;Angle;Particles", 180, 0, 180);

    std::fstream file("Output/test.txt", std::ios::in);
    double value;

    while (file >> value) {
        hist->Fill(value);
    }
    file.close();

    // -----------------------------
    //   2. Apply global style
    // -----------------------------
    gStyle->SetOptStat(0);        // No statistics box
    gStyle->SetPadGridX(kTRUE);   // Grid (X)
    gStyle->SetPadGridY(kTRUE);   // Grid (Y)
    gStyle->SetPadTickX(1);       // Ticks on all sides
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kLightTemperature);

    // -----------------------------
    //   3. Histogram style
    // -----------------------------
    hist->SetLineColor(kRed);
    hist->SetLineWidth(2);
    hist->SetFillColorAlpha(kRed, 0.3); // transparent fill
    hist->SetStats(0);

    // -----------------------------
    //   4. Canvas
    // -----------------------------
    TCanvas *c1 = new TCanvas("c1", "Histogram", 900, 600);
    //c1->SetLogy(); // optional, remove if not needed

    // -----------------------------
    //   5. Legend
    // -----------------------------
    //TLegend *leg = new TLegend(0.75, 0.80, 0.90, 0.92); 
    //leg->SetTextSize(0.03);                              // smaller font
    //leg->AddEntry(hist, "Entries", "lF");
    //leg->SetBorderSize(0);                               // without borders
    //leg->SetFillStyle(0);                                // transparent background


    // -----------------------------
    //   6. Draw
    // -----------------------------
    hist->Draw("HIST");
    //leg->Draw();
}