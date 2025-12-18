void his1() {
    // -----------------------------
    //   1. Create histogram and fill
    // -----------------------------
    TH1F *hist = new TH1F("hist", "Histogram 3 Coincidences 15 degree;Angle;Particles", 12, 0, 180);

    std::fstream file("Output/dane_4/processed/coin_4_3.txt", std::ios::in);
    double value;
    while (file >> value) {
        hist->Fill(value);
    }
    file.close();

    // -----------------------------
    //   CALCULATION OF MEAN_Y AND SIGMA_Y
    // -----------------------------
    
    // Set manually the sigma of the histogram (e.g., 5 sigmas)
    double manual_sigma = 1.0;

    // Calculate mean Y value (number of particles/counts)
    // Create a temporary histogram to easily retrieve Y statistics
    TH1F *h_temp = (TH1F*)hist->Clone("h_temp");
    
    // In ROOT, GetMean() on a TH1F returns the mean of the X-axis.
    // To get the mean of the counts (Y-axis), we must average the content of all bins.
    double sum_of_entries = hist->GetEntries();
    double sum_of_content = hist->Integral(); // Sum of content of all bins
    
    // Mean bin content (Mean_Y)
    double mean_y = sum_of_content / hist->GetNbinsX(); 
    
    // Sigma as the square root of mean bin content (Poisson Statistics)
    double sigma_y = manual_sigma * TMath::Sqrt(mean_y);
    
    // Significance level
    double sigma = manual_sigma * sigma_y;

    delete h_temp; // Remove the temporary histogram
    
    // -----------------------------
    //   2. Apply global style
    // -----------------------------
    gStyle->SetOptStat(0);       // No statistics box
    gStyle->SetPadGridX(kTRUE);  // Grid (X)
    gStyle->SetPadGridY(kTRUE);  // Grid (Y)
    gStyle->SetPadTickX(1);      // Ticks on all sides
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kLightTemperature);

    // -----------------------------
    //   3. Histogram style
    // -----------------------------
    hist->SetLineColor(kRed);
    hist->SetLineWidth(2);
    hist->SetFillColorAlpha(kRed, 0.3); // Transparent fill
    hist->SetStats(0);

    // -----------------------------
    //   4. Canvas
    // -----------------------------
    TCanvas *c1 = new TCanvas("c1", "Histogram", 900, 600);

    // -----------------------------
    //   6. Draw (first to establish the axes)
    // -----------------------------
    hist->Draw("HIST");
    
    // Get minimum and maximum X values for the full width of the plot
    double x_min = hist->GetXaxis()->GetXmin();
    double x_max = hist->GetXaxis()->GetXmax();
    
    // Ensure the Y range is sufficient to display the lines
    double y_max_needed = TMath::Max(hist->GetMaximum() * 1.05, mean_y + sigma * 1.1);
    hist->SetMaximum(y_max_needed);
    hist->SetMinimum(0); // Histograms usually start at 0


    // -----------------------------
    //   7. Mean Y line (Baseline)
    // -----------------------------
    // Draw a horizontal line at Mean_Y height
    TLine *mean_y_line = new TLine(x_min, mean_y, x_max, mean_y);
    mean_y_line->SetLineColor(kBlue);
    mean_y_line->SetLineWidth(3);
    mean_y_line->SetLineStyle(2); // Dotted line
    mean_y_line->Draw("SAME");
    
    // Text showing the mean Y value
    TLatex *mean_y_text = new TLatex();
    mean_y_text->SetTextSize(0.04);
    mean_y_text->SetTextColor(kBlue);
    mean_y_text->SetTextAlign(31); // Align Right, Bottom
    mean_y_text->DrawLatex(x_max * 0.95, mean_y * 1.05, Form("Mean (Y): %.2f", mean_y));


    // -----------------------------
    //   8. Significance lines (+- 5*sigma Y)
    // -----------------------------
    
    // Upper significance level (+5*sigma)
    double y_high = mean_y + sigma;
    TLine *sig_high = new TLine(x_min, y_high, x_max, y_high);
    sig_high->SetLineColor(kBlack);
    sig_high->SetLineWidth(2);
    sig_high->SetLineStyle(3); // Dashed line
    sig_high->Draw("SAME");

    // Lower significance level (-5*sigma)
    double y_low = mean_y - sigma;
    // Limit the bottom line to 0, as counts cannot be negative
    if (y_low < 0) y_low = 0; 
    TLine *sig_low = new TLine(x_min, y_low, x_max, y_low);
    sig_low->SetLineColor(kBlack);
    sig_low->SetLineWidth(2);
    sig_low->SetLineStyle(3); // Dashed line
    sig_low->Draw("SAME");

    // Text for the sigma range
    TLatex *sigma_text = new TLatex();
    sigma_text->SetTextSize(0.035);
    sigma_text->SetTextColor(kBlack);
    sigma_text->SetTextAlign(12); // Align Left, Center
    
    // +5sigma label
    sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01, y_high * 0.98, Form("+%.1f#sigma", manual_sigma));
    // -5sigma label
    if (y_low > 0) {
        sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01, y_low * 0.9, Form("-%.1f#sigma", manual_sigma));
    }
}