void His_1D(
    const std::string& input_file,
    const std::string& output_pdf,
    const std::string& title_text,
    int bin_width_deg,
    double sigma_input
) {
    // -----------------------------
    //   1. Create histogram and fill
    // -----------------------------

    int n_bins = 180 / bin_width_deg;

    TH1F *hist = new TH1F(
        "hist",
        Form("%s %d degree;Angle;Particles", title_text.c_str(), bin_width_deg),
        n_bins, 0, 180
    );

    std::fstream file(input_file, std::ios::in);
    double value;
    while (file >> value) {
        hist->Fill(value);
    }
    file.close();

    // -----------------------------
    //   CALCULATION OF MEAN_Y AND SIGMA_Y
    // -----------------------------

    double manual_sigma = sigma_input;

    double sum_of_content = hist->Integral();
    double mean_y = sum_of_content / hist->GetNbinsX();
    double sigma_y = manual_sigma * TMath::Sqrt(mean_y);
    double sigma = sigma_y;

    // -----------------------------
    //   Apply global style
    // -----------------------------
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kLightTemperature);

    // -----------------------------
    //   Histogram style
    // -----------------------------
    hist->SetLineColor(kRed);
    hist->SetLineWidth(2);
    hist->SetFillColorAlpha(kRed, 0.3);
    hist->SetStats(0);

    // -----------------------------
    //   Canvas
    // -----------------------------
    TCanvas *c1 = new TCanvas("c1", "Histogram", 900, 600);

    hist->Draw("HIST");

    double x_min = hist->GetXaxis()->GetXmin();
    double x_max = hist->GetXaxis()->GetXmax();

    double y_max_needed = TMath::Max(hist->GetMaximum() * 1.05, mean_y + sigma * 1.1);
    hist->SetMaximum(y_max_needed);
    hist->SetMinimum(0);

    // -----------------------------
    //   Mean Y line
    // -----------------------------
    TLine *mean_y_line = new TLine(x_min, mean_y, x_max, mean_y);
    mean_y_line->SetLineColor(kBlue);
    mean_y_line->SetLineWidth(3);
    mean_y_line->SetLineStyle(2);
    mean_y_line->Draw("SAME");

    TLatex *mean_y_text = new TLatex();
    mean_y_text->SetTextSize(0.04);
    mean_y_text->SetTextColor(kBlue);
    mean_y_text->SetTextAlign(31);
    mean_y_text->DrawLatex(x_max * 0.95, mean_y * 1.05,
                           Form("Mean (Y): %.2f", mean_y));

    // -----------------------------
    //   Significance lines
    // -----------------------------
    double y_high = mean_y + sigma;
    double y_low  = TMath::Max(0.0, mean_y - sigma);

    TLine *sig_high = new TLine(x_min, y_high, x_max, y_high);
    sig_high->SetLineColor(kBlack);
    sig_high->SetLineWidth(2);
    sig_high->SetLineStyle(3);
    sig_high->Draw("SAME");

    TLine *sig_low = new TLine(x_min, y_low, x_max, y_low);
    sig_low->SetLineColor(kBlack);
    sig_low->SetLineWidth(2);
    sig_low->SetLineStyle(3);
    sig_low->Draw("SAME");

    TLatex *sigma_text = new TLatex();
    sigma_text->SetTextSize(0.035);
    sigma_text->SetTextColor(kBlack);
    sigma_text->SetTextAlign(12);

    sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01,
                          y_high * 0.98,
                          Form("+%.1f#sigma", manual_sigma));

    if (y_low > 0) {
        sigma_text->DrawLatex(x_min + (x_max - x_min) * 0.01,
                              y_low * 0.9,
                              Form("-%.1f#sigma", manual_sigma));
    }

    // -----------------------------
    //   SAVE PDF
    // -----------------------------
    std::string dir = output_pdf.substr(0, output_pdf.find_last_of('/'));
    gSystem->mkdir(dir.c_str(), kTRUE);
    c1->SaveAs(output_pdf.c_str());

    delete c1;
    delete hist;
}
