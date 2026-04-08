#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TPad.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <fstream>
#include <iostream>
#include <string>

void plot_discovery() {
    // ================= PARAMETERS =================
    std::string data_file = "D:\\MuonDensityMapping\\output\\dane_10\\processed\\coin_10_34f.txt"; 
    std::string bkg_file  = "D:\\MuonDensityMapping\\bin\\simulation\\simulation_output\\MC_one_by_one\\30.00_metrow.txt";
    
    // Set our "golden middle" – 20 bins of 9 degrees each
    int n_bins = 20; 
    double x_min = 0;
    double x_max = 180;

    // ================= HISTOGRAM CREATION =================
    TH1F *h_data = new TH1F("h_data", "", n_bins, x_min, x_max);
    TH1F *h_bkg  = new TH1F("h_bkg",  "", n_bins, x_min, x_max);

    // To draw error bars for data, ROOT must know these are Poisson uncertainties
    h_data->Sumw2(); 

    // Load Data
    std::ifstream fd(data_file);
    double val;
    if (fd.is_open()) { while (fd >> val) h_data->Fill(val); fd.close(); }

    // Load Background (30m solid rock simulation)
    std::ifstream fb(bkg_file);
    if (fb.is_open()) { while (fb >> val) h_bkg->Fill(val); fb.close(); }

    // Compute "Signal" (Data - Background)
    TH1F *h_signal = (TH1F*)h_data->Clone("h_signal");
    h_signal->Add(h_bkg, -1.0); // Subtract background from data

    // ================= STYLING =================
    gStyle->SetOptStat(0);

    // Background style (classic particle physics: blue, semi-transparent fill)
    h_bkg->SetFillColorAlpha(kBlue - 7, 0.5);
    h_bkg->SetLineColor(kBlue + 1);
    h_bkg->SetLineWidth(2);

    // Data style (classic: black points with error bars)
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.2);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineColor(kBlack);
    h_data->SetLineWidth(2);

    // Signal style (lower panel)
    h_signal->SetMarkerStyle(21);
    h_signal->SetMarkerSize(1.2);
    h_signal->SetMarkerColor(kRed + 1);
    h_signal->SetLineColor(kRed + 1);
    h_signal->SetLineWidth(2);

    // ================= CANVAS AND PANELS (Key for Discovery-style plot) =================
    TCanvas *c1 = new TCanvas("c1", "Signal Analysis", 1000, 800);

    // Create upper panel (70% of height)
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02); // Remove bottom margin to merge with lower panel
    pad1->SetLeftMargin(0.12);
    pad1->SetGridx();
    pad1->SetGridy();
    pad1->Draw();
    pad1->cd();

    // Draw on upper panel
    h_bkg->SetMaximum(std::max(h_data->GetMaximum(), h_bkg->GetMaximum()) * 1.2);
    h_bkg->GetYaxis()->SetTitle("Muon counts / 9#circ");
    h_bkg->GetYaxis()->SetTitleSize(0.05);
    h_bkg->GetYaxis()->SetTitleOffset(1.0);
    h_bkg->GetYaxis()->SetLabelSize(0.04);
    h_bkg->Draw("HIST");       // Draw background
    h_data->Draw("E1 P SAME"); // Draw data as points with errors (E1 P) on top

    // Legend
    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->AddEntry(h_data, "Observed Data", "pe");
    leg->AddEntry(h_bkg, "Background Model (Solid rock 30m)", "f");
    leg->Draw();

    // Plot title
    TLatex *title = new TLatex();
    title->SetNDC();
    title->SetTextSize(0.06);
    title->SetTextFont(62);
    title->DrawLatex(0.15, 0.85, "Muography Anomaly Search");

    // Switch to main canvas to create lower panel
    c1->cd(); 
    
    // Create lower panel (30% of height)
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.12);
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    // Axis settings for lower panel (larger fonts due to smaller panel)
    h_signal->SetTitle("");
    h_signal->GetXaxis()->SetTitle("Zenith angle [#circ]");
    h_signal->GetXaxis()->SetTitleSize(0.12);
    h_signal->GetXaxis()->SetLabelSize(0.10);
    
    h_signal->GetYaxis()->SetTitle("Data - Background");
    h_signal->GetYaxis()->SetTitleSize(0.10);
    h_signal->GetYaxis()->SetTitleOffset(0.5);
    h_signal->GetYaxis()->SetLabelSize(0.08);
    h_signal->GetYaxis()->SetNdivisions(505);

    // Symmetric Y-axis for signal so zero is perfectly centered
    double max_sig = std::abs(h_signal->GetMaximum());
    double min_sig = std::abs(h_signal->GetMinimum());
    double limit = std::max(max_sig, min_sig) * 1.2;
    h_signal->SetMaximum(limit);
    h_signal->SetMinimum(-limit);

    h_signal->Draw("E1 P");

    // Zero line (perfect agreement between data and model)
    TLine *line0 = new TLine(x_min, 0, x_max, 0);
    line0->SetLineColor(kBlack);
    line0->SetLineStyle(2);
    line0->SetLineWidth(2);
    line0->Draw("SAME");

    c1->SaveAs("D:\\MuonDensityMapping\\wykres_higgs_style.png");
    c1->SaveAs("D:\\MuonDensityMapping\\wykres_higgs_style.pdf");
    std::cout << "Generated CERN-style discovery plot!" << std::endl;
}
