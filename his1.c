void his1(){
    TH1F *hist = new TH1F("hist", "Histogram", 180, 0, 180);

    fstream file;
    file.open("C:\\Users\\Kacperakis\\Downloads\\angles.txt", ios::in);

    double value;

    while(1){
        file >> value;
        hist->Fill(value);
        if(file.eof()) break;
    }

    file.close();

    hist->GetXaxis()->SetTitle("Angle");
    hist->GetYaxis()->SetTitle("Particles");

    TCanvas *c1 = new TCanvas();
    hist->Draw();
}