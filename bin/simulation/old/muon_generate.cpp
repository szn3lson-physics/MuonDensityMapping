#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

int main() {
    ofstream fout("Output/test.txt");
    if (!fout.is_open()) {
        cerr << "Nie można otworzyć pliku!" << endl;
        return 1;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> unif(0.0, 1.0);

    // Wzmocniony kierunek (jaskinia)
    const double caveMin = 30.0 * M_PI / 180.0;
    const double caveMax = 50.0 * M_PI / 180.0;

    // JAK DUŻO WIĘCEJ? (prosty parametr)
    const double p_cave = 0.2;   // 50% mionów z kierunku jaskini

    int hours = 1000;  // więcej żeby zobaczyć statystykę

    for (int h = 0; h < hours; h++) {

        // --- LOSOWANIE AZYMUTU φ ---
        double phi;

        if (unif(gen) < p_cave) {
            // losujemy TYLKO z przedziału 30–50°
            phi = caveMin + unif(gen) * (caveMax - caveMin);
        } else {
            // losujemy z całego zakresu 0–180°
            phi = unif(gen) * M_PI;
        }

        // --- ZENIT θ ≈ poziomo ---
        uniform_real_distribution<> tdist(85.0, 95.0);
        double theta = tdist(gen) * M_PI / 180.0;

        // zapis do pliku
        int angle = phi * 180.0 / M_PI;
        fout << angle << endl;
    }

    fout.close();
    cout << "Zapisano miony do test.txt" << endl;
    return 0;
}

/*
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

int main() {
    ofstream fout("Output/test.txt");
    if (!fout.is_open()) {
        cerr << "Nie można otworzyć pliku!" << endl;
        return 1;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> unif(0.0, 1.0);

    // wymiary detektora
    const double detX = 10.0;
    const double detY = 20.0;

    // głębokość skały
    const double depth_m = 30.0;
    const double rock_density = 2.2;
    const double L = depth_m * 100.0 * rock_density;

    const double dEdx = 2e-3;

    // jaskinia 30–50° w zakresie 0–180°
    const double caveMin = 30.0 * M_PI / 180.0;
    const double caveMax = 50.0 * M_PI / 180.0;
    const double boost = 5.0;

    int hours = 1000;

    for (int h = 0; h < hours; h++) {

        // --- LOSOWANIE AZYMUTU φ TYLKO W PRZEDZIALE 0–180° ---
        double phi;
        while (true) {
            double ph = unif(gen) * M_PI;   // 0–π czyli 0°–180°

            double w = 10.0;
            if (ph >= caveMin && ph <= caveMax)
                w *= boost;

            if (unif(gen) < w) {
                phi = ph;
                break;
            }
        }

        // --- ZENIT θ ≈ 90° (poziomy mion) ---
        uniform_real_distribution<> tdist(85.0, 95.0);
        double theta = tdist(gen) * M_PI / 180.0;

        // kierunek
        double vx = sin(theta) * cos(phi);
        double vy = sin(theta) * sin(phi);
        double vz = cos(theta);

        // trafienie
        double x = (unif(gen) - 0.5) * detX;
        double y = (unif(gen) - 0.5) * detY;

        // energia
        uniform_real_distribution<> Edist(5.0, 300.0);
        double E0 = Edist(gen);

        double E_loss = dEdx * (L / fabs(cos(theta)));
        bool survives = (E0 > E_loss);

        int angle = phi * 180.0 / M_PI;

        //fout << "hour: " << h << "\n";
        //fout << "theta_deg: " << theta * 180.0 / M_PI << "\n";
        fout << angle << "\n";
        //fout << "vx: " << vx << "\n";
        //fout << "vy: " << vy << "\n";
        //fout << "vz: " << vz << "\n";
        //fout << "x_cm: " << x << "\n";
        //fout << "y_cm: " << y << "\n";
        //fout << "E0_GeV: " << E0 << "\n";
        //fout << "survives: " << (survives ? 1 : 0) << "\n";
        //fout << "---\n";
    }

    fout.close();
    cout << "Wygenerowano dane mionów dla " << hours << " godzin." << endl;
    return 0;
}
*/