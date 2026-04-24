#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>
#include <iomanip>

using namespace std;
namespace fs = std::filesystem;

struct DataPoint {
    double time = 0.0;
    long long counts = 0;
};

// Funkcja pomocnicza do przetwarzania pojedynczego pliku
bool processFile(const string& filepath, DataPoint dataArray[]) {
    ifstream file(filepath);
    if (!file.is_open()) {
        cerr << "[BŁĄD] Nie można otworzyć pliku: " << filepath << endl;
        return false;
    }

    string line;
    // Pomiń nagłówek
    getline(file, line);

    while (getline(file, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string item;
        vector<string> tokens;

        while (getline(ss, item, ';')) {
            tokens.push_back(item);
        }

        if (tokens.size() >= 3) {
            try {
                int angle = stoi(tokens[0]);
                double time = stod(tokens[1]);
                long long counts = stoll(tokens[2]);

                if (angle >= 0 && angle <= 180) {
                    dataArray[angle].time += time;
                    dataArray[angle].counts += counts;
                }
            } catch (const exception& e) {
                // Pomiń błędy konwersji w uszkodzonych liniach
                continue;
            }
        }
    }
    file.close();
    return true;
}

int main() {
    // Ścieżki wejściowe
    string file1 = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_series/series_1_2_3_4_5_7/coin_34f.txt";
    string file2 = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_rotations_900/6_series/51-54_rotations.txt";
    
    // Ścieżka wyjściowa
    string outputDir = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_rot_series/";
    string outputFile = outputDir + "series_1_2_3_4_5_6-51-52_7_9-1-13.txt";

    // Tablica na zsumowane dane (0-180)
    DataPoint totalData[181];

    cout << "Rozpoczynam sumowanie danych..." << endl;

    // Przetwarzanie pierwszego pliku
    if (processFile(file1, totalData)) {
        cout << "[OK] Przetworzono plik 1/2" << endl;
    }

    // Przetwarzanie drugiego pliku
    if (processFile(file2, totalData)) {
        cout << "[OK] Przetworzono plik 2/2" << endl;
    }

    // Tworzenie katalogu jeśli nie istnieje
    if (!fs::exists(outputDir)) {
        fs::create_directories(outputDir);
        cout << "[INFO] Utworzono katalog: " << outputDir << endl;
    }

    // Zapis do pliku wynikowego
    ofstream out(outputFile);
    if (out.is_open()) {
        out << "Angle;Time_Minutes;Counts\n";
        for (int i = 0; i <= 180; ++i) {
            // Zapisujemy tylko te kąty, które mają jakiekolwiek dane
            if (totalData[i].time > 0 || totalData[i].counts > 0) {
                out << i << ";" 
                    << fixed << setprecision(4) << totalData[i].time << ";" 
                    << totalData[i].counts << "\n";
            }
        }
        out.close();
        cout << "[SUKCES] Dane zostały zsumowane i zapisane w: " << outputFile << endl;
    } else {
        cerr << "[BŁĄD] Nie można utworzyć pliku wyjściowego!" << endl;
    }

    return 0;
}