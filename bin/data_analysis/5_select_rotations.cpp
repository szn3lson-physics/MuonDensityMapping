#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <sstream>
#include <algorithm>

using namespace std;
namespace fs = std::filesystem;

// =========================================================================
// FUNKCJE POMOCNICZE DO PARSOWANIA ZAKRESÓW (np. "1-5, 7, 10-12")
// =========================================================================

// Funkcja rozwija string z zakresami do wektora pojedynczych intów
vector<int> parse_ranges(const string& input) {
    vector<int> result;
    stringstream ss(input);
    string token;

    while (getline(ss, token, ',')) {
        // Usuwanie spacji z tokena
        token.erase(remove_if(token.begin(), token.end(), ::isspace), token.end());
        if (token.empty()) continue;

        size_t dash_pos = token.find('-');
        if (dash_pos != string::npos) {
            // Znaleziono zakres (np. "14-23")
            int start = stoi(token.substr(0, dash_pos));
            int end = stoi(token.substr(dash_pos + 1));
            for (int i = start; i <= end; ++i) {
                result.push_back(i);
            }
        } else {
            // Pojedyncza liczba (np. "27")
            result.push_back(stoi(token));
        }
    }
    return result;
}

// Funkcja zamienia przecinki na podkreślniki, tworząc bezpieczną nazwę pliku
string get_safe_filename(string input) {
    input.erase(remove_if(input.begin(), input.end(), ::isspace), input.end()); // usuń spacje
    replace(input.begin(), input.end(), ',', '_'); // zamień przecinki na _
    return input;
}

int main() {
    // =========================================================================
    // KONFIGURACJA: WYBIERZ SERIĘ I ZAKRES ROTACJI DO ZSUMOWANIA
    // =========================================================================
    int target_series = 9; // Numer serii, z której pobierasz rotacje
    
    // Wpisz tutaj swój zakres tak jak lubisz: przecinki oddzielają, myślniki to zakresy
    string target_rotations_str = "14-43"; 
    // =========================================================================

    // Automatyczne rozwinięcie stringa do wektora {14, 15, ..., 23, 27, 29, 30, 31...}
    vector<int> target_rotations = parse_ranges(target_rotations_str);

    string INPUT_BASE_DIR = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/single_rotation/";
    string OUTPUT_BASE_DIR = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_rotations_900/";

    // Określenie folderów dla danej serii
    string series_input_dir = INPUT_BASE_DIR + to_string(target_series) + "_series/";
    string series_output_dir = OUTPUT_BASE_DIR + to_string(target_series) + "_series/";

    if (!fs::exists(series_output_dir)) {
        fs::create_directories(series_output_dir);
        cout << "Utworzono katalog docelowy: " << series_output_dir << endl;
    }

    // Globalne tablice do sumowania czasu i zliczeń dla wybranej grupy rotacji
    double total_time[181] = {0.0};
    long long total_counts[181] = {0};

    cout << "Rozpoczynam agregacje dla serii " << target_series << "..." << endl;
    cout << "Wybrane rotacje: " << target_rotations_str << endl;
    cout << "Lacznie rotacji do sprawdzenia: " << target_rotations.size() << "\n" << endl;

    // Pętla po wszystkich wskazanych przez Ciebie rotacjach
    for (int rot : target_rotations) {
        string input_filepath = series_input_dir + to_string(rot) + "_rotation.txt";
        ifstream infile(input_filepath);

        if (!infile.is_open()) {
            cout << "[OSTRZEZENIE] Nie znaleziono pliku rotacji: " << rot << "_rotation.txt. Pomijam." << endl;
            continue;
        }

        string line;
        // Pominięcie pierwszego wiersza z nagłówkami (Angle;Time_Minutes;Counts_34f)
        getline(infile, line);

        while (getline(infile, line)) {
            if (line.empty()) continue;

            stringstream ss(line);
            string token;
            int angle = -1;
            double time_min = 0.0;
            long long counts = 0;

            // Wczytywanie w formacie: Angle ; Time_Minutes ; Counts_34f
            if (getline(ss, token, ';')) angle = stoi(token);
            if (getline(ss, token, ';')) time_min = stod(token);
            if (getline(ss, token, ';')) counts = stoll(token);

            // Sumowanie wartości dla danego kąta
            if (angle >= 0 && angle <= 180) {
                total_time[angle] += time_min;
                total_counts[angle] += counts;
            }
        }
        infile.close();
        cout << "-> Zsumowano dane z rotacji: " << rot << endl;
    }

    // Wygenerowanie eleganckiej nazwy pliku na podstawie wpisanego stringa
    string output_filename = series_output_dir + get_safe_filename(target_rotations_str) + "_rotations.txt";
    ofstream outfile(output_filename, ios::trunc);

    // Zapis zsumowanych wyników (zawsze wypisujemy od 0 do 180)
    outfile << "Angle;Time_Minutes;Counts_34f\n";
    for (int i = 0; i <= 180; ++i) {
        outfile << i << ";" << total_time[i] << ";" << total_counts[i] << "\n";
    }
    outfile.close();

    cout << "\n[SUKCES] Zapisano zagregowane dane do pliku:\n" << output_filename << endl;

    return 0;
}