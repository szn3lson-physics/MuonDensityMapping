#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstring> 

using namespace std;
namespace fs = std::filesystem;

// --- Deklaracje funkcji pomocniczych ---
string read_line(int position, int line_size, string line);
int read_bit(int position, int line_size, string line);
long long read_long_long(int position, int line_size, string line);

// Funkcja zapisująca zsumowane wyniki do pliku
void save_aggregated_to_file(string filepath, double time_array[], long long count_array[]) {
    ofstream out(filepath, ios::trunc);
    out << "Angle;Time_Minutes;Counts\n";
    for (int i = 0; i <= 180; ++i) {
        // Wypisujemy tylko kąty, dla których upłynął jakikolwiek czas lub były zliczenia
        if (time_array[i] > 0 || count_array[i] > 0) {
            out << i << ";" << time_array[i] << ";" << count_array[i] << "\n";
        }
    }
    out.close();
}

// Funkcja tworząca string z wektora (np. {7,8,9} -> "7_8_9")
string vec_to_string(const vector<int>& v) {
    string s = "";
    for (size_t i = 0; i < v.size(); ++i) {
        s += to_string(v[i]);
        if (i != v.size() - 1) s += "_";
    }
    return s;
}

int main() {
    // =========================================================================
    // KONFIGURACJA: WPISZ TUTAJ SERIE DO ZSUMOWANIA
    // =========================================================================
    vector<int> target_series = {1,2,3,4,5,7}; 
    // =========================================================================
    string INPUT_DIR = "/home/kacper/MuonDensityMapping/bin/data_analysis/data/"; 
    string OUTPUT_BASE_DIR = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_series/";

    // Tworzenie nazwy folderu wyjściowego na podstawie wektora
    string folder_name = "series_" + vec_to_string(target_series);
    string targetDir = OUTPUT_BASE_DIR + folder_name + "/";

    if (!fs::exists(targetDir)) {
        fs::create_directories(targetDir);
        cout << "Utworzono katalog docelowy: " << targetDir << endl;
    }

    // --- GLOBALNE TABLICE DLA WSZYSTKICH WYBRANYCH SERII ---
    double total_time_min[181] = {0.0};
    
    long long total_counts_3[181]   = {0};
    long long total_counts_4[181]   = {0};
    long long total_counts_34[181]  = {0};
    long long total_counts_3f[181]  = {0};
    long long total_counts_4f[181]  = {0};
    long long total_counts_34f[181] = {0};

    // Pętla po wskazanych seriach
    for (int series_id : target_series) {
        string inputFileName = INPUT_DIR + "dane_" + to_string(series_id) + ".log";

        ifstream inputFile(inputFileName);
        if (!inputFile.is_open()) {
            cout << "[OSTRZEZENIE] Nie mozna otworzyc pliku: " << inputFileName << ". Pomijam." << endl;
            continue; 
        }

        cout << "Przetwarzanie i agregowanie serii: " << series_id << "..." << endl;

        string line;
        int current_direction = -1;

        // Zmienne do śledzenia czasu i kątów dla bieżącego pliku
        long long t0_angle = 0; 
        long long t_enter = 0;  
        long long t_prev = 0;   
        int current_angle = 0;
        int tracked_angle_for_time = 0; 

        while (getline(inputFile, line)) {
            if (line.rfind("Adafruit", 0) == 0) continue;
            if (line.rfind("Podłączenie karty SD", 0) == 0) continue;
            if (line.rfind("Pod", 0) == 0) continue; // Dodatkowe zabezpieczenie
            if (line.empty()) continue;

            long long t_current = read_long_long(1, line.size(), line);
            string coin = read_line(2, line.size(), line);
            int dir = read_bit(13, line.size(), line);
            
            if (dir == -1 || coin.length() != 6) continue; // Pusta lub błędna linia

            // --- INICJALIZACJA PIERWSZEJ LINII W PLIKU ---
            if (current_direction == -1) {
                current_direction = dir;
                t0_angle = t_current;
                t_enter = t_current;
                current_angle = (dir == 1) ? 1 : 180;
                tracked_angle_for_time = current_angle;
                t_prev = t_current;
            }

            // --- ZMIANA KIERUNKU (Reset czasu, ale akumulacja trwa!) ---
            if (dir != current_direction) {
                if (tracked_angle_for_time >= 0 && tracked_angle_for_time <= 180 && t_prev >= t_enter) {
                    total_time_min[tracked_angle_for_time] += (t_prev - t_enter) / 60000.0;
                }
                current_direction = dir;
                t0_angle = t_current;
                t_enter = t_current;
                current_angle = (dir == 1) ? 1 : 180;
                tracked_angle_for_time = current_angle;
            } 
            // --- ANOMALIE CZASOWE (np. restart detektora) ---
            else if ((t_current - t_prev) > 600000 || t_current < t_prev) {
                if (tracked_angle_for_time >= 0 && tracked_angle_for_time <= 180 && t_prev >= t_enter) {
                    total_time_min[tracked_angle_for_time] += (t_prev - t_enter) / 60000.0;
                }
                t0_angle = t_current; 
                t_enter = t_current;
            }

            // --- LOGIKA AKTUALIZACJI KĄTA (co 3 minuty = 180000 ms) ---
            long long dt_angle = t_current - t0_angle;
            if (dt_angle >= 180000) {
                t0_angle = t_current - dt_angle + 180000;
                if (dir == 1) current_angle++; 
                else current_angle--;
            }

            // --- DODANIE CZASU JEŚLI KĄT SIĘ ZMIENIŁ ---
            if (current_angle != tracked_angle_for_time) {
                if (tracked_angle_for_time >= 0 && tracked_angle_for_time <= 180 && t_current >= t_enter) {
                    total_time_min[tracked_angle_for_time] += (t_current - t_enter) / 60000.0;
                }
                tracked_angle_for_time = current_angle;
                t_enter = t_current;
            }

            // --- SEGREGACJA KOINCYDENCJI ---
            if (current_angle >= 0 && current_angle <= 180) {
                bool is_4  = (coin == "000000");
                bool is_3  = (coin == "001011" || coin == "111000");
                bool is_3f = (is_3 || coin == "010101" || coin == "100110");

                if (is_3)  total_counts_3[current_angle]++;
                if (is_4)  total_counts_4[current_angle]++;
                if (is_3 || is_4) total_counts_34[current_angle]++;
                if (is_3f) total_counts_3f[current_angle]++;
                if (is_4)  total_counts_4f[current_angle]++; // 4f to po prostu 4
                if (is_3f || is_4) total_counts_34f[current_angle]++;
            }

            t_prev = t_current;
        }

        // Dodanie czasu z ostatniego segmentu w danym pliku
        if (current_direction != -1) {
            if (tracked_angle_for_time >= 0 && tracked_angle_for_time <= 180 && t_prev >= t_enter) {
                total_time_min[tracked_angle_for_time] += (t_prev - t_enter) / 60000.0;
            }
        }
        inputFile.close();
    }

    // =========================================================================
    // ZAPIS WYNIKÓW DO PLIKÓW W FOLDERZE DOCELOWYM
    // =========================================================================
    cout << "\nZapisywanie zsumowanych wynikow do folderu: " << targetDir << endl;

    save_aggregated_to_file(targetDir + "coin_3.txt",   total_time_min, total_counts_3);
    save_aggregated_to_file(targetDir + "coin_4.txt",   total_time_min, total_counts_4);
    save_aggregated_to_file(targetDir + "coin_34.txt",  total_time_min, total_counts_34);
    save_aggregated_to_file(targetDir + "coin_3f.txt",  total_time_min, total_counts_3f);
    save_aggregated_to_file(targetDir + "coin_4f.txt",  total_time_min, total_counts_4f);
    save_aggregated_to_file(targetDir + "coin_34f.txt", total_time_min, total_counts_34f);

    cout << "[SUKCES] Wygenerowano 6 plikow koincydencji!" << endl;

    return 0;
}

// ---------------------------------------------------------------------------
// Funkcje pomocnicze
// ---------------------------------------------------------------------------
string read_line(int position, int line_size, string line){
    string cell = "";
    int count = 0;
    for (int i = 0; i < line_size; i++) {
        char c = line[i];
        if (c == ';') {count++; continue;} 
        if (count == position) cell += c;  
        else if (count > position) break;  
    }
    return cell;
}

int read_bit(int position, int line_size, string line) {
    string cell = read_line(position, line_size, line);
    if(cell.empty()) return -1; 
    return cell[0] - '0';
}

long long read_long_long(int position, int line_size, string line) {
    string cell = read_line(position, line_size, line);
    if(cell.empty()) return 0;
    return stoll(cell);
}