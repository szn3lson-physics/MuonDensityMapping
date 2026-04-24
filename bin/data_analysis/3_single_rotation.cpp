#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstring> // dla memset

using namespace std;
namespace fs = std::filesystem;

// --- Deklaracje funkcji pomocniczych ---
string read_line(int position, int line_size, string line);
int read_bit(int position, int line_size, string line);
long long read_long_long(int position, int line_size, string line);

// Funkcja zapisująca wyniki jednego obrotu do pliku
void save_rotation_to_file(string filepath, double time_array[], long long count_array[]) {
    ofstream out(filepath, ios::trunc);
    out << "Angle;Time_Minutes;Counts_34f\n";
    for (int i = 0; i <= 180; ++i) {
        // Wypisujemy od 0 do 180 tak jak prosiłeś w formacie Angle;Time;Counts
        out << i << ";" << time_array[i] << ";" << count_array[i] << "\n";
    }
    out.close();
}

int main() {
    const int START_SERIES = 11;
    const int END_SERIES = 11;
    string INPUT_DIR = "/home/kacper/MuonDensityMapping/bin/data_analysis/data/"; 
    string OUTPUT_BASE_DIR = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/single_rotation/";

    for (int i = START_SERIES; i <= END_SERIES; ++i) {
        string inputFileName = INPUT_DIR + "dane_" + to_string(i) + ".log";
        string seriesDir = OUTPUT_BASE_DIR + to_string(i) + "_series/";

        ifstream inputFile(inputFileName);
        if (!inputFile.is_open()) {
            cout << "[OSTRZEZENIE] Nie mozna otworzyc pliku wejsciowego: " << inputFileName << endl;
            continue; 
        }

        if (!fs::exists(seriesDir)) {
            fs::create_directories(seriesDir);
        }

        string line;
        int current_direction = -1;
        int rotation_counter = 1;

        // Tablice do przechowywania danych dla pojedynczego obrotu
        double time_per_angle_min[181] = {0.0};
        long long count_34f_per_angle[181] = {0};

        // Zmienne do śledzenia czasu i kątów (połączona logika time_to_angle i calculate)
        long long t0_angle = 0; // Czas referencyjny dla przeskoku o 1 stopień
        long long t_enter = 0;  // Czas wejścia w dany kąt (do sumowania)
        long long t_prev = 0;   // Czas z poprzedniej linijki
        int current_angle = 0;
        int tracked_angle_for_time = 0; // Pomocnicza do zliczania czasu

        cout << "Przetwarzanie serii " << i << "... Rozpoczęto obrót " << rotation_counter << endl;

        while (getline(inputFile, line)) {
            if (line.rfind("Adafruit", 0) == 0) continue;
            if (line.rfind("Podłączenie karty SD", 0) == 0) continue;
            if (line.empty()) continue;

            long long t_current = read_long_long(1, line.size(), line);
            string coin = read_line(2, line.size(), line);
            int dir = read_bit(13, line.size(), line);
            
            if (dir == -1) continue; // Pusta lub błędna linia

            // --- INICJALIZACJA PIERWSZEJ LINII ---
            if (current_direction == -1) {
                current_direction = dir;
                t0_angle = t_current;
                t_enter = t_current;
                current_angle = (dir == 1) ? 1 : 180;
                tracked_angle_for_time = current_angle;
                t_prev = t_current;
            }

            // --- ZMIANA KIERUNKU = NOWY OBRÓT ---
            if (dir != current_direction) {
                // 1. Zakończ mierzenie czasu dla ostatniego kąta w starym obrocie
                if (tracked_angle_for_time >= 0 && tracked_angle_for_time <= 180 && t_prev >= t_enter) {
                    time_per_angle_min[tracked_angle_for_time] += (t_prev - t_enter) / 60000.0;
                }

                // 2. Zapisz stary obrót do pliku
                string outputFileName = seriesDir + to_string(rotation_counter) + "_rotation.txt";
                save_rotation_to_file(outputFileName, time_per_angle_min, count_34f_per_angle);

                // 3. Reset dla nowego obrotu
                memset(time_per_angle_min, 0, sizeof(time_per_angle_min));
                memset(count_34f_per_angle, 0, sizeof(count_34f_per_angle));
                rotation_counter++;
                current_direction = dir;
                
                t0_angle = t_current;
                t_enter = t_current;
                current_angle = (dir == 1) ? 1 : 180;
                tracked_angle_for_time = current_angle;
            } 
            // --- ANOMALIE CZASOWE (np. restart detektora) ---
            else if ((t_current - t_prev) > 600000 || t_current < t_prev) {
                if (tracked_angle_for_time >= 0 && tracked_angle_for_time <= 180 && t_prev >= t_enter) {
                    time_per_angle_min[tracked_angle_for_time] += (t_prev - t_enter) / 60000.0;
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
                    time_per_angle_min[tracked_angle_for_time] += (t_current - t_enter) / 60000.0;
                }
                tracked_angle_for_time = current_angle;
                t_enter = t_current;
            }

            // --- SPRAWDZANIE KOINCYDENCJI 34f ---
            if (current_angle >= 0 && current_angle <= 180) {
                if (coin == "000000" || 
                    coin == "001011" || 
                    coin == "111000" || 
                    coin == "010101" || 
                    coin == "100110") {
                    count_34f_per_angle[current_angle]++;
                }
            }

            t_prev = t_current;
        }

        // --- ZAPIS OSTATNIEGO OBROTU PO ZAKOŃCZENIU PLIKU ---
        if (current_direction != -1) {
            if (tracked_angle_for_time >= 0 && tracked_angle_for_time <= 180 && t_prev >= t_enter) {
                time_per_angle_min[tracked_angle_for_time] += (t_prev - t_enter) / 60000.0;
            }
            string outputFileName = seriesDir + to_string(rotation_counter) + "_rotation.txt";
            save_rotation_to_file(outputFileName, time_per_angle_min, count_34f_per_angle);
        }

        inputFile.close();
        cout << "Zakonczono serie " << i << ". Wygenerowano " << rotation_counter << " plikow obrotow.\n" << endl;
    }

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