// This function will calculate working time of detector
// Principle:
// Detector is rotating from 0 till 180 - changing 1 degree per 3 minutes
// Last digit in data set inform which direction detector was moving
// We calculate how many times detector changed direction as int count
// To calculate wokrin time W = 3 minutes x 180 x count
// We expect detectors to work for at least hours, so the final result
// should be in hours.

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int time(string file_dane);
string read_line(int position, int line_size, string line);

int main(){

    const int START = 11;   // from which data file to start
    const int END   = 11;   // where to finish

    for (int i = START; i <= END; ++i) {

        string data_file = "/home/kacper/MuonDensityMapping/bin/data_analysis/data/dane_" + to_string(i) + ".log";
        time(data_file);
    }
}

int time(string file_dane){
    fstream file;
    file.open(file_dane, ios::in);
    
    if (!file.is_open()) {
        cout << "[WARNING] Cannot open file: " << file_dane << endl;
        return -1;
    }

    string line;
    int a = -1;       // Używamy -1 jako stanu początkowego (niezainicjalizowanego)
    int count = 0;    // Zaczynamy od 0, pierwsze udane czytanie zmieni na 1

    while (getline(file, line)) {
        // Pomijamy puste i zepsute linie
        if (line.empty() || line == "\r") continue;
        if (line.rfind("Adafruit", 0) == 0) continue;
        if (line.rfind("Pod", 0) == 0) continue;
        
        string cell = read_line(13, line.size(), line);
        
        // Jeśli komórka jest pusta (linia była za krótka), ignorujemy tę linię
        if (cell.empty() || cell == "\r" || cell == " ") continue;

        try {
            // Próbujemy wyciągnąć liczbę z bezpiecznego już stringa
            int current_direction = stoi(cell);
            
            // Jeśli to pierwsza poprawna linia w dokumencie:
            if (a == -1) {
                a = current_direction;
                count = 1;
            } 
            // W przeciwnym razie sprawdzamy zmianę kierunku:
            else if (a != current_direction) {
                count += 1;
                a = current_direction;
            }
        } 
        // Zabezpieczenie na wypadek, gdyby w stringu były litery zamiast cyfr
        catch (const std::invalid_argument& e) {
            continue; 
        } 
        catch (const std::out_of_range& e) {
            continue; 
        }
    }
    file.close();

    // Jeśli plik był pusty lub miał same błędy
    if (count == 0) {
        cout << "No valid data in file." << endl;
        return 0;
    }

    // 3 minuty = 180 sekund * 180 stopni = 32400 sekund na jeden pełen ruch
    long long Working_time = 180LL * 180LL * count; 

    if (Working_time > 86400){
        int days = Working_time / 86400;
        int hours = (Working_time - (days * 86400)) / 3600;
        cout << "Detector was working aprox. " << days << " days " << hours << "h" << endl;
    }
    else if(Working_time < 32400){
        cout << "Detector experienced unexpected error (worked less than 1 full rotation)" << endl;
    }  
    else{
        cout << "Detector was working aprox. " << Working_time / 3600 << "h" << endl;
    }
    
    return 0;
}

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