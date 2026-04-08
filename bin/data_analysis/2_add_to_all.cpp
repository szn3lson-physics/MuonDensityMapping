#include <iostream>
#include <fstream>
#include <string>
#include <conio.h>
#include <vector>
#include <filesystem>

using namespace std;

void add_to_file (string output_file_1, string output_file_2);
void folder_managment(string file);
string read_line(int position, int line_size, string line);
void sum_time_and_counts(int start_idx, int end_idx, string output_file);

    string single_directory = "D:/MuonDensityMapping/output/single/dane_";
    string all_directory = "D:/MuonDensityMapping/output/all/zussamen_";

int main(){
    const int START = 1;   // from which data file to start
    const int END   = 10;   // where to finish

    //folder_managment("D:/MuonDensityMapping/output/all");
    folder_managment(all_directory + to_string(START));
    
    filesystem::copy_file(single_directory + to_string(START) + "/processed/coin_" + to_string(START) + "_3.txt", 
        all_directory + to_string(START) + "/coin_3.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file(single_directory + to_string(START) + "/processed/coin_" + to_string(START) + "_3f.txt", 
        all_directory + to_string(START) + "/coin_3f.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file(single_directory + to_string(START) + "/processed/coin_" + to_string(START) + "_4.txt", 
        all_directory + to_string(START) + "/coin_4.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file(single_directory + to_string(START) + "/processed/coin_" + to_string(START) + "_34.txt", 
        all_directory + to_string(START) + "/coin_34.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file(single_directory + to_string(START) + "/processed/coin_" + to_string(START) + "_34f.txt", 
        all_directory + to_string(START) + "/coin_34f.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file(single_directory + to_string(START) + "/processed/coin_" + to_string(START) + "_true.txt", 
        all_directory + to_string(START) + "/coin_true.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file(single_directory + to_string(START) + "/processed/coin_" + to_string(START) + "_false.txt", 
        all_directory + to_string(START) + "/coin_false.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file(single_directory + to_string(START) + "/processed/time_per_angle_" + to_string(START) + ".txt", 
        all_directory + to_string(START) + "/coin_time.txt", filesystem::copy_options::overwrite_existing);

    for (int i = START + 1; i <= END; ++i){
        folder_managment(all_directory + to_string(i));
        
        filesystem::copy_file(all_directory + to_string(i-1) + "/coin_3.txt", all_directory + to_string(i) + "/coin_3.txt", filesystem::copy_options::overwrite_existing);
        add_to_file (single_directory + to_string(i) + "/processed/coin_" + to_string(i) + "_3.txt",
                     all_directory + to_string(i) + "/coin_3.txt");
        
        filesystem::copy_file(all_directory + to_string(i-1) + "/coin_3f.txt", all_directory + to_string(i) + "/coin_3f.txt", filesystem::copy_options::overwrite_existing);
        add_to_file (single_directory + to_string(i) + "/processed/coin_" + to_string(i) + "_3f.txt",
                     all_directory + to_string(i) + "/coin_3f.txt");
        
        filesystem::copy_file(all_directory + to_string(i-1) + "/coin_4.txt", all_directory + to_string(i) + "/coin_4.txt", filesystem::copy_options::overwrite_existing);
        add_to_file (single_directory + to_string(i) + "/processed/coin_" + to_string(i) + "_4.txt",
                     all_directory + to_string(i) + "/coin_4.txt");
                    
        filesystem::copy_file(all_directory + to_string(i-1) + "/coin_34.txt", all_directory + to_string(i) + "/coin_34.txt", filesystem::copy_options::overwrite_existing);
        add_to_file (single_directory + to_string(i) + "/processed/coin_" + to_string(i) + "_34.txt",
                     all_directory + to_string(i) + "/coin_34.txt");
        
        filesystem::copy_file(all_directory + to_string(i-1) + "/coin_34f.txt", all_directory + to_string(i) + "/coin_34f.txt", filesystem::copy_options::overwrite_existing);
        add_to_file (single_directory + to_string(i) + "/processed/coin_" + to_string(i) + "_34f.txt",
                     all_directory + to_string(i) + "/coin_34f.txt");
                     
        filesystem::copy_file(all_directory + to_string(i-1) + "/coin_true.txt", all_directory + to_string(i) + "/coin_true.txt", filesystem::copy_options::overwrite_existing);
        add_to_file (single_directory + to_string(i) + "/processed/coin_" + to_string(i) + "_true.txt",
                     all_directory + to_string(i) + "/coin_true.txt");

        filesystem::copy_file(all_directory + to_string(i-1) + "/coin_false.txt", all_directory + to_string(i) + "/coin_false.txt", filesystem::copy_options::overwrite_existing);
        add_to_file (single_directory + to_string(i) + "/processed/coin_" + to_string(i) + "_false.txt",
                     all_directory + to_string(i) + "/coin_false.txt");

        sum_time_and_counts(START, i, all_directory + to_string(i) + "/coin_time.txt");
    }
}



void add_to_file (string output_file_1, string output_file_2){
    ifstream file_1(output_file_1);
    ofstream file_2(output_file_2, std::ios::app);
    file_2 << "\n";
    file_2 << file_1.rdbuf();
}

void folder_managment(string file){
    if (filesystem::exists(file)) {
        if (filesystem::is_directory(file)) {
            cout << "Folder '" << file << "' already exist." << endl;
        } else {
            cout << "Error - name conflict: File '" << file << "'exist, cannot make folder" << endl;
        }
        } else {
            // Jeśli nie istnieje, tworzymy folder
            if (filesystem::create_directory(file)) {
                cout << "Succesufully created file: '" << file << "'" << endl;
            } else {
                cout << "There was an error creating file" << endl;
            }
        }
}

string read_line(int position, int line_size, string line){
    string cell = "";
    int count = 0;
    for (int i = 0; i < line_size; i++) {
        char c = line[i];
        if (c == ';') {count++; continue;} // Count separators and skip them
        if (count == position) cell += c;  // If we reached the correct field, collect its characters
        else if (count > position) break;  // If we passed beyond the target field, stop reading
    }
    return cell;
}

// ---------------------------------------------------------------------------
// sum_time_and_counts()
// ---------------------------------------------------------------------------
// Otwiera pliki z poszczególnych serii, matematycznie sumuje czas i zliczenia
// dla każdego kąta i zapisuje wynik do jednego pliku zbiorczego.
// Wykorzystuje Twoją natywną funkcję read_line do parsowania.
// ---------------------------------------------------------------------------
void sum_time_and_counts(int start_idx, int end_idx, string output_file) {
    double total_time[181] = {0.0};
    long long total_counts[181] = {0};

    for (int i = start_idx; i <= end_idx; ++i) {
        string input_file = single_directory + to_string(i) + "/processed/time_per_angle_" + to_string(i) + ".txt";
        ifstream in(input_file);
        
        if (!in.is_open()) {
            cout << "[OSTRZEZENIE] Nie mozna otworzyc pliku do sumowania: " << input_file << endl;
            continue;
        }

        string line;
        getline(in, line); // Pomijamy nagłówek "Angle;Time_Minutes;Counts_34f"

        while (getline(in, line)) {
            if (line.empty()) continue;
            
            string ang_str = read_line(0, line.size(), line);
            string time_str = read_line(1, line.size(), line);
            string count_str = read_line(2, line.size(), line);

            if (!ang_str.empty() && !time_str.empty() && !count_str.empty()) {
                int ang = stoi(ang_str);
                double time_min = stod(time_str);
                long long counts = stoll(count_str);

                if (ang >= 0 && ang <= 180) {
                    total_time[ang] += time_min;
                    total_counts[ang] += counts;
                }
            }
        }
        in.close();
    }

    // Zapis zsumowanych wartości do pliku docelowego
    ofstream out(output_file, ios::trunc);
    out << "Angle;Time_Minutes;Counts_34f\n";
    for (int i = 0; i <= 180; ++i) {
        if (total_time[i] > 0 || total_counts[i] > 0) {
            out << i << ";" << total_time[i] << ";" << total_counts[i] << "\n";
        }
    }
    out.close();
    cout << "Zsumowano matematycznie czas i zliczenia. Plik docelowy: " << output_file << endl;
}