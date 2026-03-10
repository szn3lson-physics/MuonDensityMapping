#include <iostream>
#include <fstream>
#include <string>
#include <conio.h>
#include <vector>
#include <filesystem>

using namespace std;

void add_to_file (string output_file_1, string output_file_2);
void folder_managment(string file);

int main(){
    const int START = 1;   // from which data file to start
    const int END   = 9;   // where to finish

    folder_managment("D:/MuonDensityMapping/output/all/" + to_string(START) + "_zussamen");
    filesystem::copy_file("D:/MuonDensityMapping/output/dane_1/processed/coin_1_3.txt", "D:/MuonDensityMapping/output/all/" + to_string(START) + "_zussamen/coin__3.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file("D:/MuonDensityMapping/output/dane_1/processed/coin_1_3f.txt", "D:/MuonDensityMapping/output/all/" + to_string(START) + "_zussamen/coin__3f.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file("D:/MuonDensityMapping/output/dane_1/processed/coin_1_4.txt", "D:/MuonDensityMapping/output/all/" + to_string(START) + "_zussamen/coin__4.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file("D:/MuonDensityMapping/output/dane_1/processed/coin_1_34.txt", "D:/MuonDensityMapping/output/all/" + to_string(START) + "_zussamen/coin__34.txt", filesystem::copy_options::overwrite_existing);
    filesystem::copy_file("D:/MuonDensityMapping/output/dane_1/processed/coin_1_34f.txt", "D:/MuonDensityMapping/output/all/" + to_string(START) + "_zussamen/coin__34f.txt", filesystem::copy_options::overwrite_existing);

    for (int i = START + 1; i <= END; ++i){
        folder_managment("D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen");
        filesystem::copy_file("D:/MuonDensityMapping/output/all/" + to_string(i-1) + "_zussamen/coin__3.txt", "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__3.txt", filesystem::copy_options::overwrite_existing);
        add_to_file ("D:/MuonDensityMapping/output/dane_" + to_string(i) + "/processed/coin_" + to_string(i) + "_3.txt",
                     "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__3.txt");
        
        filesystem::copy_file("D:/MuonDensityMapping/output/all/" + to_string(i-1) + "_zussamen/coin__3f.txt", "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__3f.txt", filesystem::copy_options::overwrite_existing);
        add_to_file ("D:/MuonDensityMapping/output/dane_" + to_string(i) + "/processed/coin_" + to_string(i) + "_3f.txt",
                     "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__3f.txt");
        
        filesystem::copy_file("D:/MuonDensityMapping/output/all/" + to_string(i-1) + "_zussamen/coin__4.txt", "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__4.txt", filesystem::copy_options::overwrite_existing);
        add_to_file ("D:/MuonDensityMapping/output/dane_" + to_string(i) + "/processed/coin_" + to_string(i) + "_4.txt",
                     "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__4.txt");
                    
        filesystem::copy_file("D:/MuonDensityMapping/output/all/" + to_string(i-1) + "_zussamen/coin__34.txt", "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__34.txt", filesystem::copy_options::overwrite_existing);
        add_to_file ("D:/MuonDensityMapping/output/dane_" + to_string(i) + "/processed/coin_" + to_string(i) + "_34.txt",
                     "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__34.txt");
        
        filesystem::copy_file("D:/MuonDensityMapping/output/all/" + to_string(i-1) + "_zussamen/coin__34f.txt", "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__34f.txt", filesystem::copy_options::overwrite_existing);
        add_to_file ("D:/MuonDensityMapping/output/dane_" + to_string(i) + "/processed/coin_" + to_string(i) + "_34f.txt",
                     "D:/MuonDensityMapping/output/all/" + to_string(i) + "_zussamen/coin__34f.txt");
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