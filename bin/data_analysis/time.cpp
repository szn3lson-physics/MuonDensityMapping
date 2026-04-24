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

int time (string file_dane);
string read_line(int position, int line_size, string line);

int main(){

    const int START = 11;   // from which data file to start
    const int END   = 12;   // where to finish

    for (int i = START; i <= END; ++i) {

        string data_file = "/home/kacper/MuonDensityMapping/data/dane_" + to_string(i) + ".log";
        time(data_file);
    }
}

int time (string file_dane){
    //we need to read the file obviosuly
    fstream file;
    file.open(file_dane, ios::in);
    //cout << "Access to the file " << file_dane << " has been granted!" << endl;

    string line;

    getline(file, line);
    if (line.rfind("Adafruit", 0) == 0) getline(file, line);
    if (line.rfind("Pod", 0) == 0) getline(file, line);

    int a = stoi(read_line(13, line.size(), line));
    //cout << "initial a wynosi" << a << endl;
    int count = 1;

    while (getline(file, line)) {
            if (line.empty()) continue;
            if (line.rfind("Adafruit", 0) == 0) continue;
            if (line.rfind("Pod", 0) == 0) continue;
            else{
                if (a != stoi(read_line(13, line.size(), line)))
                {
                    count += 1;
                    a = stoi(read_line(13, line.size(), line));
                    //cout << "nowe a wynosi" << a << endl;
                }
            }
        }
    file.close();

    int Working_time = 180 * 180 * count;

    if (Working_time > 86400){
        int days = Working_time/86400;
        int hours = (Working_time - days*86400)/3600;
        cout << "Detector was working aprox. " << days << "days" << hours << "h"<<endl;
        return 0;
    }
    else if(Working_time < 32400){
        cout << "Detector experienced unexpected error"<<endl;
        return 0;
    }  
    else{
        cout << "Detector was working aprox. " << Working_time/3600 << "h"<<endl;
        return 0;
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
