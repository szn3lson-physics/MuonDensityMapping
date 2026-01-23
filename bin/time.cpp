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

// In the x place put the file number
string file_name = "../Data/time_data/dane_x.log";

string read_line(int position, int line_size, string line);

int main(){

    //we need to read the file obviosuly
    fstream file;
    file.open(file_name, ios::in);
    cout << "Access to the file " << file_name << " has been granted!" << endl;

    string line;

    getline(file, line);
    int a = stoi(read_line(13, line.size(), line));
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
                }
            }
        }
    file.close();

    int Working_time = 180 * 180 * count;

    if (Working_time > 86400){
        int days = Working_time/86400;
        int hours = (Working_time - days*86400)/3600;
        cout << "Detector was working aprox. " << days << "days" << hours << "h";
        return 0;
    }
    else if(Working_time < 32400){
        cout << "Detector experienced unexpected error";
        return 0;
    }  
    else{
        cout << "Detector was working aprox. " << Working_time/3600 << "h";
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
