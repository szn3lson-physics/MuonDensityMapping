#include <iostream>
#include <fstream>
#include <string>
#include <conio.h>
#include <vector>

using namespace std;

string read_line(int position, int size, string line);
int read_int(int position, int line_size, string line);
string read_binary_string(int position, int line_size, string line);
int read_bit(int position, int line_size, string line);
void time_to_angle(const string& file_angles);
void output_data(string file_dane, string file_angles);
void coindidence_4(string file_angles, string file_coindicence);
void coindidence_3_4(string file_angles, string file_coindicence);
void coindidence_3_4f(string file_angles, string file_coindicence);
void coindidence_3f(string file_angles, string file_coindicence);
void coindidence_3(string file_angles, string file_coindicence);
void process(string inputFileName, string outputFileName, int multiplier);

void coindidence_134_234(string file_angles, string file_coindicence);
void coindidence_123_124(string file_angles, string file_coindicence);


int main (){
    //excavate only nessessery data
    output_data("../Data/dane_1.log", "../Output/dane_1/angles_1.txt");
    output_data("../Data/dane_2.log", "../Output/dane_2/angles_2.txt");
    output_data("../Data/dane_3.log", "../Output/dane_3/angles_3.txt");
    output_data("../Data/dane_4.log", "../Output/dane_4/angles_4.txt");
    
    //take excaveted data and assign to each line angles for which detection was done
    time_to_angle("../Output/dane_1/angles_1.txt");
    time_to_angle("../Output/dane_2/angles_2.txt");
    time_to_angle("../Output/dane_3/angles_3.txt");
    time_to_angle("../Output/dane_4/angles_4.txt");

    //take excaveted data with angles and search for needed coincidences
    //raw means that in data there could be time where detector stopped working
    //this data for now have to be deleted manually
    coindidence_3("../Output/dane_1/angles_1.txt", "../Output/dane_1/raw/coin_1_3.txt");
    coindidence_3("../Output/dane_2/angles_2.txt", "../Output/dane_2/raw/coin_2_3.txt");
    coindidence_3("../Output/dane_3/angles_3.txt", "../Output/dane_3/raw/coin_3_3.txt");
    coindidence_3("../Output/dane_4/angles_4.txt", "../Output/dane_4/raw/coin_4_3.txt");

    coindidence_3f("../Output/dane_1/angles_1.txt", "../Output/dane_1/raw/coin_1_3f.txt");
    coindidence_3f("../Output/dane_2/angles_2.txt", "../Output/dane_2/raw/coin_2_3f.txt");
    coindidence_3f("../Output/dane_3/angles_3.txt", "../Output/dane_3/raw/coin_3_3f.txt");
    coindidence_3f("../Output/dane_4/angles_4.txt", "../Output/dane_4/raw/coin_4_3f.txt");

    coindidence_4("../Output/dane_1/angles_1.txt", "../Output/dane_1/raw/coin_1_4.txt");
    coindidence_4("../Output/dane_2/angles_2.txt", "../Output/dane_2/raw/coin_2_4.txt");
    coindidence_4("../Output/dane_3/angles_3.txt", "../Output/dane_3/raw/coin_3_4.txt");
    coindidence_4("../Output/dane_4/angles_4.txt", "../Output/dane_4/raw/coin_4_4.txt");

    coindidence_3_4("../Output/dane_1/angles_1.txt", "../Output/dane_1/raw/coin_1_34.txt");
    coindidence_3_4("../Output/dane_2/angles_2.txt", "../Output/dane_2/raw/coin_2_34.txt");
    coindidence_3_4("../Output/dane_3/angles_3.txt", "../Output/dane_3/raw/coin_3_34.txt");
    coindidence_3_4("../Output/dane_4/angles_4.txt", "../Output/dane_4/raw/coin_4_34.txt");

    coindidence_3_4f("../Output/dane_1/angles_1.txt", "../Output/dane_1/raw/coin_1_34f.txt");
    coindidence_3_4f("../Output/dane_2/angles_2.txt", "../Output/dane_2/raw/coin_2_34f.txt");
    coindidence_3_4f("../Output/dane_3/angles_3.txt", "../Output/dane_3/raw/coin_3_34f.txt");
    coindidence_3_4f("../Output/dane_4/angles_4.txt", "../Output/dane_4/raw/coin_4_34f.txt");
    
    coindidence_134_234("../Output/all/angles_all.txt", "../Output/all/coin_134_234.txt");
    coindidence_123_124("../Output/all/angles_all.txt", "../Output/all/coin_123_124.txt");

    //process data for histograms (each bar have to have 3 degrees)
    process("../Output/coin_x_x.txt", "../Output/processed_x_x.txt", 3);
    process("../Output/coin_x_x.txt", "../Output/processed_x_x.txt", 3);
    
    return 0;
}

// ---------------------------------------------------------------------------
// Extracts a specific field from a semicolon-separated line.
// position  – index of the field to extract (0-based)
// line_size – length of the input line
// line      – the full semicolon-delimited string
//
// The function scans the line character by character, counts semicolons,
// and collects characters belonging to the requested field.
// ---------------------------------------------------------------------------
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
// Reads a field from a semicolon-separated line and converts it to an integer.
// Internally uses read_line() and applies stoi() to the extracted substring.
// ---------------------------------------------------------------------------
int read_int(int position, int line_size, string line)
{
    return stoi(read_line(position, line_size, line)); //stoi - string to intiger
}

// ---------------------------------------------------------------------------
// Reads a field from a semicolon-separated line and returns it as a string.
// This is essentially a direct wrapper around read_line().
// ---------------------------------------------------------------------------
string read_binary_string(int position, int line_size, string line)
{
    return read_line(position, line_size, line);
}

// ---------------------------------------------------------------------------
// Reads a single-bit value (0 or 1) from a semicolon-separated line.
// Extracts the field at the given position using read_line(),
// takes its first character, and converts it from char to integer.
// Example: '0' → 0, '1' → 1.
// ---------------------------------------------------------------------------
int read_bit(int position, int line_size, string line)
{
    return read_line(position, line_size, line)[0] - '0';
}

// ---------------------------------------------------------------------------
// time_to_angle()
// ---------------------------------------------------------------------------
// Reads a semicolon-separated data file where each line contains:
//   - a timestamp (position 0)
//   - coincidence for ex. '011111'
//   - a rotation direction bit (position 2): 0 or 1
//
// The function calculates an angle value for each line based on:
//
//   • rotation direction:
//        v3 = 0 → angle increases from   0 to 180 degrees
//        v3 = 1 → angle decreases from 180 to   0 degrees
//
//   • time difference between measurements:
//        every 180,000 ms (3 minutes) → angle changes by 1 degree
//
// If the direction bit changes between lines, the angle resets:
//        v3 = 0 → 0 degrees
//        v3 = 1 → 180 degrees
//
// The angle is appended to each line as a new semicolon-separated field,
// and the modified data is written back to the same file.
// ---------------------------------------------------------------------------
void time_to_angle(const string& file_angles)
{
    cout << "The time_to_angle() function has been activated." << endl;

    ifstream file(file_angles);
    if (!file.good()) {
        cout << "Access to the file " << file_angles << " has been denied!" << endl;
        return;
    }

    cout << "Access to the file " << file_angles << " has been granted!" << endl;

    vector<string> lines;
    string line;

    // --- Get first line ----
    getline(file, line);
    lines.push_back(line);

    // --- Get first value from the first line - initial time & rotation
    int t0 = read_int(0, line.size(), line);        //initial time
    int v3_prev = read_bit(2, line.size(), line);   //rotation 0 or 1 determines direction

    int angle_inc = 0;     // for v3 = 0 - increasing to 180
    int angle_dec = 180;   // for v3 = 1 - decreasing to 0

    // --- Angle of the first line ---
    if (v3_prev == 0)
        lines[0] += ";" + to_string(angle_inc);
    else
        lines[0] += ";" + to_string(angle_dec);

    // ---- Get next lines ----
    while (getline(file, line)) {
        int t1 = read_int(0, line.size(), line);
        int dt = t1 - t0;
        int v3 = read_bit(2, line.size(), line);

        // --- If the direction changes → reset the counter ---
        if (v3 != v3_prev) {
            angle_inc = 0;
            angle_dec = 180;
            t0 = t1;  // reference time reset
        }

        // --- Update every 180,000 ms = 1 degree change---
        if (dt >= 180000) {
            t0 = t1;

            if (v3 == 0)
                angle_inc++;
            else
                angle_dec--;
        }

        // --- Adding an angle to the end of a line ---
        if (v3 == 0)
            line += ";" + to_string(angle_inc);
        else
            line += ";" + to_string(angle_dec);

        lines.push_back(line);

        v3_prev = v3; // remembering the direction for the next iteration
    }

    file.close();
    cout << "The file " << file_angles <<" has been closed" << endl;

    // write vector to file
    ofstream out(file_angles, ios::trunc); //delete the old file -> print the values from vector
    for (auto& l : lines)
        out << l << "\n";
    out.close();

    cout << "Processing of the file " << file_angles <<" has been completed and safed." << endl;
}

// ---------------------------------------------------------------------------
// output_data()
// ---------------------------------------------------------------------------
// Reads a source text file (file_data), processes its lines, and writes
// selected semicolon-separated fields into a new output file (file_angles).
//
//   • Lines beginning with "Adafruit" or "Podłączenie karty SD" are skipped.
//   • For all other lines, the function extracts:
//         - field at position 1
//         - field at position 2
//         - field at position 13
//     using read_line(), and writes them to the output file in the format:
//         field1;field2;field13
// ---------------------------------------------------------------------------
void output_data(string file_data, string file_angles){
    fstream data;
    fstream file;
    data.open(file_data, ios::in);

    if(data.good()){
        cout << "Access to the file " << file_data << " has been granted!" << endl;
        string linia;
        
        file.open(file_angles, ios::out);
        cout << "Output file" << file_angles << " has been created." << endl;
        
        while (getline(data, linia)) {
            if (linia.rfind("Adafruit", 0) == 0) continue;
            if (linia.rfind("Podłączenie karty SD", 0) == 0) continue;
            else{
                file << read_line(1, linia.size(), linia) << ";";
                file << read_line(2, linia.size(), linia) << ";";                
                file << read_line(13, linia.size(), linia) << endl;
            }
        }
        file.close();
        cout << "Data saved successfully in "<< file_angles << endl;
    }
    else {
        cout << "Access to the file has been denied!" << endl;
    }
    data.close();
}

// ---------------------------------------------------------------------------
// coincidence()
// ---------------------------------------------------------------------------
// Reads an input file (file_angles) and filters the lines based on field 1.
// Only lines where field 1 equals "000000" are written to the output file
// (file_coindicence). For those lines, the function outputs fields
// 0, 1, 2, and 3 in semicolon-separated format.
//
// Example output line:
//     value0 
//     value1
//     value2
//     value3
// ---------------------------------------------------------------------------
void coindidence_4(string file_angles, string file_coindicence){
    fstream data;
    data.open(file_angles, ios::in);

    if(data.good()){
        cout << "Acess to file has been granted!" << endl;
        
        string line;
        fstream file(file_coindicence, ios::out);

        while (!data.eof()) {
            getline(data, line);

            // If field 1 equals "000000", write selected fields
            if (read_line(1, line.size(), line) == "000000"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
        }
    } 
    else {
        cout << "Access to file was denided!" << endl;
    }

    data.close();
}

// ---------------------------------------------------------------------------
// coincidence_3_4()
// ---------------------------------------------------------------------------
// Reads an input file (file_angles) and filters the lines based on field 1.
// Only lines where field 1 are "000000" or "001011" or "111000" reffered as 
// true three coincidences (only combinations where muon pass through detector 
// 1-2, or 3-4, all hits all of them), then are written to the output file 
// (file_coindicence). For those lines, the function outputs fields 0, 1, 2, 
// and 3 in semicolon-separated format.
//
// Example output line:
//     value0 
//     value1
//     value2
//     value3
// ---------------------------------------------------------------------------
void coindidence_3_4(string file_angles, string file_coindicence){
    fstream data;
    data.open(file_angles, ios::in);

    if(data.good()){
        cout << "Acess to file has been granted!" << endl;
        
        string line;
        fstream file(file_coindicence, ios::out);

        while (!data.eof()) {
            getline(data, line);

            // If field 1 equals "000000", write selected fields
            if (read_line(1, line.size(), line) == "000000"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "001011"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "111000"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
        }
    } 
    else {
        cout << "Access to file was denided!" << endl;
    }

    data.close();
}

// ---------------------------------------------------------------------------
// coincidence_3_4f()
// ---------------------------------------------------------------------------
// Reads an input file (file_angles) and filters the lines based on field 1.
// Only lines where field 1 are all possible three and four coincidences, even the
// ones that are not physical - i.e. detector 1 and 3. They are taken in coside-
// ration due to effectiveness of detectors which is below 100%. Then are written 
// to the output file (file_coindicence). For those lines, the function outputs 
// fields 0, 1, 2, and 3 in semicolon-separated format.
//
// Example output line:
//     value0 
//     value1
//     value2
//     value3
// ---------------------------------------------------------------------------
void coindidence_3_4f(string file_angles, string file_coindicence){
    fstream data;
    data.open(file_angles, ios::in);

    if(data.good()){
        cout << "Acess to file has been granted!" << endl;
        
        string line;
        fstream file(file_coindicence, ios::out);

        while (!data.eof()) {
            getline(data, line);

            // If field 1 equals "000000", write selected fields
            if (read_line(1, line.size(), line) == "000000"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "001011"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "111000"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "010101"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "100110"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
        }
    } 
    else {
        cout << "Access to file was denided!" << endl;
    }

    data.close();
}

// ---------------------------------------------------------------------------
// coincidence_3f()
// ---------------------------------------------------------------------------
// Reads an input file (file_angles) and filters the lines based on field 1.
// Only lines where field 1 are all possible three three coincidences, even the
// ones that are not physical - i.e. detector 1, 3 and 4. They are taken in 
// cosideration due to effectiveness of detectors which is below 100%. Then are 
// written to the output file (file_coindicence). For those lines, the function 
// outputs fields 0, 1, 2, and 3 in semicolon-separated format.
//
// Example output line:
//     value0 
//     value1
//     value2
//     value3
// ---------------------------------------------------------------------------
void coindidence_3f(string file_angles, string file_coindicence){
    fstream data;
    data.open(file_angles, ios::in);

    if(data.good()){
        cout << "Acess to file has been granted!" << endl;
        
        string line;
        fstream file(file_coindicence, ios::out);

        while (!data.eof()) {
            getline(data, line);
            if (read_line(1, line.size(), line) == "001011"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "111000"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "010101"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "100110"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
        }
    } 
    else {
        cout << "Access to file was denided!" << endl;
    }

    data.close();
}

// ---------------------------------------------------------------------------
// coincidence_3()
// ---------------------------------------------------------------------------
// Reads an input file (file_angles) and filters the lines based on field 1.
// Only lines where field 1 are all possible three three coincidences - ONLY
// THE PHYSICAL ONCES so detectors 1, 2 and 3 or 2, 3 and 4- without f - false 
// coincidences. Then are written to  the output file (file_coindicence). For 
// those lines, the function outputs fields 0, 1, 2, and 3 in semicolon-separa-
// ted format.
//
// Example output line:
//     value0 
//     value1
//     value2
//     value3
// ---------------------------------------------------------------------------
void coindidence_3(string file_angles, string file_coindicence){
fstream data;
    data.open(file_angles, ios::in);

    if(data.good()){
        cout << "Acess to file has been granted!" << endl;
        
        string line;
        fstream file(file_coindicence, ios::out);

        while (!data.eof()) {
            getline(data, line);
            if (read_line(1, line.size(), line) == "001011"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "111000"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
        }
    } 
    else {
        cout << "Access to file was denided!" << endl;
    }

    data.close();
}

// ---------------------------------------------------------------------------
// process()
// ---------------------------------------------------------------------------
// Reads an input text file containing integer values, one per line. Each value
// is transformed by rounding it to the nearest multiple of (3 * multiplier).
// 
// The step size used for rounding is:
//     step = 3 * multiplier
//
// The transformation applied is:
//     output = step * floor((input + step/2) / step)
// where integer division is used. This results in rounding each number to the
// nearest multiple of the step size.
//
// Examples for multiplier = 1  (step = 3):
//   input: 1  → output: 0
//   input: 2  → output: 3
//   input: 4  → output: 3
//   input: 5  → output: 6
//
// Examples for multiplier = 10 (step = 30):
//   input: 17  → output: 0
//   input: 29  → output: 30
//   input: 44  → output: 30
//   input: 46  → output: 60
//
// Parameters:
//   - inputFileName:   path to the input file containing integers
//   - outputFileName:  path to the output file where results will be saved
//   - multiplier:      scaling factor; the rounding step becomes 3 * multiplier
//
// Side effects:
//   - Writes processed values to the output file
//   - Prints status and error messages to stdout/stderr
//
// The function completes after all integers from the input file are processed.
// ---------------------------------------------------------------------------
void process(string inputFileName, string outputFileName, int multiplier){
    // 1. Opening the input file
    ifstream inputFile(inputFileName); 
    if (!inputFile.is_open()) {
        cerr << "Error: Cannot open file " << inputFileName << " for reading." << endl;
    }

    // 2. Opening the output file
    ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        cerr << "Error: Cannot open file " << outputFileName << " for writing." << endl;
        inputFile.close();
    }

    int inputValue;
    int outputValue;

    // Step size computed as 3 * multiplier
    int step = 3 * multiplier;

    // 3. Reading, processing, and writing
    while (inputFile >> inputValue) {
        // Transformation logic: step * floor((inputValue + step/2) / step)
        // integer division is used
        outputValue = step * ((inputValue + step/2) / step);

        // Writing to file
        outputFile << outputValue << endl;
    }

    // 4. Closing the files
    inputFile.close();
    outputFile.close();

    cout << "Processing completed with step " << step 
         << ". Results saved to file " << outputFileName << endl;
}

// ---------------------------------------------------------------------------
// coincidence_123_124()
// ---------------------------------------------------------------------------
// Reads an input file (file_angles) and filters the lines based on field 1.
// The function is trying to obtain specific coincidences for analysis. Here
// coincidences between detectors 1, 2, 3 and 1, 2, 4 are gathered.
//
// Example output line:
//     value0 
//     value1
//     value2
//     value3
// ---------------------------------------------------------------------------
void coindidence_123_124(string file_angles, string file_coindicence){
    fstream data;
    data.open(file_angles, ios::in);

    if(data.good()){
        cout << "Acess to file has been granted!" << endl;
        
        string line;
        fstream file(file_coindicence, ios::out);

        while (!data.eof()) {
            getline(data, line);

            if (read_line(1, line.size(), line) == "001011"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "010101"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
        }
    } 
    else {
        cout << "Access to file was denided!" << endl;
    }
    data.close();
}

// ---------------------------------------------------------------------------
// coincidence_134_234()
// ---------------------------------------------------------------------------
// Reads an input file (file_angles) and filters the lines based on field 1.
// The function is trying to obtain specific coincidences for analysis. Here
// coincidences between detectors 1, 3, 4 and 2, 3, 4 are gathered.
//
// Example output line:
//     value0 
//     value1
//     value2
//     value3
// ---------------------------------------------------------------------------
void coindidence_134_234(string file_angles, string file_coindicence){
    fstream data;
    data.open(file_angles, ios::in);

    if(data.good()){
        cout << "Acess to file has been granted!" << endl;
        
        string line;
        fstream file(file_coindicence, ios::out);

        while (!data.eof()) {
            getline(data, line);

            if (read_line(1, line.size(), line) == "100110"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
            else if (read_line(1, line.size(), line) == "111000"){
            //file << read_line(0, line.size(), line) << ";";
            //file << read_line(1, line.size(), line) << ";";
            //file << read_line(2, line.size(), line) << ";";
            file << read_line(3, line.size(), line) << endl;
            }
        }
    } 
    else {
        cout << "Access to file was denided!" << endl;
    }
    data.close();
}
