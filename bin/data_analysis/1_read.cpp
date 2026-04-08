#include <iostream>
#include <fstream>
#include <string>
#include <conio.h>
#include <vector>
#include <filesystem>

using namespace std;

void folder_managment(string folder);
string read_line(int position, int size, string line);
int read_int(int position, int line_size, string line);
string read_binary_string(int position, int line_size, string line);
int read_bit(int position, int line_size, string line);
void time_to_angle(const string& file_angles);
void delete_wrong_angles(string file_angles, string file_angles_processed);
void one_direction(string file_data, string file_angles, bool what_direction);
void output_data(string file_dane, string file_angles);
void coindidence_4(string file_angles, string file_coindicence);
void coindidence_3_4(string file_angles, string file_coindicence);
void coindidence_3_4f(string file_angles, string file_coindicence);
void coindidence_3f(string file_angles, string file_coindicence);
void coindidence_3(string file_angles, string file_coindicence);
void process(string inputFileName, string outputFileName, int multiplier);
void coindidence_134_234(string file_angles, string file_coindicence);
void coindidence_123_124(string file_angles, string file_coindicence);
long long read_long_long(int position, int line_size, string line);
void calculate_time_per_angle(string file_angles, string file_time_output);

int main (){

    const int START = 1;   // from which data file to start
    const int END   = 10;   // where to finish

        // Lambda function for safe file copying, which checks if the destination 
        // file exists and removes it before copying to avoid MinGW errors.
        auto safe_copy = [](const string& src, const string& dst) {
            if (filesystem::exists(dst)) {
                filesystem::remove(dst);
            }
            filesystem::copy_file(src, dst);
        };

        for (int i = START; i <= END; ++i) {

            //this folder have to exist bc user is putting there data manually
            string data_file =      "D:/MuonDensityMapping/data/dane_" + to_string(i) + ".log";
            string angle_path =     "D:/MuonDensityMapping/output/single/dane_" + to_string(i);
            string raw_path =       "D:/MuonDensityMapping/output/single/dane_" + to_string(i) + "/raw/";
            string processed_path = "D:/MuonDensityMapping/output/single/dane_" + to_string(i) + "/processed/";
            
            //this folder doesnt have to exist so first to check
            folder_managment(angle_path);
            folder_managment(processed_path);

            // 1. excavate only nessessery data
            output_data(data_file, angle_path + "/angles_raw_" + to_string(i) + ".txt");

            // 2. take excaveted data and assign to each line angles for which detection was done
            time_to_angle(angle_path + "/angles_raw_" + to_string(i) + ".txt");

            // 3. delete the lines where the angles are below 0 or above 180, which are not physical
            //    and then save the file with only those lines where angles are between 0 and 180
            delete_wrong_angles(angle_path + "/angles_raw_" + to_string(i) + ".txt", angle_path + "/angles_processed_" + to_string(i) + ".txt");
            
            
            // 4. create file with time stamp on each angle
            string time_output_file = processed_path + "time_per_angle_" + to_string(i) + ".txt";
            calculate_time_per_angle(angle_path + "/angles_processed_" + to_string(i) + ".txt", time_output_file);

            // 5. Test - if we want to process only one direction of rotation, for example only increasing or only decreasing.
            string angles_file_true = angle_path + "/angles_processed_true" + to_string(i) + ".txt";
            string angles_file_false = angle_path + "/angles_processed_false" + to_string(i) + ".txt";

            one_direction(angle_path + "/angles_processed_" + to_string(i) + ".txt", angle_path + "/angles_processed_true" + to_string(i) + ".txt", true);
            one_direction(angle_path + "/angles_processed_" + to_string(i) + ".txt", angle_path + "/angles_processed_false" + to_string(i) + ".txt", false);
            
            // Coincidence 3_4f - N -> S
            string coin34f_proc_true = processed_path + "coin_" + to_string(i) + "_true.txt";
            coindidence_3_4f(angles_file_true, coin34f_proc_true);

            // Coincidence 3_4f - S -> N
            string coin34f_proc_false = processed_path + "coin_" + to_string(i) + "_false.txt";
            coindidence_3_4f(angles_file_false, coin34f_proc_false);
            
            // 6. take excaveted data with angles and search for needed coincidences
            string angles_file = angle_path + "/angles_processed_" + to_string(i) + ".txt";
            
            // Coincidence 3
            string coin3_proc = processed_path + "coin_" + to_string(i) + "_3.txt";
            coindidence_3(angles_file, coin3_proc);

            // Coincidence 3f
            string coin3f_proc = processed_path + "coin_" + to_string(i) + "_3f.txt";
            coindidence_3f(angles_file, coin3f_proc);

            // Coincidence 4
            string coin4_proc = processed_path + "coin_" + to_string(i) + "_4.txt";
            coindidence_4(angles_file, coin4_proc);

            // Coincidence 3_4
            string coin34_proc = processed_path + "coin_" + to_string(i) + "_34.txt";
            coindidence_3_4(angles_file, coin34_proc);

            // Coincidence 3_4f
            string coin34f_proc = processed_path + "coin_" + to_string(i) + "_34f.txt";
            coindidence_3_4f(angles_file, coin34f_proc);
        }

    
    //additional coincidances if needed
        //coindidence_134_234("../Output/all/angles_all.txt", "../Output/all/coin_134_234.txt");
        //coindidence_123_124("../Output/all/angles_all.txt", "../Output/all/coin_123_124.txt");

    //process data for histograms (each bar have to have 3 degrees)
        //process("../Output/coin_x_x.txt", "../Output/processed_x_x.txt", 3);
        //process("../Output/coin_x_x.txt", "../Output/processed_x_x.txt", 3);


    return 0;
}

// ---------------------------------------------------------------------------
// Manages with foldres in order to automatise the process
// folder  – indicates the places of folder
// The function checks whatever the designated folder have to be created
// or its already existing
// ---------------------------------------------------------------------------
void folder_managment(string folder){
    if (filesystem::exists(folder)) {
        if (filesystem::is_directory(folder)) {
            cout << "Folder '" << folder << "' already exist." << endl;
        } else {
            cout << "Error - name conflict: File '" << folder << "'exist, cannot make folder" << endl;
        }
        } else {
            // If the folder does not exist, attempt to create it and report success or failure.
            if (filesystem::create_directory(folder)) {
                cout << "Succesufully created file: '" << folder << "'" << endl;
            } else {
                cout << "There was an error creating file" << endl;
            }
        }
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

    // --- error time test ---
    int t_err = read_int(0, line.size(), line);
    int t1_err = read_int(0, line.size(), line);
    int delta_t;
    int t_prev = t0;

    int angle_inc = 1;     // for v3 = 1 - increasing to 180
    int angle_dec = 180;   // for v3 = 0 - decreasing to 0

    // --- Angle of the first line ---
    if (v3_prev == 1)
        lines[0] += ";" + to_string(angle_inc);
    else
        lines[0] += ";" + to_string(angle_dec);

    // ---- Get next lines ----
    while (getline(file, line)) {
        int t1 = read_int(0, line.size(), line);
        int dt = t1 - t0;
        int v3 = read_bit(2, line.size(), line);
        t1_err = t_prev;                         // --error time test ---

        // --- If the direction changes → reset the counter ---
        if (v3 != v3_prev) {
            angle_inc = 1;
            angle_dec = 180;
            delta_t = t1_err - t_err - 32400000; // --error time test ---
            t0 = t1;  // reference time reset                          
            t_err = t1;                          // --error time test ---
            //printf("%d \n", delta_t);          // --error time test ---
            if (delta_t < -269071 || delta_t > -200293){
                //printf("delta %d \n", delta_t);
                printf("error line - time: %d \n", t0);
            }
        }

        // --- Update every 180,000 ms = 1 degree change---
        if (dt >= 180000) {
            t0 = t1 - dt + 180000;
            if (v3 == 1)
                angle_inc++;
            else
                angle_dec--;
        }

        // --- Adding an angle to the end of a line ---
        if (v3 == 1)
            line += ";" + to_string(angle_inc);
        else
            line += ";" + to_string(angle_dec);

        lines.push_back(line);

        v3_prev = v3; // remembering the direction for the next iteration
        t_prev = t1;                             // --error time test ---
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
// deletes the lines where the angles are below 0 or above 180, which are not physical, 
//and then save the file with only those lines where angles are between 0 and 180
void delete_wrong_angles(string file_angles, string file_angles_processed){
    cout << "delete_wrong_angles() function has been activated." << endl;

    ifstream file(file_angles);
    if (!file.good()) {
        cout << "Access to the file " << file_angles << " has been denied!" << endl;
        return;
    }
    cout << "Access to the file " << file_angles << " has been granted!" << endl;

    vector<string> valid_lines;   // This will hold all lines that are valid and can be saved to the output file
    vector<string> current_block; // This is a temporary storage for lines belonging to the current block of data 
                                  // with the same direction. We will only move these lines to valid_lines if the entire block is valid.
    
    int current_direction = -1;   // Flag for tracking the current direction (-1 at start)
    bool is_block_valid = true;   // We assume the new data block is valid until proven otherwise

    string line;
    while (getline(file, line)) {
        // Extract the direction and angle from the line
        int direction = read_bit(2, line.size(), line); 
        int angle = read_int(3, line.size(), line);

        // Set the initial direction for the first line
        if (current_direction == -1) {
            current_direction = direction;
        }

        // Check if the direction has changed (e.g., from 1 to 0)
        if (direction != current_direction) {
            
            // Before we start processing the new block, we need to check if the previous block was valid. 
            // If it was, we move all lines from current_block to valid_lines. If it wasn't, we discard the 
            // entire block and print a message.
            if (is_block_valid) {
                // If the block is valid, we add all lines from current_block to valid_lines
                valid_lines.insert(valid_lines.end(), current_block.begin(), current_block.end());
            } else {
                cout << "Deleted invalid block for direction: " << current_direction << endl;
            }

            // Now we can start processing the new block. We clear the current_block, update the current_direction, 
            //and reset the is_block_valid flag for the new block.
            current_block.clear();          // Clear the current block for the new direction
            current_direction = direction;  // Update the current direction to the new one
            is_block_valid = true;          // Assume the new block is valid until we find an invalid angle in it
        }

        // Now we are processing a line that belongs to the current block. We add it to the current_block vector. 
        // However, if the angle is invalid (below 0 or above 180), we mark the entire block as invalid by setting is_block_valid to false. 
        // This means that when we encounter a direction change, we will discard all lines in this block instead of just the single invalid line.
        current_block.push_back(line);

        // Check if the angle is valid. If it's not, we mark the entire block as invalid. We will handle the actual discarding 
        // of lines when we encounter a direction change or at the end of the file.
        if (angle < 0 || angle > 180) {
            is_block_valid = false;
        }
    }

    // After we finish reading the file, we need to check the last block of data. If it was valid, we add its lines to valid_lines. 
    //If it wasn't, we discard it and print a message.
    if (!current_block.empty()) {
        if (is_block_valid) {
            valid_lines.insert(valid_lines.end(), current_block.begin(), current_block.end());
        } else {
            cout << "Deleted invalid block at the end of the file." << endl;
        }
    }

    file.close();
    cout << "The file " << file_angles <<" has been closed" << endl;

    // write vector to file
    ofstream out(file_angles_processed, ios::trunc); // delete the old file -> print the values from vector
    for (auto& l : valid_lines) {
        out << l << "\n";
    }
    out.close();

    cout << "Processing of the file " << file_angles_processed <<" has been completed and saved." << endl;
}

// ---------------------------------------------------------------------------
// Bool true means that the direction is increasing to 180 (North to South), false - decreasing to 0 (South to North). 
// The function is used when we want to process only one direction of rotation, for example only increasing or only decreasing. 
void one_direction(string file_data, string file_angles, bool what_direction){
    cout << "one_direction() function has been activated." << endl;

    ifstream file(file_data);
    if (!file.good()) {
        cout << "Access to the file " << file_data << " has been denied!" << endl;
        return;
    }

    cout << "Access to the file " << file_data << " has been granted!" << endl;

    vector<string> lines;
    string line;

    while (getline(file, line)) {
        int v3 = read_bit(2, line.size(), line);
        if ((what_direction && v3 == 1) || (!what_direction && v3 == 0)) {
            lines.push_back(line);
        }
    }

    file.close();
    cout << "The file " << file_data <<" has been closed" << endl;

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

// ---------------------------------------------------------------------------
// A helper function for safely reading large numbers (e.g. time in ms)
// ---------------------------------------------------------------------------
long long read_long_long(int position, int line_size, string line) {
    string cell = read_line(position, line_size, line);
    if(cell.empty()) return 0;
    return stoll(cell);
}

// ----------------------------------------------------------------------------
// calculate_time_per_angle()
// ------- --------------------------------------------------------------------
// Iterates through a file containing assigned angles and calculates the actual time 
// the detector spends at a given angle. It sums the time differences (delta t) 
// and saves the result to a new file in the format: Angle;Time_in_Minutes.
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// calculate_time_per_angle()
// ---------------------------------------------------------------------------
// Iterates through the file with assigned angles and calculates the actual time 
// the detector spends at a given angle, and counts the number of 34f coincidences.
// Saves the result to a new file in the format: Angle;Time_in_Minutes;Number_of_Counts.
// ---------------------------------------------------------------------------
void calculate_time_per_angle(string file_angles, string file_time_output) {
    cout << "calculate_time_per_angle() function has been activated." << endl;

    ifstream file(file_angles);
    if (!file.good()) {
        cout << "Access to the file " << file_angles << " has been denied!" << endl;
        return;
    }

    // arrays storing the total time (ms) and the number of counts for each angle (0–180)
    long long time_per_angle_ms[181] = {0}; 
    long long count_34f_per_angle[181] = {0}; // <--- NOWA TABLICA
    
    string line;
    if (!getline(file, line)) return;

    // Initialising variables based on the first line
    long long t_enter = read_long_long(0, line.size(), line);
    string coin_current = read_line(1, line.size(), line); // POBIERANIE ZNAKÓW KOINCYDENCJI
    int current_dir = read_bit(2, line.size(), line);
    int current_angle = read_int(3, line.size(), line);
    long long t_prev = t_enter;

    // Checking the 34f coincidence for the first line
    if (current_angle >= 0 && current_angle <= 180) {
        if (coin_current == "000000" || 
            coin_current == "001011" || 
            coin_current == "111000" || 
            coin_current == "010101" || 
            coin_current == "100110") {
            count_34f_per_angle[current_angle]++;
        }
    }

    while (getline(file, line)) {
        long long t_current = read_long_long(0, line.size(), line);
        string coin_string = read_line(1, line.size(), line);
        int dir_current = read_bit(2, line.size(), line);
        int angle_current = read_int(3, line.size(), line);

        // Ignore angle errors
        if (angle_current < 0 || angle_current > 180) continue;

        if (coin_string == "000000" || 
            coin_string == "001011" || 
            coin_string == "111000" || 
            coin_string == "010101" || 
            coin_string == "100110") {
            count_34f_per_angle[angle_current]++;
        }

        // 1. A change of direction or a massive time gap or restart
        if (dir_current != current_dir || (t_current - t_prev) > 600000 || t_current < t_prev) {
            if (current_angle >= 0 && current_angle <= 180) {
                // Dodaj czas tylko jeśli nie ma głupot przed restartem
                if (t_prev >= t_enter) {
                    time_per_angle_ms[current_angle] += (t_prev - t_enter);
                }
            }
            current_dir = dir_current;
            current_angle = angle_current;
            t_enter = t_current; // Prawidłowy reset stopera
        } 
        // 2. Normal transition to the next angle
        else if (angle_current != current_angle) {
            if (current_angle >= 0 && current_angle <= 180) {
                if (t_current >= t_enter) {
                    time_per_angle_ms[current_angle] += (t_current - t_enter);
                }
            }
            current_angle = angle_current;
            t_enter = t_current;
        }
        
        t_prev = t_current;
    }
    
    // Adding the final segment after exiting the loop
    if (current_angle >= 0 && current_angle <= 180) {
        time_per_angle_ms[current_angle] += (t_prev - t_enter);
    }

    file.close();

    // Save the total times and 34f counts to a file
    ofstream out(file_time_output, ios::trunc);
    out << "Angle;Time_Minutes;Counts\n"; //
    
    for (int i = 0; i <= 180; ++i) {
        if (time_per_angle_ms[i] > 0 || count_34f_per_angle[i] > 0) { 
            double minutes = time_per_angle_ms[i] / 60000.0;
            out << i << ";" << minutes << ";" << count_34f_per_angle[i] << "\n"; 
        }
    }
    out.close();

    cout << "Exposure times and 34f counts successfully saved to: " << file_time_output << endl;
}