#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <chrono>
#include <filesystem>

// === Required to load PNG images in C++ ===
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace fs = std::filesystem;

// ================= SIMULATION PARAMETERS =================
const std::string IMAGE_PATH = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_data/converted_black_white.png"; 
const std::string DATA_PATH = "/home/kacper/MuonDensityMapping/bin/data_analysis/output/selected_rot_series/series_900.txt"; 

// --- NOWA ŚCIEŻKA WYJŚCIOWA ---
const std::string OUTPUT_FOLDER = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/900/MC/";
const std::string OUTPUT_FILE = OUTPUT_FOLDER + "pos_2_900_bkg_instant.txt";

const double PIXEL_SIZE = 0.13333;

//const double DETECTOR_X = 1000.0;
//const double DETECTOR_Y = 850.0;

const double DETECTOR_X = 916.0;
const double DETECTOR_Y = 1163.0;

const double MAX_RADIUS_M = 48.0;   // R0 dla wygaszania
const double VIEW_RADIUS_M = 100.0; // Maksymalny zasięg promieni (wcześniej skanowało tylko do MAX_RADIUS_M)

const double W = 0.20;  
const double H = 0.10;  
const double L = 1.0;   

const double SURFACE_Z = 15.0;       
const double STEP_SIZE_M = 0.1;
const double ROCK_DENSITY = 2.65; 

// =======================================================
// --- PRZEŁĄCZNIK MODELU WYGASZANIA PUSTEJ PRZESTRZENI ---
// 0 = Hard Cutoff (Natychmiastowe obcięcie po R0)
// 1 = Exponential (Zanik eksponencjalny e^(-r/R0))
// 2 = Gaussian (Zanik gaussowski e^(-(r/R0)^2))
// =======================================================
const int ATTENUATION_MODE = 0; 

// ---- MONTE CARLO SCHEME ----
const int BACKGROUND_NOISE_PER_BIN = 10;  
const double MC_STAT_MULTIPLIER = 1.0; // Ustawione na 1.0 dla naturalnej symulacji, chyba że chcesz więcej statystyki

const int ASCII_BIN_SIZE = 5;  

const double A1 = -11.22, A2 = -0.00262, A3 = -14.10, A4 = -0.001213;

double crouch_mu(double h_mwe) {
    return std::exp(A1 + A2 * h_mwe) + std::exp(A3 + A4 * h_mwe);
}

void load_data_and_time(const std::string& filename, std::vector<int>& counts, std::vector<double>& times, int& total_events) {
    counts.assign(180, 0);
    times.assign(180, 0.0);
    total_events = 0;
    
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "[ERROR] Missing data file: " << filename << "\n";
        exit(1);
    }
    
    std::string line;
    std::getline(file, line); 
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string item;
        
        std::getline(ss, item, ';'); int angle = std::stoi(item);
        std::getline(ss, item, ';'); double time = std::stod(item);
        std::getline(ss, item, ';'); int count = std::stoi(item);
        
        int idx = angle - 1; 
        if (idx >= 0 && idx < 180) {
            counts[idx] = count;
            times[idx] = time;
            total_events += count;
        }
    }
}

int main() {
    std::string mode_str = (ATTENUATION_MODE == 0) ? "Hard Cutoff" : (ATTENUATION_MODE == 1) ? "Exponential" : "Gaussian";
    std::cout << "============================================================\n";
    std::cout << "--- MONTE CARLO SIMULATOR | R0 = " << MAX_RADIUS_M << "m | Mode: " << mode_str << " ---\n";
    std::cout << "============================================================\n";

    std::cout << "1. Loading cave map and real data...\n";
    int width, height, channels, total_observed = 0;
    
    std::vector<int> real_counts;
    std::vector<double> time_data;
    load_data_and_time(DATA_PATH, real_counts, time_data, total_observed);

    unsigned char *img_data = stbi_load(IMAGE_PATH.c_str(), &width, &height, &channels, 1);
    if (!img_data) {
        std::cerr << "[ERROR] Cannot load file " << IMAGE_PATH << ".\n";
        return 1;
    }

    double max_az = std::atan(W / L);
    double max_el = std::atan(H / L);
    
    std::vector<double> az_valid, el_valid, weights_valid;

    // POPRAWKA: Stała rozdzielczość próbkowania oparta o pełen widok
    double step_az_el = PIXEL_SIZE / VIEW_RADIUS_M;
    
    for (double el = 0.0; el <= max_el + 1e-6; el += step_az_el) {
        for (double az = -max_az; az <= max_az + 1e-6; az += step_az_el) {
            double w_h = 1.0 - (L / W) * std::abs(std::tan(az));
            double w_v = 1.0 - (L / H) * std::abs(std::tan(el));
            double weight = w_h * w_v;
            
            if (weight > 0.0) {
                az_valid.push_back(az);
                el_valid.push_back(el);
                weights_valid.push_back(weight);
            }
        }
    }

    std::cout << "2. Phase 1/2: Starting analytical space scanning (0-179 deg)...\n";
    
    std::vector<double> raw_flux_360(360, 0.0);
    int bar_width = 40; 
    int max_map_steps = static_cast<int>(std::round(VIEW_RADIUS_M / STEP_SIZE_M));

    for (int angle_deg = 0; angle_deg < 180; ++angle_deg) {
        for (int side = 0; side < 2; ++side) {
            double rad = (angle_deg + side * 180) * (M_PI / 180.0);
            double eff_intensity_num = 0.0, weight_sum = 0.0;

            for (size_t i = 0; i < az_valid.size(); ++i) {
                double ray_az = rad + az_valid[i];
                double dir_x = std::sin(ray_az);
                double dir_y = -std::cos(ray_az);
                double tan_el = std::tan(el_valid[i]);
                
                double horiz_rock_m = 0.0;

                for (int s = 0; s <= max_map_steps; ++s) {
                    double r = s * STEP_SIZE_M;
                    int xi = static_cast<int>(std::round(DETECTOR_X + (r / PIXEL_SIZE) * dir_x));
                    int yi = static_cast<int>(std::round(DETECTOR_Y + (r / PIXEL_SIZE) * dir_y));
                    double z = r * tan_el;

                    bool is_rock = true;
                    if (xi >= 0 && xi < width && yi >= 0 && yi < height) {
                        is_rock = (img_data[yi * width + xi] < 128 && z <= SURFACE_Z);
                    }

                    if (is_rock) {
                        horiz_rock_m += STEP_SIZE_M;
                    } else {
                        // APLIKOWANIE MODELU WYGASZANIA
                        double visibility_weight = 0.0;
                        if (MAX_RADIUS_M > 0.0) {
                            if (ATTENUATION_MODE == 0) {
                                visibility_weight = (r <= MAX_RADIUS_M) ? 1.0 : 0.0;
                            } else if (ATTENUATION_MODE == 1) {
                                visibility_weight = std::exp(-r / MAX_RADIUS_M);
                            } else if (ATTENUATION_MODE == 2) {
                                visibility_weight = std::exp(-(r * r) / (MAX_RADIUS_M * MAX_RADIUS_M));
                            }
                        }
                        horiz_rock_m += STEP_SIZE_M * (1.0 - visibility_weight);
                    }
                }

                // POPRAWKA: Usunięcie cos^2(theta), poprawny slant depth
                double h_mwe = (horiz_rock_m / std::cos(el_valid[i])) * ROCK_DENSITY;
                double intensity = crouch_mu(h_mwe);

                eff_intensity_num += intensity * weights_valid[i];
                weight_sum += weights_valid[i];
            }

            raw_flux_360[angle_deg + side * 180] = eff_intensity_num / weight_sum;
        }

        float progress = (float)(angle_deg + 1) / 180.0f;
        int pos = bar_width * progress;
        std::cout << "\r[";
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << std::setw(3) << int(progress * 100.0) << "% " << std::flush;
    }
    
    std::cout << "\nRay-tracing finished successfully.\n";
    stbi_image_free(img_data);

    // ================= PHASE 2: NORMALIZATION AND MONTE CARLO =================
    std::cout << "\n3. Phase 2/2: Generating MC events (Poisson distribution)...\n";

    std::vector<double> expected_total_per_bin(180, 0.0);

    // POPRAWKA: Brak sztucznego rozmywania Bartlettem, sumujemy po prostu strony 0-180 i 180-360
    for (int i = 0; i < 180; ++i) {
        expected_total_per_bin[i] = raw_flux_360[i] + raw_flux_360[i + 180];
    }

    double sum_time_flux = 0.0;
    for (int i = 0; i < 180; ++i) {
        sum_time_flux += expected_total_per_bin[i] * time_data[i];
    }

    for (int i = 0; i < 180; ++i) {
        expected_total_per_bin[i] = (expected_total_per_bin[i] * time_data[i] / sum_time_flux) * (total_observed * MC_STAT_MULTIPLIER);
        expected_total_per_bin[i] += (BACKGROUND_NOISE_PER_BIN * MC_STAT_MULTIPLIER);
    }

    unsigned int master_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(master_seed); 
    
    std::vector<long long> mc_counts(181, 0); 
    long long total_simulated_events = 0;

    // Losowanie zdarzeń (Poisson)
    for (int i = 0; i < 180; ++i) {
        double mean = std::max(1e-6, expected_total_per_bin[i]);
        std::poisson_distribution<> d(mean);
        // Indeksujemy od 1 do 180, żeby pasowało do kątów w pliku
        mc_counts[i + 1] = d(gen);
        total_simulated_events += mc_counts[i + 1];
    }

    // ================= PLIK WYJŚCIOWY =================
    if (!fs::exists(OUTPUT_FOLDER)) {
        fs::create_directories(OUTPUT_FOLDER);
    }

    std::ofstream out(OUTPUT_FILE);
    out << "Angle;Time_Minutes;Counts_34f\n";

    // Pętla po kątach 0 - 180 (0 jest z reguły puste, ale trzymamy format)
    for (int i = 0; i <= 180; ++i) {
        double t_min = 0.0;
        if (i >= 1 && i <= 180) {
            t_min = time_data[i - 1]; // W time_data indeksy lecą 0-179 dla kątów 1-180
        }
        out << i << ";" << t_min << ";" << mc_counts[i] << "\n";
    }
    out.close();

    std::cout << "\nSaved formatted output file: " << OUTPUT_FILE << "\n";
    std::cout << "Total generated MC events: " << total_simulated_events << "\n";

    // ================= ASCII HISTOGRAM =================
    std::cout << "\n================== MONTE CARLO HISTOGRAM ==================\n";
    int num_bins_ascii = (180 + ASCII_BIN_SIZE - 1) / ASCII_BIN_SIZE;
    std::vector<int> binned_counts_ascii(num_bins_ascii, 0);

    for (int i = 0; i < 180; ++i) {
        binned_counts_ascii[i / ASCII_BIN_SIZE] += mc_counts[i + 1];
    }

    int max_bin_count_ascii = 0;
    for (int count : binned_counts_ascii) {
        if (count > max_bin_count_ascii) max_bin_count_ascii = count;
    }

    const int MAX_STARS = 50;

    for (int b = 0; b < num_bins_ascii; ++b) {
        int start_deg = b * ASCII_BIN_SIZE + 1;
        int end_deg = std::min(start_deg + ASCII_BIN_SIZE - 1, 180);
        int count = binned_counts_ascii[b];
        int num_stars = max_bin_count_ascii > 0 ? (count * MAX_STARS / max_bin_count_ascii) : 0;
        
        std::cout << std::setw(3) << start_deg << "-" << std::setw(3) << end_deg << " deg | " 
                  << std::setw(5) << count << " | " << std::string(num_stars, '*') << "\n";
    }
    std::cout << "===========================================================\n";

    return 0;
}