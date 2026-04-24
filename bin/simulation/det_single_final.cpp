#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>

// === Required to load PNG images in C++ ===
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace fs = std::filesystem;

// ================= SIMULATION PARAMETERS =================
const std::string IMAGE_PATH = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_data/converted_black_white_glory_hole.png"; 

const std::string OUTPUT_FOLDER = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/MC_single/";
const std::string OUTPUT_FILE = OUTPUT_FOLDER + "det_4_with_hole.txt";

const double PIXEL_SIZE = 0.13333;

const double DETECTOR_X = 963.0;
const double DETECTOR_Y = 1247.0;

const double MAX_RADIUS_M = 48.0;   // R0 dla wygaszania
const double VIEW_RADIUS_M = 100.0; // Maksymalny zasięg promieni

const double W = 0.20;  
const double H = 0.10;  
const double L = 1.0;   

const double SURFACE_Z = 15.0;       
const double STEP_SIZE_M = 0.1;
const double ROCK_DENSITY = 2.65; 

// =======================================================
const double VIRTUAL_TIME_PER_ANGLE_MIN = 500.0; 
const int VIRTUAL_TOTAL_EVENTS = 10000;        
// =======================================================

// 0 = Hard Cutoff, 1 = Exponential, 2 = Gaussian
const int ATTENUATION_MODE = 0; 

const int ASCII_BIN_SIZE = 5;  

const double A1 = -11.22, A2 = -0.00262, A3 = -14.10, A4 = -0.001213;

double crouch_mu(double h_mwe) {
    return std::exp(A1 + A2 * h_mwe) + std::exp(A3 + A4 * h_mwe);
}

int main() {
    std::string mode_str = (ATTENUATION_MODE == 0) ? "Hard Cutoff" : (ATTENUATION_MODE == 1) ? "Exponential" : "Gaussian";
    std::cout << "============================================================\n";
    std::cout << "--- PURE DETERMINISTIC SIMULATOR | R0 = " << MAX_RADIUS_M << "m | Mode: " << mode_str << " ---\n";
    std::cout << "============================================================\n";

    std::cout << "1. Loading theoretical cave map (test_1.png)...\n";
    int width, height, channels;
    
    std::vector<double> time_data(180, VIRTUAL_TIME_PER_ANGLE_MIN);

    unsigned char *img_data = stbi_load(IMAGE_PATH.c_str(), &width, &height, &channels, 1);
    if (!img_data) {
        std::cerr << "[ERROR] Cannot load file " << IMAGE_PATH << ".\n";
        return 1;
    }

    double max_az = std::atan(W / L);
    double max_el = std::atan(H / L);
    
    std::vector<double> az_valid, el_valid, weights_valid;
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

    // ================= PHASE 2: NORMALIZATION (DETERMINISTIC) =================
    std::cout << "\n3. Phase 2/2: Normalizing theoretical events...\n";

    std::vector<double> expected_total_per_bin(180, 0.0);
    for (int i = 0; i < 180; ++i) {
        expected_total_per_bin[i] = raw_flux_360[i] + raw_flux_360[i + 180];
    }

    double sum_time_flux = 0.0;
    for (int i = 0; i < 180; ++i) {
        sum_time_flux += expected_total_per_bin[i] * time_data[i];
    }

    double total_simulated_events = 0;
    std::vector<double> det_counts(181, 0.0); 

    for (int i = 0; i < 180; ++i) {
        expected_total_per_bin[i] = (expected_total_per_bin[i] * time_data[i] / sum_time_flux) * VIRTUAL_TOTAL_EVENTS;
        det_counts[i + 1] = expected_total_per_bin[i];
        total_simulated_events += det_counts[i + 1];
    }

    if (!fs::exists(OUTPUT_FOLDER)) fs::create_directories(OUTPUT_FOLDER);
    std::ofstream out(OUTPUT_FILE);
    out << "Angle;Time_Minutes;Counts_34f\n";
    for (int i = 0; i <= 180; ++i) {
        double t_min = (i >= 1 && i <= 180) ? time_data[i - 1] : 0.0;
        out << i << ";" << t_min << ";" << std::fixed << std::setprecision(6) << det_counts[i] << "\n";
    }
    out.close();

    // ================= DYNAMIC ASCII HISTOGRAM =================
    std::cout << "\n================== THEORETICAL HISTOGRAM ==================\n";
    int num_bins_ascii = (180 + ASCII_BIN_SIZE - 1) / ASCII_BIN_SIZE;
    std::vector<double> binned_counts_ascii(num_bins_ascii, 0.0);

    for (int i = 0; i < 180; ++i) {
        binned_counts_ascii[i / ASCII_BIN_SIZE] += det_counts[i + 1];
    }

    // Szukamy minimum i maksimum do auto-zooma
    double min_bin = binned_counts_ascii[0];
    double max_bin = binned_counts_ascii[0];
    for (double count : binned_counts_ascii) {
        if (count < min_bin) min_bin = count;
        if (count > max_bin) max_bin = count;
    }

    const int MAX_STARS = 50;
    double range = max_bin - min_bin;
    if (range < 1e-6) range = 1.0; // Zabezpieczenie

    for (int b = 0; b < num_bins_ascii; ++b) {
        int start_deg = b * ASCII_BIN_SIZE + 1;
        int end_deg = std::min(start_deg + ASCII_BIN_SIZE - 1, 180);
        double count = binned_counts_ascii[b];
        
        // Obliczamy gwiazdki po odjęciu tła (minimum)
        int num_stars = std::round(((count - min_bin) / range) * MAX_STARS);
        
        std::cout << std::setw(3) << start_deg << "-" << std::setw(3) << end_deg << " deg | " 
                  << std::fixed << std::setprecision(1) << std::setw(8) << count << " | " 
                  << std::string(num_stars, '*') << "\n";
    }
    std::cout << "===========================================================\n";

    return 0;
}