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

// === Required to load PNG images in C++ ===
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// ================= SIMULATION PARAMETERS =================
//const std::string IMAGE_PATH = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_data/converted_black_white.png";
//const std::string IMAGE_PATH = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_data/converted_black_white_4.png"; //for background MC
const std::string IMAGE_PATH = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_data/converted_black_white_glory_hole.png"; //for hole 
const std::string DATA_PATH = "/home/kacper/MuonDensityMapping/output/all/zussamen_10/coin_time.txt"; 


const std::string FILE1_DETAILS = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/3_hole_sim_details_C.txt";
const std::string FILE2_MC_SUMMARY = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/3_hole_sim_summary_C.txt";
const std::string FILE3_MC_EVENTS = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/3_hole_sim_root_C.txt";

const double PIXEL_SIZE = 0.13333;  

const double DETECTOR_X = 916.0;      
const double DETECTOR_Y = 1163.0;     

const double MAX_RADIUS_M = 30.0;   
const double VIEW_RADIUS_M = 100.0; 

const double W = 0.20;  
const double H = 0.10;  
const double L = 1.0;   

const double CAVE_FLOOR_Z = -1.0;    
const double CAVE_CEILING_Z = 1.5;   
const double SURFACE_Z = 15.0;       

const double STEP_SIZE_M = 0.1;
const double ROCK_DENSITY = 2.65; 

// ---- MONTE CARLO SCHEME ----
const int BACKGROUND_NOISE_PER_BIN = 10;  

// ---- HISTOGRAM CONTROL ----
const int ASCII_BIN_SIZE = 5;  

const double A1 = -11.22, A2 = -0.00262, A3 = -14.10, A4 = -0.001213;

// ================= AUXILIARY STRUCTURES =================
struct RayData {
    double horiz_m;
    double vert_m;
    double rock_3d_m;
    double intensity;
};

struct AngleData {
    int angle_deg;
    RayData fwd;
    RayData bwd;
    double total_intensity;
};

// Analytical function (Crouch)
double crouch_mu(double h_mwe) {
    return std::exp(A1 + A2 * h_mwe) + std::exp(A3 + A4 * h_mwe);
}

// Function generating Bartlett window (consistent with np.bartlett)
std::vector<double> get_bartlett_kernel(int M) {
    std::vector<double> w(M);
    double sum = 0.0;
    for (int n = 0; n < M; ++n) {
        w[n] = (2.0 / (M - 1.0)) * ((M - 1.0) / 2.0 - std::abs(n - (M - 1.0) / 2.0));
        sum += w[n];
    }
    for (int n = 0; n < M; ++n) w[n] /= sum; // Normalization to 1.0
    return w;
}

// Function loading both counts and exposure time simultaneously
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
    std::getline(file, line); // Skip header (Angle;Time_Minutes;Counts_34f)
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string item;
        
        std::getline(ss, item, ';'); int angle = std::stoi(item);
        std::getline(ss, item, ';'); double time = std::stod(item);
        std::getline(ss, item, ';'); int count = std::stoi(item);
        
        int idx = angle - 1; // Map angles 1-180 to indices 0-179
        if (idx >= 0 && idx < 180) {
            counts[idx] = count;
            times[idx] = time;
            total_events += count;
        }
    }
}

int main() {
    std::cout << "1. Loading cave map and real data...\n";
    int width, height, channels, total_observed = 0;
    
    // Load real data and times
    std::vector<int> real_counts;
    std::vector<double> time_data;
    load_data_and_time(DATA_PATH, real_counts, time_data, total_observed);

    // Load image in forced 1-channel (grayscale) mode
    unsigned char *img_data = stbi_load(IMAGE_PATH.c_str(), &width, &height, &channels, 1);
    if (!img_data) {
        std::cerr << "[ERROR] Cannot load file " << IMAGE_PATH << ". Reason: " << stbi_failure_reason() << "\n";
        return 1;
    }

    // Calculating mesh grid and masks
    double max_az = std::atan(W / L);
    double max_el = std::atan(H / L);
    
    std::vector<double> az_valid;
    std::vector<double> el_valid;
    std::vector<double> weights_valid;

    double step_az_el = PIXEL_SIZE / MAX_RADIUS_M;
    
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

    auto scan_sector = [&](double central_azimuth_rad) -> RayData {
        double eff_intensity_num = 0.0;
        double eff_horiz_num = 0.0;
        double eff_vert_num = 0.0;
        double weight_sum = 0.0;

        for (size_t i = 0; i < az_valid.size(); ++i) {
            double ray_az = central_azimuth_rad + az_valid[i];
            double dir_x = std::sin(ray_az);
            double dir_y = -std::cos(ray_az);
            double tan_el = std::tan(el_valid[i]);
            
            double horiz_rock_m = 0.0;

            // Raycasting
            for (double r = 0.0; r <= MAX_RADIUS_M + 1e-6; r += STEP_SIZE_M) {
                double X = DETECTOR_X + (r / PIXEL_SIZE) * dir_x;
                double Y = DETECTOR_Y + (r / PIXEL_SIZE) * dir_y;
                double Z = r * tan_el;

                int X_idx = static_cast<int>(std::round(X));
                int Y_idx = static_cast<int>(std::round(Y));

                bool in_bounds = (X_idx >= 0 && X_idx < width && Y_idx >= 0 && Y_idx < height);
                int X_safe = std::clamp(X_idx, 0, width - 1);
                int Y_safe = std::clamp(Y_idx, 0, height - 1);

                unsigned char pixel_val = img_data[Y_safe * width + X_safe];
                
                bool is_white_pixel = (pixel_val >= 128) && in_bounds;
                bool is_rock = true;
                
                if (is_white_pixel) is_rock = false;
                if (Z > SURFACE_Z) is_rock = false;

                if (is_rock) horiz_rock_m += STEP_SIZE_M;
            }

            double vert_rock_m = horiz_rock_m * tan_el;
            double total_rock_3d = std::sqrt(horiz_rock_m * horiz_rock_m + vert_rock_m * vert_rock_m);
            double h_mwe = total_rock_3d * ROCK_DENSITY;

            double zenith_rad = (M_PI / 2.0) - el_valid[i];
            double initial_flux = std::cos(zenith_rad) * std::cos(zenith_rad);
            double intensity = initial_flux * crouch_mu(h_mwe);

            eff_intensity_num += intensity * weights_valid[i];
            eff_horiz_num += horiz_rock_m * weights_valid[i];
            eff_vert_num += vert_rock_m * weights_valid[i];
            weight_sum += weights_valid[i];
        }

        RayData res;
        res.horiz_m = eff_horiz_num / weight_sum;
        res.vert_m = eff_vert_num / weight_sum;
        res.rock_3d_m = std::sqrt(res.horiz_m * res.horiz_m + res.vert_m * res.vert_m);
        res.intensity = eff_intensity_num / weight_sum;
        return res;
    };

    std::cout << "2. Phase 1/2: Starting analytical space scanning (0-179 deg)...\n";
    
    std::vector<AngleData> sim_data;
    std::vector<double> raw_flux_360(360, 0.0);
    double sum_of_all_intensities = 0.0;

    int bar_width = 40; // Szerokość paska ładowania

    for (int angle_deg = 0; angle_deg < 180; ++angle_deg) {
        double fwd_rad = angle_deg * (M_PI / 180.0);
        double bwd_rad = (angle_deg + 180) * (M_PI / 180.0);

        RayData fwd = scan_sector(fwd_rad);
        RayData bwd = scan_sector(bwd_rad);

        raw_flux_360[angle_deg] = fwd.intensity;
        raw_flux_360[angle_deg + 180] = bwd.intensity;

        double total_i = fwd.intensity + bwd.intensity;
        sum_of_all_intensities += total_i;

        AngleData ad = { angle_deg, fwd, bwd, total_i };
        sim_data.push_back(ad);

        // --- DYNAMIC PROGRESS BAR ---
        float progress = (float)(angle_deg + 1) / 180.0f;
        int pos = bar_width * progress;
        
        std::cout << "\r[";
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << std::setw(3) << int(progress * 100.0) << "% "
                  << "| Angle: " << std::setw(3) << angle_deg 
                  << " | I: " << std::scientific << std::setprecision(2) << total_i << std::flush;
    }
    
    std::cout << "\nRay-tracing finished successfully.\n";

    // Zwolnienie pamięci obrazu
    stbi_image_free(img_data);

    // ================= PHASE 2: SMOOTHING AND MONTE CARLO =================
    std::cout << "\n3. Phase 2/2: Applying Bartlett filter and generating events (Real observed total: " 
              << total_observed << ", Bkg noise/bin: " << BACKGROUND_NOISE_PER_BIN << ")...\n";

    std::vector<double> kernel = get_bartlett_kernel(23);
    std::vector<double> expected_total_per_bin(180, 0.0);
    int half_k = kernel.size() / 2; // 11

    // Cyclic smoothing
    for (int i = 0; i < 360; ++i) {
        double smoothed = 0.0;
        for (int j = 0; j < (int)kernel.size(); ++j) {
            int idx = i + j - half_k;
            idx = (idx % 360 + 360) % 360; 
            smoothed += raw_flux_360[idx] * kernel[j];
        }
        if (i < 180) expected_total_per_bin[i] = smoothed;
        else expected_total_per_bin[i - 180] += smoothed; // Summing opposite directions
    }

    // SCALING TO ACTUAL EXPOSURE TIME
    double sum_time_flux = 0.0;
    for (int i = 0; i < 180; ++i) {
        sum_time_flux += expected_total_per_bin[i] * time_data[i];
    }

    for (int i = 0; i < 180; ++i) {
        expected_total_per_bin[i] = (expected_total_per_bin[i] * time_data[i] / sum_time_flux) * total_observed;
        expected_total_per_bin[i] += BACKGROUND_NOISE_PER_BIN;
    }

    // Monte Carlo
    unsigned int master_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(master_seed); 
    std::vector<int> mc_counts(180, 0);
    int total_simulated = 0;

    for (int i = 0; i < 180; ++i) {
        double mean = std::max(1e-6, expected_total_per_bin[i]);
        std::poisson_distribution<> d(mean);
        mc_counts[i] = d(gen);
        total_simulated += mc_counts[i];
    }

    // File Writing (Phase 1 & 2)
    std::ofstream f1(FILE1_DETAILS);
    std::ofstream f2(FILE2_MC_SUMMARY);
    
    f1 << "Angle_deg Fwd_Horiz_Rock_m Fwd_Vert_Rock_m Fwd_Rock_3D_m Bwd_Horiz_Rock_m Bwd_Vert_Rock_m Bwd_Rock_3D_m Total_Expected_Intensity_Raw\n";
    f2 << "Angle_deg MC_Simulated_Counts\n";

    for (int ang = 0; ang < 180; ++ang) {
        AngleData& d = sim_data[ang];
        f1 << std::fixed << std::setprecision(2)
           << ang << " " 
           << d.fwd.horiz_m << " " << d.fwd.vert_m << " " << d.fwd.rock_3d_m << " "
           << d.bwd.horiz_m << " " << d.bwd.vert_m << " " << d.bwd.rock_3d_m << " "
           << std::scientific << std::setprecision(6) << d.total_intensity << "\n";
           
        f2 << ang << " " << mc_counts[ang] << "\n";
    }

    // File 3: Root Events
    std::vector<int> raw_data_stream;
    for (int angle = 0; angle < 180; ++angle) {
        for (int c = 0; c < mc_counts[angle]; ++c) {
            raw_data_stream.push_back(angle);
        }
    }

    std::shuffle(raw_data_stream.begin(), raw_data_stream.end(), gen);

    std::ofstream f3(FILE3_MC_EVENTS);
    for (int val : raw_data_stream) {
        f3 << val << "\n";
    }

    std::cout << "Saved files: " << FILE1_DETAILS << ", " << FILE2_MC_SUMMARY << ", " << FILE3_MC_EVENTS << "\n";
    std::cout << "Total generated MC events (Muons + Noise): " << total_simulated << "\n";

    // ================= PHASE 3: TERMINAL ASCII HISTOGRAM =================
    std::cout << "\n================== MONTE CARLO HISTOGRAM ==================\n";
    int num_bins_ascii = (180 + ASCII_BIN_SIZE - 1) / ASCII_BIN_SIZE;
    std::vector<int> binned_counts_ascii(num_bins_ascii, 0);

    for (int i = 0; i < 180; ++i) {
        binned_counts_ascii[i / ASCII_BIN_SIZE] += mc_counts[i];
    }

    int max_bin_count_ascii = 0;
    for (int count : binned_counts_ascii) {
        if (count > max_bin_count_ascii) max_bin_count_ascii = count;
    }

    const int MAX_STARS = 50;

    for (int b = 0; b < num_bins_ascii; ++b) {
        int start_deg = b * ASCII_BIN_SIZE;
        int end_deg = std::min(start_deg + ASCII_BIN_SIZE - 1, 179);
        int count = binned_counts_ascii[b];
        int num_stars = max_bin_count_ascii > 0 ? (count * MAX_STARS / max_bin_count_ascii) : 0;
        
        std::cout << std::setw(3) << start_deg << "-" << std::setw(3) << end_deg << " deg | " 
                  << std::setw(5) << count << " | " << std::string(num_stars, '*') << "\n";
    }
    std::cout << "===========================================================\n";

    return 0;
}