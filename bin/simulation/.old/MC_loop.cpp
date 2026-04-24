#define _USE_MATH_DEFINES
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <map>
#include <tuple>
#include <chrono> 

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// ================= PARAMETERS =================
const std::string IMAGE_PATH = "D:\\MuonDensityMapping\\converted_black_white.png";
const std::string DATA_PATH = "D:\\MuonDensityMapping\\output\\all\\zussamen_10\\coin_time.txt"; 
const std::string OUTPUT_RESULTS = "D:\\MuonDensityMapping\\bin\\simulation\\simulation_output\\distance_loop_analysis.txt";
const std::string OUTPUT_FOLDER = "D:\\MuonDensityMapping\\bin\\simulation\\simulation_output\\";

const double MIN_SCAN_RADIUS_M = 10.0;  
const double MAX_SCAN_RADIUS_M = 50.0; 
const double SCAN_STEP_M = 10.0;        

const double PIXEL_SIZE = 0.13333;  
const double DETECTOR_X = 915.0;      
const double DETECTOR_Y = 1116.0;     
const double W = 0.20;  
const double H = 0.10;  
const double L = 1.0;   
const double SURFACE_Z = 15.0;       
const double STEP_SIZE_M = 0.1;
const double ROCK_DENSITY = 2.65; 

const int BACKGROUND_NOISE_PER_BIN = 0;  
const int COMPARE_BIN_SIZE = 15; 

const double A1 = -11.22, A2 = -0.00262, A3 = -14.10, A4 = -0.001213;

double crouch_mu(double h_mwe) {
    return std::exp(A1 + A2 * h_mwe) + std::exp(A3 + A4 * h_mwe);
}

std::vector<double> get_bartlett_kernel(int M) {
    std::vector<double> w(M);
    double sum = 0.0;
    for (int n = 0; n < M; ++n) {
        w[n] = (2.0 / (M - 1.0)) * ((M - 1.0) / 2.0 - std::abs(n - (M - 1.0) / 2.0));
        sum += w[n];
    }
    for (int n = 0; n < M; ++n) w[n] /= sum; 
    return w;
}

void load_data_and_time(const std::string& filename, std::vector<int>& counts, std::vector<double>& times, int& total_events) {
    counts.assign(180, 0);
    times.assign(180, 0.0);
    total_events = 0;
    
    std::ifstream file(filename);
    if (!file) {
        printf("[ERROR] Missing data file: %s\n", filename.c_str());
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

double igamc(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 0.0;
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int n = 1; n <= 100; ++n) {
            ap += 1.0; del *= x / ap; sum += del;
            if (std::abs(del) < std::abs(sum) * 3.0e-7) break;
        }
        return 1.0 - (sum * std::exp(-x + a * std::log(x) - std::lgamma(a)));
    } else {
        double b = x + 1.0 - a, c = 1.0 / 1.0e-30, d = 1.0 / b, h = d;
        for (int i = 1; i <= 100; ++i) {
            double an = -i * (i - a); b += 2.0; d = an * d + b;
            if (std::abs(d) < 1.0e-30) d = 1.0e-30;
            c = b + an / c;
            if (std::abs(c) < 1.0e-30) c = 1.0e-30;
            d = 1.0 / d;
            double del = d * c; 
            h *= del;
            if (std::abs(del - 1.0) < 3.0e-7) break;
        }
        return h * std::exp(-x + a * std::log(x) - std::lgamma(a));
    }
}

struct TestResult {
    double chi2;
    int df;
    double p_value;
    double red_chi2;
};

// ================= MAIN SIMULATION ENGINE =================
TestResult evaluate_radius(double current_radius_m, const std::vector<int>& real_data, const std::vector<double>& time_data, int total_observed, unsigned char* img_data, int width, int height, unsigned int master_seed) {
    
    double max_az = std::atan(W / L);
    double max_el = std::atan(H / L);
    
    std::vector<double> az_valid, el_valid, weights_valid;
    
    if (current_radius_m > 0.0) {
        double step_az_el = PIXEL_SIZE / current_radius_m;
        for (double el = 0.0; el <= max_el + 1e-6; el += step_az_el) {
            for (double az = -max_az; az <= max_az + 1e-6; az += step_az_el) {
                double weight = (1.0 - (L / W) * std::abs(std::tan(az))) * (1.0 - (L / H) * std::abs(std::tan(el)));
                if (weight > 0.0) {
                    az_valid.push_back(az);
                    el_valid.push_back(el);
                    weights_valid.push_back(weight);
                }
            }
        }
    } else {
        az_valid.push_back(0.0);
        el_valid.push_back(0.0);
        weights_valid.push_back(1.0);
    }

    std::vector<double> raw_flux_360(360, 0.0);
    int num_steps = static_cast<int>(std::round(current_radius_m / STEP_SIZE_M));

    int bar_width = 50;

    for (int angle_deg = 0; angle_deg < 180; ++angle_deg) {
        double side_h[2], side_3d[2], side_i[2];
        for (int side = 0; side < 2; ++side) {
            double rad = (angle_deg + side * 180) * (M_PI / 180.0);
            double eff_intensity_num = 0.0, eff_horiz_num = 0.0, eff_3d_num = 0.0, weight_sum = 0.0;

            for (size_t i = 0; i < az_valid.size(); ++i) {
                double ray_az = rad + az_valid[i];
                double dir_x = std::sin(ray_az), dir_y = -std::cos(ray_az), tan_el = std::tan(el_valid[i]);
                double horiz_rock_m = 0.0;

                if (current_radius_m > 0.0) {
                    for (int s = 0; s <= num_steps; ++s) {
                        double r = s * STEP_SIZE_M;
                        int xi = static_cast<int>(std::round(DETECTOR_X + (r / PIXEL_SIZE) * dir_x));
                        int yi = static_cast<int>(std::round(DETECTOR_Y + (r / PIXEL_SIZE) * dir_y));

                        if (xi >= 0 && xi < width && yi >= 0 && yi < height) {
                            if (img_data[yi * width + xi] < 128 && r * tan_el <= SURFACE_Z) {
                                horiz_rock_m += STEP_SIZE_M;
                            }
                        }
                    }
                }

                double h_3d = horiz_rock_m / std::cos(el_valid[i]);
                double h_mwe = h_3d * ROCK_DENSITY;
                double zenith_rad = (M_PI / 2.0) - el_valid[i];
                double intensity = std::pow(std::cos(zenith_rad), 2) * crouch_mu(h_mwe);
                
                eff_intensity_num += intensity * weights_valid[i];
                eff_horiz_num += horiz_rock_m * weights_valid[i];
                eff_3d_num += h_3d * weights_valid[i];
                weight_sum += weights_valid[i];
            }
            side_h[side] = eff_horiz_num / weight_sum;
            side_3d[side] = eff_3d_num / weight_sum;
            side_i[side] = eff_intensity_num / weight_sum;
            raw_flux_360[angle_deg + side * 180] = side_i[side];
        }

        // --- DYNAMIC PROGRESS BAR ---
        float progress = (float)(angle_deg + 1) / 180.0f;
        int pos = bar_width * progress;
        printf("\r[");
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) printf("=");
            else if (i == pos) printf(">");
            else printf(" ");
        }
        printf("] %3d%% (Angle: %3d/180) | R = %.2fm", int(progress * 100.0), angle_deg + 1, current_radius_m);
        fflush(stdout);
    }
    
    printf("\n"); // Nowa linia po zakończeniu paska ładowania
    fflush(stdout);

    // Smoothing
    std::vector<double> kernel = get_bartlett_kernel(23);
    std::vector<double> expected_1deg(180, 0.0);
    
    for (int i = 0; i < 180; ++i) {
        double smoothed = 0.0;
        for (int j = 0; j < (int)kernel.size(); ++j) {
            smoothed += raw_flux_360[(i + j - 11 + 360) % 360] * kernel[j];
            smoothed += raw_flux_360[(i + 180 + j - 11 + 360) % 360] * kernel[j];
        }
        expected_1deg[i] = smoothed;
    }

    double sum_time_flux = 0.0;
    for (int i = 0; i < 180; ++i) {
        sum_time_flux += expected_1deg[i] * time_data[i];
    }

    for (int i = 0; i < 180; ++i) {
        expected_1deg[i] = (expected_1deg[i] * time_data[i] / sum_time_flux) * total_observed;
        expected_1deg[i] += BACKGROUND_NOISE_PER_BIN;
    }

    std::mt19937 gen(master_seed); 
    std::vector<int> raw_data_stream;
    
    int num_compare_bins = (180 + COMPARE_BIN_SIZE - 1) / COMPARE_BIN_SIZE; 
    std::vector<double> mc_binned(num_compare_bins, 0.0);

    for (int i = 0; i < 180; ++i) {
        double mean = std::max(1e-6, expected_1deg[i]);
        std::poisson_distribution<> d(mean);
        int count = d(gen);
        
        mc_binned[i / COMPARE_BIN_SIZE] += count;
        for (int c = 0; c < count; ++c) raw_data_stream.push_back(i);
    }
    std::shuffle(raw_data_stream.begin(), raw_data_stream.end(), gen);

    char filename[512];
    snprintf(filename, sizeof(filename), "%sMC_one_by_one\\%.2f_metrow.txt", OUTPUT_FOLDER.c_str(), current_radius_m);
    
    FILE* fout = fopen(filename, "w");
    if (fout) {
        for (int val : raw_data_stream) fprintf(fout, "%d\n", val);
        fclose(fout);
    }

    double chi2 = 0.0;
    int df = 0;

    // === TUTAJ PRZYWRÓCONA TABELKA ===
    printf("---> Chi-Square analysis in %d degree bins <---\n", COMPARE_BIN_SIZE);
    printf("%-10s | %-12s | %-12s | %-10s\n", "Bin No.", "Observed(O)", "Expected(E)", "Contrib");
    
    for (int b = 0; b < num_compare_bins; ++b) {
        double O_bin = 0.0, E_bin = 0.0;
        for (int i = 0; i < COMPARE_BIN_SIZE; ++i) {
            if (b * COMPARE_BIN_SIZE + i < 180) { 
                O_bin += real_data[b * COMPARE_BIN_SIZE + i];
            }
        }
        E_bin = mc_binned[b]; 
        
        double contrib = 0.0;
        if (E_bin > 0) {
            contrib = ((O_bin - E_bin) * (O_bin - E_bin)) / E_bin;
            chi2 += contrib;
            df++;
        }
        // UROCZYSTY POWRÓT WYDRUKU:
        printf("Bin %2d    | %-12.0f | %-12.2f | %-10.4f\n", b, O_bin, E_bin, contrib);
    }
    
    df = std::max(1, df - 1); 
    double red_chi2 = chi2 / df;
    double p_value = igamc(df / 2.0, chi2 / 2.0);

    TestResult res;
    res.chi2 = chi2; res.df = df; res.p_value = p_value; res.red_chi2 = red_chi2;
    return res;
}

int main() {
    printf("1. Loading cave map and detector data...\n");
    int width, height, channels, total_observed = 0;
    unsigned char *img_data = stbi_load(IMAGE_PATH.c_str(), &width, &height, &channels, 1);
    if (!img_data) {
        printf("[ERROR] Cannot load file %s\n", IMAGE_PATH.c_str());
        return 1;
    }

    std::vector<int> real_data;
    std::vector<double> time_data;
    load_data_and_time(DATA_PATH, real_data, time_data, total_observed);

    unsigned int master_seed = std::chrono::system_clock::now().time_since_epoch().count();
    
    printf("\n==========================================================\n");
    printf("[MONTE CARLO] Generated MASTER SEED for this run: %u\n", master_seed);
    printf("==========================================================\n\n");

    std::map<double, TestResult> cache_results;
    
    auto get_eval = [&](double r) {
        if (cache_results.find(r) != cache_results.end()) return cache_results[r];
        TestResult res = evaluate_radius(r, real_data, time_data, total_observed, img_data, width, height, master_seed);
        cache_results[r] = res;
        return res;
    };

    // ================= PHASE 1: COARSE SCAN =================
    printf("\n2. Starting coarse scan from %.1fm to %.1fm...\n", MIN_SCAN_RADIUS_M, MAX_SCAN_RADIUS_M);
    double best_r_coarse = MIN_SCAN_RADIUS_M;
    double min_red_chi2 = 1e12;

    for (double r = MIN_SCAN_RADIUS_M; r <= MAX_SCAN_RADIUS_M + 1e-6; r += SCAN_STEP_M) {
        TestResult res = get_eval(r);
        printf("=> SCAN RESULT %.1fm | Chi^2: %.2f | DF: %d | Red. Chi^2: %.4f | P-Value: %e <=\n", 
               r, res.chi2, res.df, res.red_chi2, res.p_value);
        
        double dist = std::abs(res.red_chi2 - 1.0);
        if (dist < min_red_chi2) {
            min_red_chi2 = dist;
            best_r_coarse = r;
        }
    }

    double L = std::max(MIN_SCAN_RADIUS_M, best_r_coarse - SCAN_STEP_M);
    double R = std::min(MAX_SCAN_RADIUS_M, best_r_coarse + SCAN_STEP_M);

    printf("\nBest coarse interval is: [%.1f, %.1f]\n", L, R);

    // ================= PHASE 2: 5 BISECTION STEPS =================
    printf("\n3. Densifying interval (Bisection - 5 iterations)...\n");
    for (int step = 0; step < 5; ++step) {
        double mid = (L + R) / 2.0;
        printf("\nStep %d: Dividing in half. Checking R = %.4fm...\n", step + 1, mid);
        
        TestResult res_mid = get_eval(mid);
        double dist_L = std::abs(get_eval(L).red_chi2 - 1.0);
        double dist_R = std::abs(get_eval(R).red_chi2 - 1.0);

        if (dist_L < dist_R) {
            R = mid;
        } else {
            L = mid;
        }
    }

    // ================= SAVING AND SUMMARY =================
    double final_result = std::round((L + R) / 2.0);
    get_eval(final_result); 
    stbi_image_free(img_data);

    FILE* out = fopen(OUTPUT_RESULTS.c_str(), "w");
    if (out) {
        fprintf(out, "Radius_m Chi2 DF Reduced_Chi2 P_Value\n");
        for(auto const& pair : cache_results) {
            fprintf(out, "%.4f %.2f %d %.4f %e\n", 
                    pair.first, pair.second.chi2, pair.second.df, pair.second.red_chi2, pair.second.p_value);
        }
        fclose(out);
    }

    std::vector<std::pair<double, TestResult>> sorted_results(cache_results.begin(), cache_results.end());
    std::sort(sorted_results.begin(), sorted_results.end(), [](const auto& a, const auto& b) {
        return std::abs(a.second.red_chi2 - 1.0) < std::abs(b.second.red_chi2 - 1.0);
    });

    printf("\n============================================================\n");
    printf("FINAL ANSWER: OPTIMAL RADII (Closest to 1.0)\n");
    printf("============================================================\n");
    for(int i = 0; i < 3 && i < sorted_results.size(); ++i) {
        printf("%d. Radius: %5.2f m | Red. Chi2 = %.4f | P-value = %.4f\n", 
               i+1, sorted_results[i].first, sorted_results[i].second.red_chi2, sorted_results[i].second.p_value);
    }
    printf("------------------------------------------------------------\n");
    printf("Full scan history saved to file: %s\n", OUTPUT_RESULTS.c_str());

    return 0;
}