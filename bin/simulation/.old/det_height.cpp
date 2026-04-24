#define _USE_MATH_DEFINES
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <tuple>
#include <chrono>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// ================= PARAMETERS =================
const std::string IMAGE_PATH = "bin/simulation/simulation_data/converted_black_white.png";
const std::string DATA_PATH = "output/all/zussamen_10/coin_time.txt"; 
const std::string OUTPUT_RESULTS = "bin/simulation/simulation_output/height_analysis.txt";
const std::string OUTPUT_BEST_HEIGHT_ANALYSIS = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/height_analysis.txt";

const double MIN_HEIGHT_M = 0.0;  
const double MAX_HEIGHT_M = 15.0; 
const double HEIGHT_STEP_M = 0.1;        

const double SCAN_RADIUS_M = 100.0; 

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
const int COMPARE_BIN_SIZE = 9; 

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
    if (!file) { printf("[ERROR] Missing data file: %s\n", filename.c_str()); exit(1); }
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
        if (idx >= 0 && idx < 180) { counts[idx] = count; times[idx] = time; total_events += count; }
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
            double del = d * c; h *= del;
            if (std::abs(del - 1.0) < 3.0e-7) break;
        }
        return h * std::exp(-x + a * std::log(x) - std::lgamma(a));
    }
}

struct TestResult {
    double chi2;
    int df;
    double p_value;
    std::vector<double> expected_1deg;
};

// ================= MAIN SIMULATION ENGINE =================
TestResult evaluate_height(double current_height_m, const std::vector<int>& real_data, const std::vector<double>& time_data, int total_observed, unsigned char* img_data, int width, int height) {
    printf("\n============================================================\n");
    printf("--- Physics scanning: Height H = %.2f m ---\n", current_height_m);
    printf("============================================================\n");
    fflush(stdout); 
    
    double max_az = std::atan(W / L);
    double max_el = std::atan(H / L);
    std::vector<double> az_valid, el_valid, weights_valid;
    double step_az_el = PIXEL_SIZE / SCAN_RADIUS_M; 
    
    for (double el = 0.0; el <= max_el + 1e-6; el += step_az_el) {
        for (double az = -max_az; az <= max_az + 1e-6; az += step_az_el) {
            double weight = (1.0 - (L / W) * std::abs(std::tan(az))) * (1.0 - (L / H) * std::abs(std::tan(el)));
            if (weight > 0.0) {
                az_valid.push_back(az); el_valid.push_back(el); weights_valid.push_back(weight);
            }
        }
    }

    std::vector<double> raw_flux_360(360, 0.0);
    double MAX_RAY_DIST = 500.0; 
    int bar_width = 50;
    
    for (int angle_deg = 0; angle_deg < 180; ++angle_deg) {
        for (int side = 0; side < 2; ++side) {
            double rad = (angle_deg + side * 180) * (M_PI / 180.0);
            double eff_intensity_num = 0.0, weight_sum = 0.0;

            for (size_t i = 0; i < az_valid.size(); ++i) {
                double ray_az = rad + az_valid[i];
                double dir_x = std::sin(ray_az), dir_y = -std::cos(ray_az), tan_el = std::tan(el_valid[i]);
                double horiz_rock_m = 0.0;

                // KLUCZOWA ZMIANA: pętla r idzie do SURFACE_Z, nie kończy na SCAN_RADIUS_M
                for (double r = 0.0; r <= MAX_RAY_DIST; r += STEP_SIZE_M) {
                    double z = r * tan_el;
                    if (z > SURFACE_Z) break; 

                    bool is_rock = true;
                    // Jaskinię sprawdzamy tylko tam, gdzie mamy mapę
                    if (r <= SCAN_RADIUS_M) {
                        int xi = static_cast<int>(std::round(DETECTOR_X + (r / PIXEL_SIZE) * dir_x));
                        int yi = static_cast<int>(std::round(DETECTOR_Y + (r / PIXEL_SIZE) * dir_y));
                        if (xi >= 0 && xi < width && yi >= 0 && yi < height) {
                            if (img_data[yi * width + xi] >= 128 && z <= current_height_m) is_rock = false;
                        }
                    } 
                    if (is_rock) horiz_rock_m += STEP_SIZE_M;
                }

                double h_mwe = (horiz_rock_m / std::cos(el_valid[i])) * ROCK_DENSITY;
                double intensity = std::pow(std::cos((M_PI / 2.0) - el_valid[i]), 2) * crouch_mu(h_mwe);
                eff_intensity_num += intensity * weights_valid[i];
                weight_sum += weights_valid[i];
            }
            raw_flux_360[angle_deg + side * 180] = eff_intensity_num / weight_sum;
        }

        // TWOJA ORYGINALNA PASEK POSTEPU
        float progress = (float)(angle_deg + 1) / 180.0f;
        int pos = bar_width * progress;
        printf("\r[");
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) printf("="); else if (i == pos) printf(">"); else printf(" ");
        }
        printf("] %3d%% (Angle: %3d/180) | Height = %.2fm", int(progress * 100.0), angle_deg + 1, current_height_m);
        fflush(stdout); 
    }
    printf("\nRay-tracing finished.\n");

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
    for (int i = 0; i < 180; ++i) sum_time_flux += expected_1deg[i] * time_data[i];
    for (int i = 0; i < 180; ++i) {
        expected_1deg[i] = (expected_1deg[i] * time_data[i] / sum_time_flux) * total_observed + BACKGROUND_NOISE_PER_BIN;
    }

    double chi2 = 0.0; int df_bins = 0;
    int num_compare_bins = (180 + COMPARE_BIN_SIZE - 1) / COMPARE_BIN_SIZE; 
    for (int b = 0; b < num_compare_bins; ++b) {
        double O_bin = 0.0, E_bin = 0.0;
        for (int i = 0; i < COMPARE_BIN_SIZE; ++i) {
            int idx = b * COMPARE_BIN_SIZE + i;
            if (idx < 180) { O_bin += real_data[idx]; E_bin += expected_1deg[idx]; }
        }
        if (E_bin > 0) { chi2 += ((O_bin - E_bin) * (O_bin - E_bin)) / E_bin; df_bins++; }
    }
    
    TestResult res; res.chi2 = chi2; res.df = std::max(1, df_bins - 1); 
    res.p_value = igamc(res.df / 2.0, chi2 / 2.0); res.expected_1deg = expected_1deg;
    return res;
}

int main() {
    printf("1. Loading cave map and detector data...\n");
    int width, height, channels, total_observed = 0;
    unsigned char *img_data = stbi_load(IMAGE_PATH.c_str(), &width, &height, &channels, 1);
    if (!img_data) return 1;

    std::vector<int> real_data; std::vector<double> time_data;
    load_data_and_time(DATA_PATH, real_data, time_data, total_observed);

    std::map<double, TestResult> cache_results;
    for (double h = MIN_HEIGHT_M; h <= MAX_HEIGHT_M + 1e-6; h += HEIGHT_STEP_M) {
        cache_results[h] = evaluate_height(h, real_data, time_data, total_observed, img_data, width, height);
        printf("=> RESULT H=%.1fm | Chi^2: %.2f | P-Value: %e <=\n", h, cache_results[h].chi2, cache_results[h].p_value);
    }
    stbi_image_free(img_data);

    std::vector<std::pair<double, TestResult>> sorted(cache_results.begin(), cache_results.end());
    std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) { return a.second.chi2 < b.second.chi2; });

    FILE* out = fopen(OUTPUT_RESULTS.c_str(), "w");
    if (out) {
        fprintf(out, "Height_m Chi2 DF P_Value\n");
        for(auto const& [h, res] : cache_results) fprintf(out, "%.4f %.2f %d %e\n", h, res.chi2, res.df, res.p_value);
        fclose(out);
    }
    
    printf("\nOPTIMAL HEIGHT: %.2f m (Chi2: %.2f)\n", sorted[0].first, sorted[0].second.chi2);
    return 0;
}