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
const std::string IMAGE_PATH = "/home/kacper/MuonDensityMapping/converted_black_white.png";
//const std::string IMAGE_PATH = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_data/converted_black_white_4.png";
const std::string DATA_PATH = "/home/kacper/MuonDensityMapping/output/all/zussamen_10/coin_time.txt"; 
const std::string OUTPUT_RESULTS = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/distance_analysis/distance_one_by_one_analysis.txt";
const std::string OUTPUT_FOLDER = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/distance_analysis_exp/";
const std::string OUTPUT_BEST_BIN_ANALYSIS = "/home/kacper/MuonDensityMapping/bin/simulation/simulation_output/distance_analysis/best_radius_bin_analysis.txt";

const double MIN_SCAN_RADIUS_M = 0.0;  
const double MAX_SCAN_RADIUS_M = 100.0; 
const double SCAN_STEP_M = 1.0;        

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
const int COMPARE_BIN_SIZE = 5; 

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
TestResult evaluate_radius(double current_radius_m, const std::vector<int>& real_data, const std::vector<double>& time_data, int total_observed, unsigned char* img_data, int width, int height) {
    printf("\n============================================================\n");
    printf("--- Physics scanning and Chi2 test for Exp Radius R0 = %.2f m ---\n", current_radius_m);
    printf("============================================================\n");
    fflush(stdout);
    
    double max_az = std::atan(W / L);
    double max_el = std::atan(H / L);
    
    std::vector<double> az_valid, el_valid, weights_valid;
    
    if (current_radius_m > 0.0) {
        // Poprawka: Rozdzielczość kątowa wyliczana zawsze z MAX_SCAN_RADIUS_M
        double step_az_el = PIXEL_SIZE / MAX_SCAN_RADIUS_M;
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
    
    // Poprawka: Idziemy zawsze do końca mapy (obszaru analizy)
    int max_map_steps = static_cast<int>(std::round(MAX_SCAN_RADIUS_M / STEP_SIZE_M));

    printf("Generating ray matrix (Theoretical Flux)...\n");
    fflush(stdout);
    
    int bar_width = 50; 

    for (int angle_deg = 0; angle_deg < 180; ++angle_deg) {
        for (int side = 0; side < 2; ++side) {
            double rad = (angle_deg + side * 180) * (M_PI / 180.0);
            double eff_intensity_num = 0.0, weight_sum = 0.0;

            for (size_t i = 0; i < az_valid.size(); ++i) {
                double ray_az = rad + az_valid[i];
                double dir_x = std::sin(ray_az), dir_y = -std::cos(ray_az), tan_el = std::tan(el_valid[i]);
                double horiz_rock_m = 0.0;

                if (current_radius_m > 0.0) {
                    for (int s = 0; s <= max_map_steps; ++s) {
                        double r = s * STEP_SIZE_M;
                        int xi = static_cast<int>(std::round(DETECTOR_X + (r / PIXEL_SIZE) * dir_x));
                        int yi = static_cast<int>(std::round(DETECTOR_Y + (r / PIXEL_SIZE) * dir_y));

                        if (xi >= 0 && xi < width && yi >= 0 && yi < height) {
                            bool is_rock = (img_data[yi * width + xi] < 128 && r * tan_el <= SURFACE_Z);
                            if (is_rock) {
                                horiz_rock_m += STEP_SIZE_M;
                            } else {
                                // Model płynnego zanikania wpływu pustych przestrzeni
                                double visibility_weight = std::exp(-r / current_radius_m); 
                                horiz_rock_m += STEP_SIZE_M * (1.0 - visibility_weight);
                            }
                        } else {
                            // Poza mapą = lita skała
                            horiz_rock_m += STEP_SIZE_M; 
                        }
                    }
                }

                double h_mwe = (horiz_rock_m / std::cos(el_valid[i])) * ROCK_DENSITY;
                double zenith_rad = (M_PI / 2.0) - el_valid[i];
                double intensity = std::pow(std::cos(zenith_rad), 2) * crouch_mu(h_mwe);
                
                eff_intensity_num += intensity * weights_valid[i];
                weight_sum += weights_valid[i];
            }
            raw_flux_360[angle_deg + side * 180] = eff_intensity_num / weight_sum;
        }

        float progress = (float)(angle_deg + 1) / 180.0f;
        int pos = bar_width * progress;
        printf("\r[");
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) printf("=");
            else if (i == pos) printf(">");
            else printf(" ");
        }
        printf("] %3d%% (Angle: %3d/180) | R0 = %.2fm", int(progress * 100.0), angle_deg + 1, current_radius_m);
        fflush(stdout);
    }
    
    printf("\nRay-tracing finished.\n");

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

    // Normalization to observed data
    double sum_time_flux = 0.0;
    for (int i = 0; i < 180; ++i) {
        sum_time_flux += expected_1deg[i] * time_data[i];
    }

    for (int i = 0; i < 180; ++i) {
        expected_1deg[i] = (expected_1deg[i] * time_data[i] / sum_time_flux) * total_observed;
        expected_1deg[i] += BACKGROUND_NOISE_PER_BIN;
    }

    // Save theoretical binned data
    char filename[512];
    snprintf(filename, sizeof(filename), "%s%.2f_metrow_exp.txt", OUTPUT_FOLDER.c_str(), current_radius_m);
    
    FILE* fout = fopen(filename, "w");
    if (fout) {
        fprintf(fout, "Angle_deg Theoretical_Counts\n");
        for (int i = 0; i < 180; ++i) fprintf(fout, "%d %.6f\n", i + 1, expected_1deg[i]);
        fclose(fout);
    }

    // --- Statistical Comparison ---
    double chi2 = 0.0;
    int df = 0;
    int num_compare_bins = (180 + COMPARE_BIN_SIZE - 1) / COMPARE_BIN_SIZE; 

    printf("---> Chi-Square analysis (Observed vs Theory) in %d degree bins <---\n", COMPARE_BIN_SIZE);
    printf("%-10s | %-12s | %-12s | %-10s\n", "Bin No.", "Observed(O)", "Theory(E)", "Contrib");
    
    for (int b = 0; b < num_compare_bins; ++b) {
        double O_bin = 0.0, E_bin = 0.0;
        for (int i = 0; i < COMPARE_BIN_SIZE; ++i) {
            int idx = b * COMPARE_BIN_SIZE + i;
            if (idx < 180) { 
                O_bin += real_data[idx];
                E_bin += expected_1deg[idx];
            }
        }
        
        double contrib = 0.0;
        if (E_bin > 0) {
            contrib = ((O_bin - E_bin) * (O_bin - E_bin)) / E_bin;
            chi2 += contrib;
            df++;
        }
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

    std::map<double, TestResult> cache_results;
    
    printf("\n2. Starting direct theoretical scan for R0 from %.1fm to %.1fm with step %.1fm...\n", MIN_SCAN_RADIUS_M, MAX_SCAN_RADIUS_M, SCAN_STEP_M);

    for (double r = MIN_SCAN_RADIUS_M; r <= MAX_SCAN_RADIUS_M + 1e-6; r += SCAN_STEP_M) {
        TestResult res = evaluate_radius(r, real_data, time_data, total_observed, img_data, width, height);
        cache_results[r] = res;
        
        printf("=> SCAN RESULT %.1fm | Chi^2: %.2f | DF: %d | Red. Chi^2: %.4f | P-Value: %e <=\n", 
               r, res.chi2, res.df, res.red_chi2, res.p_value);
    }

    // ZAPIS GŁÓWNYCH WYNIKÓW
    FILE* out = fopen(OUTPUT_RESULTS.c_str(), "w");
    if (out) {
        fprintf(out, "RadiusR0_m Chi2 DF Reduced_Chi2 P_Value\n");
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
    printf("FINAL ANSWER: OPTIMAL EXPONENTIAL RADII (R0)\n");
    printf("============================================================\n");
    for(int i = 0; i < 3 && i < (int)sorted_results.size(); ++i) {
        printf("%d. Radius R0: %5.2f m | Red. Chi2 = %.4f | P-value = %.4f\n", 
               i+1, sorted_results[i].first, sorted_results[i].second.red_chi2, sorted_results[i].second.p_value);
    }

    // ================= DETAILED ANALYSIS FOR BEST RADIUS =================
    double best_radius = sorted_results[0].first;
    TestResult best_res = sorted_results[0].second;

    printf("\n============================================================\n");
    printf("DETAILED BIN-BY-BIN ANALYSIS FOR BEST EXP RADIUS: %.2f m\n", best_radius);
    printf("============================================================\n");

    // Zaktualizowany blok z nową logiką rozdzielczości i płynnego zanikania
    double max_az2 = std::atan(W / L);
    double max_el2 = std::atan(H / L);
    std::vector<double> az_valid2, el_valid2, weights_valid2;
    
    double step_az_el = PIXEL_SIZE / MAX_SCAN_RADIUS_M; // Poprawka rozdzielczości
    
    for (double el = 0.0; el <= max_el2 + 1e-6; el += step_az_el) {
        for (double az = -max_az2; az <= max_az2 + 1e-6; az += step_az_el) {
            double weight = (1.0 - (L / W) * std::abs(std::tan(az))) * (1.0 - (L / H) * std::abs(std::tan(el)));
            if (weight > 0.0) {
                az_valid2.push_back(az); el_valid2.push_back(el); weights_valid2.push_back(weight);
            }
        }
    }

    std::vector<double> raw_flux_360b(360, 0.0);
    
    // Zawsze skanujemy do samego końca mapy
    int max_map_steps2 = static_cast<int>(std::round(MAX_SCAN_RADIUS_M / STEP_SIZE_M));

    printf("Recomputing ray-tracing for R0 = %.2f m...\n", best_radius);
    for (int angle_deg = 0; angle_deg < 180; ++angle_deg) {
        for (int side = 0; side < 2; ++side) {
            double rad = (angle_deg + side * 180) * (M_PI / 180.0);
            double eff_intensity_num = 0.0, weight_sum = 0.0;
            
            for (size_t i = 0; i < az_valid2.size(); ++i) {
                double ray_az = rad + az_valid2[i];
                double dir_x = std::sin(ray_az), dir_y = -std::cos(ray_az), tan_el = std::tan(el_valid2[i]);
                double horiz_rock_m = 0.0;
                
                if (best_radius > 0.0) {
                    for (int s = 0; s <= max_map_steps2; ++s) {
                        double r = s * STEP_SIZE_M;
                        int xi = static_cast<int>(std::round(DETECTOR_X + (r / PIXEL_SIZE) * dir_x));
                        int yi = static_cast<int>(std::round(DETECTOR_Y + (r / PIXEL_SIZE) * dir_y));
                        
                        if (xi >= 0 && xi < width && yi >= 0 && yi < height) {
                            bool is_rock = (img_data[yi * width + xi] < 128 && r * tan_el <= SURFACE_Z);
                            if (is_rock) {
                                horiz_rock_m += STEP_SIZE_M;
                            } else {
                                // Tłumienie eksponencjalne oparte na najlepszym dopasowanym parametrze
                                double visibility_weight = std::exp(-r / best_radius); 
                                horiz_rock_m += STEP_SIZE_M * (1.0 - visibility_weight);
                            }
                        } else {
                            horiz_rock_m += STEP_SIZE_M;
                        }
                    }
                }
                
                double h_mwe = (horiz_rock_m / std::cos(el_valid2[i])) * ROCK_DENSITY;
                eff_intensity_num += std::pow(std::cos((M_PI/2.0)-el_valid2[i]), 2) * crouch_mu(h_mwe) * weights_valid2[i];
                weight_sum += weights_valid2[i];
            }
            raw_flux_360b[angle_deg + side * 180] = eff_intensity_num / weight_sum;
        }
    }

    // Smoothing & Normalization
    std::vector<double> kernel2 = get_bartlett_kernel(23);
    std::vector<double> expected_1deg2(180, 0.0);
    for (int i = 0; i < 180; ++i) {
        double smoothed = 0.0;
        for (int j = 0; j < (int)kernel2.size(); ++j) {
            smoothed += raw_flux_360b[(i + j - 11 + 360) % 360] * kernel2[j];
            smoothed += raw_flux_360b[(i + 180 + j - 11 + 360) % 360] * kernel2[j];
        }
        expected_1deg2[i] = smoothed;
    }
    double sum_f = 0; 
    for(int i=0; i<180; i++) sum_f += expected_1deg2[i] * time_data[i];
    for(int i=0; i<180; i++) expected_1deg2[i] = (expected_1deg2[i] * time_data[i] / sum_f) * total_observed;

    // Wydruk i zapis binów
    FILE* fdetail = fopen(OUTPUT_BEST_BIN_ANALYSIS.c_str(), "w");
    printf("\n%-8s | %-12s | %-12s | %-10s\n", "Bin No.", "Observed(O)", "Theory(E)", "Contrib");
    if (fdetail) fprintf(fdetail, "Bin AngleFrom AngleTo Obs Theory Contrib\n");

    double chi2_final = 0;
    for (int b = 0; b < (180/COMPARE_BIN_SIZE); ++b) {
        double O_bin = 0, E_bin = 0;
        for (int i = 0; i < COMPARE_BIN_SIZE; ++i) {
            O_bin += real_data[b * COMPARE_BIN_SIZE + i];
            E_bin += expected_1deg2[b * COMPARE_BIN_SIZE + i];
        }
        double con = (E_bin > 0) ? ((O_bin-E_bin)*(O_bin-E_bin)/E_bin) : 0;
        chi2_final += con;
        printf("Bin %2d    | %-12.0f | %-12.2f | %-10.4f\n", b, O_bin, E_bin, con);
        if (fdetail) fprintf(fdetail, "%d %d %d %.0f %.2f %.4f\n", b, b*9+1, (b+1)*9, O_bin, E_bin, con);
    }
    printf("------------------------------------------------------------\n");
    printf("FINAL Chi2: %.4f\n", chi2_final);

    if (fdetail) fclose(fdetail);
    stbi_image_free(img_data); 

    return 0;
}