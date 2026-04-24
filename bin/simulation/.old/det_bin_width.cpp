#define _USE_MATH_DEFINES
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// ================= PARAMETERS =================
const std::string IMAGE_PATH = "bin/simulation/simulation_data/converted_black_white.png";
const std::string DATA_PATH = "output/all/zussamen_10/coin_time.txt"; 
const std::string OUTPUT_RESULTS = "bin/simulation/simulation_output/bin_analysis.txt";
const std::string OUTPUT_BEST_ANGLE_ANALYSIS = "bin/simulation/simulation_output/best_binsize_angle_analysis.txt";

const double TEST_RADIUS_M = 30.0;  

const double PIXEL_SIZE = 0.13333;  
const double DETECTOR_X = 915.0;      
const double DETECTOR_Y = 1116.0;     
const double W = 0.20;  
const double H = 0.10;  
const double L = 1.0;   
const double SURFACE_Z = 15.0;       
const double STEP_SIZE_M = 0.1;
const double ROCK_DENSITY = 2.65; 

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
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * 3.0e-7) break;
        }
        return 1.0 - (sum * std::exp(-x + a * std::log(x) - std::lgamma(a)));
    } else {
        double b = x + 1.0 - a, c = 1.0 / 1.0e-30, d = 1.0 / b, h = d;
        for (int i = 1; i <= 100; ++i) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
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

int main() {
    printf("============================================================\n");
    printf("PHASE 1: THEORETICAL CALCULATIONS (R = %.1fm)\n", TEST_RADIUS_M);
    printf("============================================================\n");
    
    int width, height, channels;
    unsigned char *img_data = stbi_load(IMAGE_PATH.c_str(), &width, &height, &channels, 1);
    if (!img_data) {
        printf("[ERROR] Cannot load image: %s\n", IMAGE_PATH.c_str());
        return 1;
    }

    int total_observed = 0;
    std::vector<int> real_counts_1deg;
    std::vector<double> time_data_1deg;
    load_data_and_time(DATA_PATH, real_counts_1deg, time_data_1deg, total_observed);
    printf("Loaded %d events from the detector.\n\n", total_observed);

    double max_az = std::atan(W / L);
    double max_el = std::atan(H / L);
    std::vector<double> az_valid, el_valid, weights_valid;
    
    if (TEST_RADIUS_M > 0.0) {
        double step_az_el = PIXEL_SIZE / TEST_RADIUS_M;
        for (double el = 0.0; el <= max_el + 1e-6; el += step_az_el) {
            for (double az = -max_az; az <= max_az + 1e-6; az += step_az_el) {
                double w_h = 1.0 - (L / W) * std::abs(std::tan(az));
                double w_v = 1.0 - (L / H) * std::abs(std::tan(el));
                double weight = w_h * w_v;
                if (weight > 0.0) {
                    az_valid.push_back(az); el_valid.push_back(el); weights_valid.push_back(weight);
                }
            }
        }
    } else {
        az_valid.push_back(0.0); el_valid.push_back(0.0); weights_valid.push_back(1.0);
    }

    std::vector<double> raw_flux_360(360, 0.0);
    int num_steps = static_cast<int>(std::round(TEST_RADIUS_M / STEP_SIZE_M));
    int bar_width = 50;

    printf("Generating ray matrix (Deterministic Model)...\n");
    fflush(stdout);

    for (int angle_deg = 0; angle_deg < 180; ++angle_deg) {
        for (int side = 0; side < 2; ++side) {
            double rad = (angle_deg + side * 180) * (M_PI / 180.0);
            double eff_i = 0, weight_sum = 0;
            
            for (size_t i = 0; i < az_valid.size(); ++i) {
                double ray_az = rad + az_valid[i];
                double dir_x = std::sin(ray_az), dir_y = -std::cos(ray_az), tan_el = std::tan(el_valid[i]);
                double h_horiz = 0.0;

                if (TEST_RADIUS_M > 0.0) {
                    for (int s = 0; s <= num_steps; ++s) {
                        double r = s * STEP_SIZE_M;
                        int xi = (int)std::round(DETECTOR_X + (r / PIXEL_SIZE) * dir_x);
                        int yi = (int)std::round(DETECTOR_Y + (r / PIXEL_SIZE) * dir_y);
                        if (xi >= 0 && xi < width && yi >= 0 && yi < height) {
                            if (img_data[yi * width + xi] < 128 && r * tan_el <= SURFACE_Z) h_horiz += STEP_SIZE_M;
                        }
                    }
                }

                double h_3d = h_horiz / std::cos(el_valid[i]);
                double h_mwe = h_3d * ROCK_DENSITY;
                double zenith_rad = (M_PI / 2.0) - el_valid[i];
                double intensity = std::pow(std::cos(zenith_rad), 2) * crouch_mu(h_mwe);

                eff_i += intensity * weights_valid[i];
                weight_sum += weights_valid[i];
            }
            raw_flux_360[angle_deg + side * 180] = eff_i / weight_sum;
        }

        float progress = (float)(angle_deg + 1) / 180.0f;
        int pos = bar_width * progress;
        printf("\r[");
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) printf("="); else if (i == pos) printf(">"); else printf(" ");
        }
        printf("] %3d%% (Angle: %3d/180)", int(progress * 100.0), angle_deg + 1);
        fflush(stdout);
    }
    
    printf("\nRay-tracing finished.\n");
    stbi_image_free(img_data);

    // --- SMOOTHING ---
    std::vector<double> kernel = get_bartlett_kernel(23);
    std::vector<double> expected_1deg(180, 0.0);
    
    for (int i = 0; i < 180; ++i) {
        double smoothed = 0;
        for (int j = 0; j < 23; ++j) {
            smoothed += raw_flux_360[(i + j - 11 + 360) % 360] * kernel[j];
            smoothed += raw_flux_360[(i + 180 + j - 11 + 360) % 360] * kernel[j];
        }
        expected_1deg[i] = smoothed;
    }

    // --- NORMALIZATION TO EXPOSURE TIME AND TOTAL EVENTS ---
    double sum_time_flux = 0.0;
    for (int i = 0; i < 180; ++i) {
        sum_time_flux += expected_1deg[i] * time_data_1deg[i];
    }
    for (int i = 0; i < 180; ++i) {
        expected_1deg[i] = (expected_1deg[i] * time_data_1deg[i] / sum_time_flux) * total_observed;
    }

    // ================= TESTING BIN SIZES (DIRECT COMPARISON) =================
    printf("\n============================================================\n");
    printf("PHASE 2: CHI-SQUARE STATISTICAL ANALYSIS (DETERMINISTIC)\n");
    printf("============================================================\n");

    std::vector<int> bin_sizes = {1, 2, 3, 4, 5, 6, 9, 10, 12, 15, 18, 20, 30, 36, 45, 60, 90, 180};
    
    struct Summary { int size; double chi2; int df; double red_chi2; double p_value; };
    std::vector<Summary> final_report;

    FILE* out = fopen(OUTPUT_RESULTS.c_str(), "w");
    if (out) fprintf(out, "BinSize_deg Chi2 DF Reduced_Chi2 P_Value\n");

    for (int bw : bin_sizes) {
        printf("\n>>> CHECKING BIN WIDTH: %d degrees <<<\n", bw);
        printf("%-10s | %-12s | %-12s | %-10s\n", "Bin No.", "Observed(O)", "Expected(E)", "Contrib chi2");
        
        double chi2 = 0.0;
        int df_bins = 0;
        int num_bins = (180 + bw - 1) / bw;

        for (int b = 0; b < num_bins; ++b) {
            double O_bin = 0.0;
            double E_bin = 0.0;
            for (int i = 0; i < bw; ++i) {
                int idx = b * bw + i;
                if (idx < 180) { 
                    O_bin += real_counts_1deg[idx];
                    E_bin += expected_1deg[idx];
                }
            }
            
            if (E_bin > 0) {
                double contrib = ((O_bin - E_bin) * (O_bin - E_bin)) / E_bin;
                chi2 += contrib;
                df_bins++;
            }
            printf("Bin %2d    | %-12.0f | %-12.2f | %-10.4f\n", b, O_bin, E_bin, (E_bin > 0 ? ((O_bin - E_bin) * (O_bin - E_bin)) / E_bin : 0.0));
        }
        
        int df = std::max(1, df_bins - 1); 
        double reduced_chi2 = chi2 / df;
        double p_value = igamc(df / 2.0, chi2 / 2.0);

        final_report.push_back({bw, chi2, df, reduced_chi2, p_value});
        if (out) fprintf(out, "%d %.2f %d %.4f %e\n", bw, chi2, df, reduced_chi2, p_value);
        
        printf("RESULT FOR %d deg: Chi2 = %.2f, DF = %d, RedChi2 = %.4f, P-Value = %e\n", bw, chi2, df, reduced_chi2, p_value);
    }
    if (out) fclose(out);

    // ================= FINAL SUMMARY =================
    printf("\n============================================================\n");
    printf("FINAL SUMMARY\n");
    printf("============================================================\n");
    printf("%-10s | %-10s | %-5s | %-15s | %-10s\n", "Bin (deg)", "Chi^2", "DF", "Reduced Chi^2", "P-Value");
    for (auto& s : final_report) {
        printf("%-10d | %-10.2f | %-5d | %-15.4f | %e\n", s.size, s.chi2, s.df, s.red_chi2, s.p_value);
    }

    std::sort(final_report.begin(), final_report.end(), [](const Summary& a, const Summary& b) {
        return std::abs(a.red_chi2 - 1.0) < std::abs(b.red_chi2 - 1.0);
    });

    printf("\n============================================================\n");
    printf("FINAL ANSWER: OPTIMAL BINS (Closest to 1.0)\n");
    printf("============================================================\n");
    for(int i = 0; i < 3 && i < (int)final_report.size(); ++i) {
        printf("%d. Width: %d degrees (Chi2/df = %.4f, p-value = %.3f)\n", 
               i+1, final_report[i].size, final_report[i].red_chi2, final_report[i].p_value);
    }

    // ================= DETAILED ANALYSIS FOR BEST BIN SIZE =================
    int best_bw        = final_report[0].size;
    double best_chi2   = final_report[0].chi2;
    int best_df        = final_report[0].df;
    double best_red    = final_report[0].red_chi2;
    double best_pval   = final_report[0].p_value;

    printf("\n============================================================\n");
    printf("DETAILED BIN-BY-BIN ANALYSIS FOR BEST BIN SIZE: %d degrees\n", best_bw);
    printf("============================================================\n");
    printf("%-8s | %-10s | %-10s | %-12s | %-12s | %-10s | %-10s\n",
           "Bin No.", "Angle_from", "Angle_to", "Observed(O)", "Theory(E)", "O-E", "Contrib");
    printf("%s\n", std::string(85, '-').c_str());

    FILE* fdetail = fopen(OUTPUT_BEST_ANGLE_ANALYSIS.c_str(), "w");
    if (!fdetail) {
        printf("[ERROR] Nie mozna utworzyc pliku: %s\n", OUTPUT_BEST_ANGLE_ANALYSIS.c_str());
        printf("[ERROR] Sprawdz czy folder istnieje!\n");
    }
    if (fdetail) {
        fprintf(fdetail, "# Detailed bin-by-bin Chi2 analysis for best bin size\n");
        fprintf(fdetail, "# Radius: %.1f m | Best bin size: %d deg | Chi2: %.4f | DF: %d | Red.Chi2: %.4f | P-value: %e\n",
                TEST_RADIUS_M, best_bw, best_chi2, best_df, best_red, best_pval);
        fprintf(fdetail, "#\n");
        fprintf(fdetail, "%-8s %-10s %-10s %-12s %-12s %-12s %-10s\n",
                "Bin", "Angle_from", "Angle_to", "Observed_O", "Theory_E", "O_minus_E", "Chi2_contrib");
    }

    double chi2_recheck = 0.0;
    int num_bins_best = (180 + best_bw - 1) / best_bw;

    for (int b = 0; b < num_bins_best; ++b) {
        double O_bin = 0.0, E_bin = 0.0;
        int angle_from = b * best_bw + 1;
        int angle_to   = std::min((b + 1) * best_bw, 180);
        for (int i = 0; i < best_bw; ++i) {
            int idx = b * best_bw + i;
            if (idx < 180) {
                O_bin += real_counts_1deg[idx];
                E_bin += expected_1deg[idx];
            }
        }
        double contrib = 0.0;
        if (E_bin > 0) {
            contrib = ((O_bin - E_bin) * (O_bin - E_bin)) / E_bin;
            chi2_recheck += contrib;
        }
        double O_minus_E = O_bin - E_bin;

        printf("Bin %2d   | %4d deg    | %4d deg    | %-12.0f | %-12.2f | %-+10.2f | %-10.4f\n",
               b, angle_from, angle_to, O_bin, E_bin, O_minus_E, contrib);

        if (fdetail) {
            fprintf(fdetail, "%-8d %-10d %-10d %-12.0f %-12.2f %-+12.2f %-10.4f\n",
                    b, angle_from, angle_to, O_bin, E_bin, O_minus_E, contrib);
        }
    }

    printf("%s\n", std::string(85, '-').c_str());
    printf("TOTAL Chi2 = %.4f  (cross-check vs stored: %.4f)\n", chi2_recheck, best_chi2);

    if (fdetail) {
        fprintf(fdetail, "#\n");
        fprintf(fdetail, "# TOTAL Chi2 = %.4f\n", chi2_recheck);
        fclose(fdetail);
        printf("Detailed angle analysis saved to: %s\n", OUTPUT_BEST_ANGLE_ANALYSIS.c_str());
    }

    return 0;
}