#include "His_1D.C"

void generate_all_histograms() {

    std::vector<int> datasets = {7}; // ex. {5,6,7}
    std::vector<std::string> coins = {"3", "3f", "4", "34", "34f"};
    std::vector<int> bin_widths = {3, 9, 15};

    for (int x : datasets) {
        for (const auto& coin : coins) {
            for (int bw : bin_widths) {

                
                //put here your absolute path to the respository in your PC
                std::string input =
                    "D:/MuonDensityMapping/Output/dane_" + std::to_string(x) + "/processed/coin_" + std::to_string(x) + "_" + coin + ".txt";

                std::string output =
                "D:/MuonDensityMapping/results/Dane_" + std::to_string(x) + "/coin_" + coin + "/hist_" + coin + "_" + std::to_string(bw) + "deg.pdf";
                

                /*
                // if needed to generate data for whole
                std::string input =
                    "D:/MuonDensityMapping/Output/all/4.january/coin_" + coin + ".txt";

                std::string output =
                "D:/MuonDensityMapping/results/all/4.january/coin_" + coin + "/hist_" + coin + "_" + std::to_string(bw) + "deg.pdf";
                */

                His_1D(
                    input,output,"Coin " + coin + " / Dane " + std::to_string(x), bw,  3
                );
            }
        }
    }
}
