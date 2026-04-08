#include "his_1D.C"

void generate_all_histograms() {

    std::vector<int> datasets = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; // ex. {5,6,7}
    std::vector<std::string> coins = {"3", "3f", "4", "34", "34f", "false", "true"};
    std::vector<int> bin_widths = {3, 6, 9, 12, 15, 18};

    for (int x : datasets) {
        for (const auto& coin : coins) {
            for (int bw : bin_widths) {

/*                
                //put here your absolute path to the respository in your PC
                std::string input =
                    "D:/MuonDensityMapping/output/dane_" + std::to_string(x) + "/processed/coin_" + std::to_string(x) + "_" + coin + ".txt";

                std::string output =
                    "D:/MuonDensityMapping/results/single/dane_" + std::to_string(x) + "/coin_" + coin + "/hist_" + coin + "_" + std::to_string(bw) + "deg.pdf";
*/                


                // if needed to generate data for whole
                std::string input =
                    "D:/MuonDensityMapping/output/all/zussamen_" + std::to_string(x) + "/coin_" + coin + ".txt";

                std::string output =
                    "D:/MuonDensityMapping/results/all/zussamen_" + std::to_string(x) + "/coin_" + coin + "/hist_" + coin + "_" + std::to_string(bw) + "deg.pdf";


                his_1D(
                    input,output,"Coin " + coin + " / Dane " + std::to_string(x), bw,  3
                );
            }
        }
    }
}
