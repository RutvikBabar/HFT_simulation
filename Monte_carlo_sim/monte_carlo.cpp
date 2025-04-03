#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>

// Process class: Simulates one asset price path.
class Process {
public:
    double drift;
    double volatility;
    double delta_t;
    double current_price;
    std::vector<double> asset_price;

    Process(double drift, double volatility, double delta_t, double initial_price)
        : drift(drift), volatility(volatility), delta_t(delta_t), current_price(initial_price) {
        asset_price.push_back(initial_price);
    }

    // Advance the process one time step.
    // Pass in a random number generator and a standard normal distribution.
    void time_step(std::default_random_engine& generator, std::normal_distribution<double>& dist) {
        // Generate a normally distributed number with mean=0 and std=1, then scale by sqrt(delta_t)
        double dw = dist(generator) * std::sqrt(delta_t);
        double ds = drift * delta_t * current_price + volatility * current_price * dw;
        current_price += ds;
        asset_price.push_back(current_price);
    }
};

// StockSim class: Simulates many paths and computes the discounted average final price.
class StockSim {
public:
    double price;

    StockSim(int n_paths, double initial_asset_price, double drift, double delta_t,
        double volatility, double tte, double rfr) {
        std::vector<Process> processes;
        processes.reserve(n_paths);

        // Create n_paths processes.
        for (int i = 0; i < n_paths; i++) {
            processes.emplace_back(drift, volatility, delta_t, initial_asset_price);
        }

        int n_steps = static_cast<int>(tte / delta_t);

        // Setup random number generation.
        std::random_device rd;
        std::default_random_engine generator(rd());
        std::normal_distribution<double> dist(0.0, 1.0);

        // Simulate each path over n_steps.
        for (auto& process : processes) {
            for (int j = 0; j < n_steps; j++) {
                process.time_step(generator, dist);
            }
        }

        // Compute the average of the final prices.
        double sum_final_prices = 0.0;
        for (const auto& process : processes) {
            sum_final_prices += process.asset_price.back();
        }
        double avg_final_price = sum_final_prices / n_paths;

        // Discount the average final price back to present value.
        price = avg_final_price * std::exp(-tte * rfr);
    }
};

int main() {
    // Parameters from the S&P 500 example:
    int n_paths = 10000;
    double initial_asset_price = 5447.40; // Current index level
    double drift = 0.0457;               // ~4.57% annual drift
    double volatility = 0.18;            // ~18% annual volatility (assumption)
    double delta_t = 1.0 / 252.0;          // Daily steps (assuming 252 trading days)
    double tte = 1.0;                    // 1 year to expiry
    double rfr = 0.04501;                // 4.501% risk-free rate

    StockSim simulation(n_paths, initial_asset_price, drift, delta_t, volatility, tte, rfr);
    std::cout << "Simulated current fair value for the S&P 500: "
        << std::fixed << std::setprecision(2) << simulation.price << std::endl;

    return 0;
}
