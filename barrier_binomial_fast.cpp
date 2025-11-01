//barrier_binomial_fast.cpp
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <cassert>
#include <algorithm>
#include <omp.h>
#include <iomanip>

using namespace std;

using Clock = chrono::high_resolution_clock;
using TimePoint = Clock::time_point;

// ------------------------------------------------------------------
// 🌟 Fast Binomial Summation Formulas for computing some barrier option prices created by Dr. Burak YILDIZ, Ankara/Türkiye
//     Email: tr.burak.yildiz@gmail.com
//     Google Scholar: https://scholar.google.com/citations?user=DFV7KBQAAAAJ&hl=tr
//     Example How to compile:
//     g++-15 barrier_binomial_fast.cpp -O3 -march=native -ffast-math -funroll-loops -fopenmp -o barrier_binomial_fast
//     Example run:
//     ./barrier_binomial_fast 2000
// ------------------------------------------------------------------
// Google gemini and ChatGPT is used for code optimization
// ------------------------------------------------------------------
// 🌟 Log Binomial Precomputation (Stable)
// ------------------------------------------------------------------
// Calculates ln(C(n, j)) for j = 0 to n using log factorials
static vector<double> precompute_log_C(int n) {
    // 1. Pre-calculate ln(i!) for i = 0 to n
    vector<double> log_factorial(n + 1);
    log_factorial[0] = 0.0; // ln(0!) = ln(1) = 0
    for (int i = 1; i <= n; ++i) {
        log_factorial[i] = log_factorial[i - 1] + std::log(double(i));
    }

    // 2. Calculate ln(C(n, j)) = ln(n!) - ln(j!) - ln((n-j)!)
    vector<double> log_C(n + 1);
    for (int j = 0; j <= n; ++j) {
        log_C[j] = log_factorial[n] - log_factorial[j] - log_factorial[n - j];
    }
    return log_C;
}

// ------------------------------------------------------------------
// ⚡ Optimized & Stable Vanilla Put
// ------------------------------------------------------------------
double myputprice(int n, double S, double K, double r, double v, double T, const vector<double>& log_C) {
    double dt = T / n;
    double u = exp(v * sqrt(dt));
    double d = exp(-v * sqrt(dt));
    double p = (exp(r * dt) - d) / (u - d);
    double q = 1.0 - p;
    double disc = exp(-r * T);
    
    double ln_p = std::log(p);
    double ln_q = std::log(q);
    double ln_u = std::log(u);
    double ln_d = std::log(d);
    double ln_S = std::log(S);
    
    double beta = (std::log(K / S) - n * ln_d) / (ln_u - ln_d);
    
    double sum = 0.0;
    int end_j = (int)ceil(beta);
    if (end_j < 0) end_j = 0;
    if (end_j > n) end_j = n;

    #pragma omp parallel for reduction(+:sum)
    for (int j = 0; j <= end_j; ++j) {
        double log_coeff = log_C[j];
        
        // Calculate log(Probability)
        double log_prob = log_coeff + j * ln_p + (n - j) * ln_q;
        
        // Calculate Stock Price (ST)
        double log_ST = ln_S + j * ln_u + (n - j) * ln_d;
        double ST = std::exp(log_ST);
        
        if (ST < K) {
            double prob = std::exp(log_prob);
            sum += prob * (K - ST);
        }
    }
    return disc * sum;
}

// ------------------------------------------------------------------
// ⚡ Optimized & Stable Vanilla Call
// ------------------------------------------------------------------
double mycallprice(int n, double S, double K, double r, double v, double T, const vector<double>& log_C) {
    double dt = T / n;
    double u = exp(v * sqrt(dt));
    double d = exp(-v * sqrt(dt));
    double p = (exp(r * dt) - d) / (u - d);
    double q = 1.0 - p;
    double disc = exp(-r * T);

    double ln_p = std::log(p);
    double ln_q = std::log(q);
    double ln_u = std::log(u);
    double ln_d = std::log(d);
    double ln_S = std::log(S);

    double alfa = (n / 2.0) - (std::log(S / K) / (2.0 * v * sqrt(dt)));
    
    double sum = 0.0;
    int start_j = (int)ceil(alfa);
    if (start_j < 0) start_j = 0;
    if (start_j > n) return 0.0;

    #pragma omp parallel for reduction(+:sum)
    for (int j = start_j; j <= n; ++j) {
        double log_coeff = log_C[j];
        
        // Calculate log(Probability)
        double log_prob = log_coeff + j * ln_p + (n - j) * ln_q;
        
        // Calculate Stock Price (ST)
        double log_ST = ln_S + j * ln_u + (n - j) * ln_d;
        double ST = std::exp(log_ST);
        
        if (ST > K) {
            double prob = std::exp(log_prob);
            sum += prob * (ST - K);
        }
    }
    return disc * sum;
}

// ------------------------------------------------------------------
// 🌟 ENHANCED STABILITY: Down-and-Out Put
// ------------------------------------------------------------------
double downoutputprice(int n, double Barrier, double S, double K, double r, double v, double T, const vector<double>& log_C) {
    double dt = T / n;
    double u = exp(v * sqrt(dt));
    double d = exp(-v * sqrt(dt));
    double p = (exp(r * dt) - d) / (u - d);
    double q = 1.0 - p;
    double disc = exp(-r * T);

    double ln_p = std::log(p);
    double ln_q = std::log(q);
    double ln_u = std::log(u);
    double ln_d = std::log(d);
    double ln_S = std::log(S);

    double beta = (std::log(K / S) - n * ln_d) / (ln_u - ln_d);

    // Reflection index H
    int H = 0;
    double Sd = S * d;
    for (int i = 1; i <= n; ++i, Sd *= d)
        if (Sd < Barrier) { H = -i; break; }

    double sum = 0.0;
    int end_j = (int)ceil(beta);
    if (end_j < 0) end_j = 0;
    if (end_j > n) end_j = n;

    #pragma omp parallel for reduction(+:sum)
    for (int j = 0; j <= end_j; ++j) {
        double log_coeff_main = log_C[j];
        
        // Reflected index for down-and-out options (j - H)
        int idx_sub = j - H;
        double log_coeff_sub = (idx_sub >= 0 && idx_sub <= n) ? log_C[idx_sub] : -INFINITY;

        // Apply Reflection Principle in LOG space for stability
        double log_coeff_ratio = log_coeff_sub - log_coeff_main; // ln(B/A)
        
        if (log_coeff_ratio < 0.0) { // log(B) < log(A) implies A > B (valid path: A - B > 0)
            
            // Stable calculation of ln(A - B) = ln(A) + ln(1 - exp(ln(B) - ln(A)))
            double log_effective_coeff = log_coeff_main + std::log(1.0 - std::exp(log_coeff_ratio));

            // Calculate log(Probability) (excluding C)
            double log_prob_no_C = j * ln_p + (n - j) * ln_q;
            
            // Combined log-probability
            double log_final_prob = log_effective_coeff + log_prob_no_C;
            
            // Calculate Stock Price (ST)
            double log_ST = ln_S + j * ln_u + (n - j) * ln_d;
            double ST = std::exp(log_ST);
            
            if (ST < K) {
                double final_prob = std::exp(log_final_prob);
                sum += final_prob * (K - ST);
            }
        }
    }
    return disc * sum;
}

// ------------------------------------------------------------------
// 🌟 ENHANCED STABILITY: Down-and-Out Call
// ------------------------------------------------------------------
double downoutcallprice(int n, double Barrier, double S, double K, double r, double v, double T, const vector<double>& log_C) {
    double dt = T / n;
    double u = exp(v * sqrt(dt));
    double d = exp(-v * sqrt(dt));
    double p = (exp(r * dt) - d) / (u - d);
    double q = 1.0 - p;
    double disc = exp(-r * T);

    double ln_p = std::log(p);
    double ln_q = std::log(q);
    double ln_u = std::log(u);
    double ln_d = std::log(d);
    double ln_S = std::log(S);
    
    double alfa = (n / 2.0) - (std::log(S / K) / (2.0 * v * sqrt(dt)));

    // Reflection index H
    int H = 0;
    double Sd = S * d;
    for (int i = 1; i <= n; ++i, Sd *= d)
        if (Sd < Barrier) { H = -i; break; }

    double sum = 0.0;
    int start_j = (int)ceil(alfa);
    if (start_j < 0) start_j = 0;
    if (start_j > n) return 0.0;

    #pragma omp parallel for reduction(+:sum)
    for (int j = start_j; j <= n; ++j) {
        double log_coeff_main = log_C[j];
        
        // Reflected index for down-and-out options (j - H)
        int idx_sub = j - H;
        double log_coeff_sub = (idx_sub >= 0 && idx_sub <= n) ? log_C[idx_sub] : -INFINITY;

        // Apply Reflection Principle in LOG space for stability
        double log_coeff_ratio = log_coeff_sub - log_coeff_main; // ln(B/A)
        
        if (log_coeff_ratio < 0.0) { // log(B) < log(A) implies A > B (valid path)
            
            // Stable calculation of ln(A - B)
            double log_effective_coeff = log_coeff_main + std::log(1.0 - std::exp(log_coeff_ratio));

            // Calculate log(Probability) (excluding C)
            double log_prob_no_C = j * ln_p + (n - j) * ln_q;
            
            // Combined log-probability
            double log_final_prob = log_effective_coeff + log_prob_no_C;
            
            // Calculate Stock Price (ST)
            double log_ST = ln_S + j * ln_u + (n - j) * ln_d;
            double ST = std::exp(log_ST);
            
            if (ST > K) {
                double final_prob = std::exp(log_final_prob);
                sum += final_prob * (ST - K);
            }
        }
    }
    return disc * sum;
}

// ------------------------------------------------------------------
// 🌟 ENHANCED STABILITY: Up-and-Out Call
// ------------------------------------------------------------------
double upoutcallprice(int n, double Barrier, double S, double K, double r, double v, double T, const vector<double>& log_C) {
    double dt = T / n;
    double u = exp(v * sqrt(dt));
    double d = exp(-v * sqrt(dt));
    double p = (exp(r * dt) - d) / (u - d);
    double q = 1.0 - p;
    double disc = exp(-r * T);

    double ln_p = std::log(p);
    double ln_q = std::log(q);
    double ln_u = std::log(u);
    double ln_d = std::log(d);
    double ln_S = std::log(S);

    double alfa = (n / 2.0) - (std::log(S / K) / (2.0 * v * sqrt(dt)));

    // Reflection index H
    int H = 0;
    double Su = S * u;
    for (int i = 1; i <= n; ++i, Su *= u)
        if (Su > Barrier) { H = i; break; }

    double sum = 0.0;
    int start_j = (int)ceil(alfa);
    if (start_j < 0) start_j = 0;
    if (start_j > n) return 0.0;

    #pragma omp parallel for reduction(+:sum)
    for (int j = start_j; j <= n; ++j) {
        // Main coefficient is C(n, n-j)
        int idx_main = n - j;
        double log_coeff_main = (idx_main >= 0 && idx_main <= n) ? log_C[idx_main] : -INFINITY;

        // Reflected coefficient is C(n, n-j + H)
        int idx_sub = n - j + H;
        double log_coeff_sub = (idx_sub >= 0 && idx_sub <= n) ? log_C[idx_sub] : -INFINITY;

        // Apply Reflection Principle in LOG space for stability
        double log_coeff_ratio = log_coeff_sub - log_coeff_main; // ln(B/A)
        
        if (log_coeff_ratio < 0.0) { // log(B) < log(A) implies A > B (valid path)
            
            // Stable calculation of ln(A - B)
            double log_effective_coeff = log_coeff_main + std::log(1.0 - std::exp(log_coeff_ratio));

            // Calculate log(Probability) (excluding C)
            double log_prob_no_C = j * ln_p + (n - j) * ln_q;
            
            // Combined log-probability
            double log_final_prob = log_effective_coeff + log_prob_no_C;
            
            // Calculate Stock Price (ST)
            double log_ST = ln_S + j * ln_u + (n - j) * ln_d;
            double ST = std::exp(log_ST);
            
            if (ST > K) {
                double final_prob = std::exp(log_final_prob);
                sum += final_prob * (ST - K);
            }
        }
    }
    return disc * sum;
}

// ------------------------------------------------------------------
// 🌟 ENHANCED STABILITY: Up-and-Out Put
// ------------------------------------------------------------------
double upoutputprice(int n, double Barrier, double S, double K, double r, double v, double T, const vector<double>& log_C) {
    double dt = T / n;
    double u = exp(v * sqrt(dt));
    double d = exp(-v * sqrt(dt));
    double p = (exp(r * dt) - d) / (u - d);
    double q = 1.0 - p;
    double disc = exp(-r * T);

    double ln_p = std::log(p);
    double ln_q = std::log(q);
    double ln_u = std::log(u);
    double ln_d = std::log(d);
    double ln_S = std::log(S);

    double beta = (std::log(K / S) - n * ln_d) / (ln_u - ln_d);

    // Reflection index H
    int H = 0;
    double Su = S * u;
    for (int i = 1; i <= n; ++i, Su *= u)
        if (Su > Barrier) { H = i; break; }

    double sum = 0.0;
    int end_j = (int)ceil(beta);
    if (end_j < 0) end_j = 0;
    if (end_j > n) end_j = n;

    #pragma omp parallel for reduction(+:sum)
    for (int j = 0; j <= end_j; ++j) {
        // Main coefficient is C(n, n-j)
        int idx_main = n - j;
        double log_coeff_main = (idx_main >= 0 && idx_main <= n) ? log_C[idx_main] : -INFINITY;

        // Reflected coefficient is C(n, n-j + H)
        int idx_sub = n - j + H;
        double log_coeff_sub = (idx_sub >= 0 && idx_sub <= n) ? log_C[idx_sub] : -INFINITY;

        // Apply Reflection Principle in LOG space for stability
        double log_coeff_ratio = log_coeff_sub - log_coeff_main; // ln(B/A)
        
        if (log_coeff_ratio < 0.0) { // log(B) < log(A) implies A > B (valid path)
            
            // Stable calculation of ln(A - B)
            double log_effective_coeff = log_coeff_main + std::log(1.0 - std::exp(log_coeff_ratio));

            // Calculate log(Probability) (excluding C)
            double log_prob_no_C = j * ln_p + (n - j) * ln_q;
            
            // Combined log-probability
            double log_final_prob = log_effective_coeff + log_prob_no_C;
            
            // Calculate Stock Price (ST)
            double log_ST = ln_S + j * ln_u + (n - j) * ln_d;
            double ST = std::exp(log_ST);
            
            if (ST < K) {
                double final_prob = std::exp(log_final_prob);
                sum += final_prob * (K - ST);
            }
        }
    }
    return disc * sum;
}


// ------------------------------------------------------------------
// Main function with separate timing
// ------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc <= 1) { cerr << "Usage: " << argv[0] << " <n>\n"; return 2; }
    int n = atoi(argv[1]);
    assert(n > 0);
    
    // Set output precision
    cout << fixed << setprecision(8);

    double S = 90.0, K = 100.0, Bdown = 70.0, Bup = 120.0;
    double r = 0.05, sigma = 0.25, T = 0.5;

    cout << "-----------Parameters------------------------------\n";
    cout << "Underlying: " << S << "  Strike: " << K
         << "  Barrier Up: " << Bup << "  Barrier Down: " << Bdown
         << "  r: " << r << "  sigma: " << sigma << "  T: " << T 
         << "  Steps (n): " << n << "\n";
    cout << "---------------------------------------------------\n";

    // 🌟 STABILITY FIX: Precompute log(C(n, j)) ONCE
    TimePoint start_setup = Clock::now();
    auto log_C = precompute_log_C(n);
    TimePoint end_setup = Clock::now();
    
    // --- Option Pricing and Timing ---
    double c, p, cuo, puo, cdo, pdo;
    double t_c, t_p, t_cuo, t_puo, t_cdo, t_pdo;

    // Vanilla Call
    TimePoint start = Clock::now();
    c = mycallprice(n, S, K, r, sigma, T, log_C);
    TimePoint end = Clock::now();
    t_c = chrono::duration<double, milli>(end - start).count();

    // Vanilla Put
    start = Clock::now();
    p = myputprice(n, S, K, r, sigma, T, log_C);
    end = Clock::now();
    t_p = chrono::duration<double, milli>(end - start).count();

    // Up-Out Call
    start = Clock::now();
    cuo = upoutcallprice(n, Bup, S, K, r, sigma, T, log_C);
    end = Clock::now();
    t_cuo = chrono::duration<double, milli>(end - start).count();

    // Up-Out Put
    start = Clock::now();
    puo = upoutputprice(n, Bup, S, K, r, sigma, T, log_C);
    end = Clock::now();
    t_puo = chrono::duration<double, milli>(end - start).count();

    // Down-Out Call
    start = Clock::now();
    cdo = downoutcallprice(n, Bdown, S, K, r, sigma, T, log_C);
    end = Clock::now();
    t_cdo = chrono::duration<double, milli>(end - start).count();

    // Down-Out Put
    start = Clock::now();
    pdo = downoutputprice(n, Bdown, S, K, r, sigma, T, log_C);
    end = Clock::now();
    t_pdo = chrono::duration<double, milli>(end - start).count();

    // --- Results Report ---
    cout << "Option              | Price (USD) | Elapsed (ms) \n";
    cout << "--------------------|-------------|--------------\n";
    cout << "Vanilla Call        | " << c     << " | " << t_c   << "\n";
    cout << "Vanilla Put         | " << p     << " | " << t_p   << "\n";
    cout << "Up-Out Call (UOC)   | " << cuo   << " | " << t_cuo << "\n";
    cout << "Up-Out Put (UOP)    | " << puo   << " | " << t_puo << "\n";
    cout << "Down-Out Call (DOC) | " << cdo   << " | " << t_cdo << "\n";
    cout << "Down-Out Put (DOP)  | " << pdo   << " | " << t_pdo << "\n";
    cout << "--------------------|-------------|--------------\n";
    cout << "Setup/Coefficient Time: " << chrono::duration<double, milli>(end_setup - start_setup).count() << " ms\n";
}
