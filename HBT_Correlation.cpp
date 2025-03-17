#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <omp.h>
#include <algorithm>
#include <limits>
#include <stdexcept>

using namespace std;

// Constants
const int TotalEvents = 25000; 
const int MixedEvents = 200;
const int MaxParticles = 30000;
const int Q_num_bins = 51;
const int kT_num_bins = 4;
const int mT_num_bins = 1;

const double Q_min = -0.2, Q_max = 0.2;
const double kT_min = 0.01, kT_max = 2.1;
const double mT_min = 0.856, mT_max = 0.940;
const double rap_min = -1.0, rap_max = 1.0, buffer_rap = 0.1;
const double Q_bin_width = (Q_max - Q_min) / Q_num_bins;
const double kT_bin_width = (kT_max - kT_min) / kT_num_bins;
const double mT_bin_width = (mT_max - mT_min) / mT_num_bins;
const double pion_plus_mass = 0.139570;
const double kaon_plus_mass = 0.493677;
const double hbarC = 0.197326980;

// Particle struct to store particle data
struct Particle {
    double e, px, py, pz, x, y, z, t;
};

// Function to get PID from UrQMD MCID
int get_PID_from_urqmd_MCID(int mcid, int iso) {
    if (mcid == 101 && iso == 2) return 211;   // pi^{+}
    if (mcid == 101 && iso == 0) return 111;   // pi^{0}
    if (mcid == 101 && iso == -2) return -211; // pi^{-}
    if (mcid == 106 && iso == 1) return 321;   // k^{+}
    if (mcid == 106 && iso == -1) return 311;  // k^{0}
    if (mcid == -106 && iso == 1) return -311; // anti k^{0}
    if (mcid == -106 && iso == -1) return -321; // k^{-}
    if (mcid == 109 && iso == 0) return 333;   // phi(1020)
    if (mcid == 102 && iso == 0) return 221;   // eta
    if (mcid == 100 && iso == 0) return 22;    // photon
    if (mcid == 1 && iso == 1) return 2212;    // p
    if (mcid == 1 && iso == -1) return 2112;   // n
    if (mcid == -1 && iso == -1) return -2212; // anti p
    if (mcid == -1 && iso == 1) return -2112;  // anti n
    if (mcid == 40 && iso == 2) return 3222;   // Sigma^{+}
    if (mcid == -40 && iso == -2) return -3222; // anti Sigma^{-}
    if (mcid == 40 && iso == 0) return 3212;   // Sigma^{0}
    if (mcid == -40 && iso == 0) return -3212; // anti Sigma^{0}
    if (mcid == 40 && iso == -2) return 3112;  // Sigma^{-}
    if (mcid == -40 && iso == 2) return -3112; // anti Sigma^{+}
    if (mcid == 49 && iso == -1) return 3312;  // Xi^{-}
    if (mcid == 49 && iso == 1) return 3322;   // Xi^{0}
    if (mcid == -49 && iso == -1) return -3322; // anti Xi^{0}
    if (mcid == -49 && iso == 1) return -3312; // anti Xi^{+}
    if (mcid == 27 && iso == 0) return 3122;   // Lambda
    if (mcid == -27 && iso == 0) return -3122; // anti Lambda
    if (mcid == 55 && iso == 0) return 3334;   // Omega
    if (mcid == -55 && iso == 0) return -3334; // anti Omega
    return 0;
}

// Function to read particle data from file
std::vector<std::vector<Particle>> readParticleData(const std::string& filename, int totalEvents) {
    std::ifstream file(filename, std::ios::binary | std::ios::in);
    if (!file) {
        throw std::runtime_error("File not found: " + filename);
    }

    std::vector<std::vector<Particle>> particles(totalEvents);
    for (int event = 0; event < totalEvents; ++event) {
        int nch, count = 0;
        file.read(reinterpret_cast<char*>(&nch), sizeof(int));
        if (file.eof()) break;

        particles[event].reserve(nch);
        for (int i = 0; i < nch; ++i) {
            float particle_array[11];
            file.read(reinterpret_cast<char*>(particle_array), sizeof(float) * 11);

            int pid = get_PID_from_urqmd_MCID(static_cast<int>(particle_array[1]), static_cast<int>(particle_array[2]));
            if (pid == -2212) {  // CHANGE THE PID FOR EVERY PARTICLE
                count++;
		Particle p;
                p.e = particle_array[7];
                p.px = particle_array[8];
                p.py = particle_array[9];
                p.pz = particle_array[10];
                p.x = particle_array[4];
                p.y = particle_array[5];
                p.z = particle_array[6];
                p.t = particle_array[3];
                particles[event].push_back(p);
            }
        } 
	if (event%100 == 0) std::cout<< "Particle count for event " << event << " is:   " << count << std::endl ;
    }
    return particles;
}

// Function to calculate correlation function
void calculateCorrelation(const std::vector<std::vector<Particle>>& particles, int totalEvents, int mixedEvents) {
    std::vector<std::vector<std::vector<std::vector<double>>>> numerator(kT_num_bins,
        std::vector<std::vector<std::vector<double>>>(Q_num_bins,
        std::vector<std::vector<double>>(Q_num_bins,
        std::vector<double>(Q_num_bins, 0.0))));

    std::vector<std::vector<std::vector<std::vector<double>>>> denominator(kT_num_bins,
        std::vector<std::vector<std::vector<double>>>(Q_num_bins,
        std::vector<std::vector<double>>(Q_num_bins,
        std::vector<double>(Q_num_bins, 0.0))));

    #pragma omp parallel for
    for (int event = 0; event < totalEvents; ++event) {
        for (int mixev = 0; mixev < mixedEvents; ++mixev) {
            if (mixev == event) continue;

            for (const auto& p1 : particles[event]) {
                for (const auto& p2 : particles[mixev]) {
                
                    double Kx = p1.px + p2.px;
                    double Ky = p1.py + p2.py;
                    double Kz = p1.pz + p2.pz;
                    double Ke = p1.e + p2.e;

                    double K_perp_sq = Kx * Kx + Ky * Ky;
                    double K_perp = std::sqrt(K_perp_sq);
                    double Mt_sq = Ke * Ke - Kz * Kz;
                    double Mt = std::sqrt(Mt_sq);

                    double ratio = Kz / Ke;

                    if (K_perp < kT_min || K_perp > kT_max) continue;

                    int Kperp_index = static_cast<int>((K_perp - kT_min) / kT_bin_width);
                    if (Kperp_index >= kT_num_bins) continue;

                    // Transform to LCMS frame
                    double tBeta = ratio;
                    double tGamma = Ke / Mt;
                    double pz1LCMS = tGamma * (p1.pz - tBeta * p1.e);
                    double e1LCMS = tGamma * (p1.e - tBeta * p1.pz);

                    double qx = p1.px - p2.px;
                    double qy = p1.py - p2.py;
                    double qz = pz1LCMS - tGamma * (p2.pz - tBeta * p2.e);
                    double qe = e1LCMS - tGamma * (p2.e - tBeta * p2.pz);

                    // Bin the relative momentum
                    int qout_idx = static_cast<int>((qx - Q_min) / Q_bin_width);
                    int qside_idx = static_cast<int>((qy - Q_min) / Q_bin_width);
                    int qlong_idx = static_cast<int>((qz - Q_min) / Q_bin_width);

                    if (qout_idx < 0 || qout_idx >= Q_num_bins || qside_idx < 0 || 
                    qside_idx >= Q_num_bins || qlong_idx < 0 || qlong_idx >= Q_num_bins) continue;

                    double cosqx = std::cos(qe * (p1.t - p2.t) / hbarC - qx * (p1.x - p2.x) 
                    / hbarC - qy * (p1.y - p2.y) / hbarC - qz * (p1.z - p2.z) / hbarC);

                    #pragma omp atomic
                    numerator[Kperp_index][qout_idx][qside_idx][qlong_idx] += (1 + cosqx);

                    #pragma omp atomic
                    denominator[Kperp_index][qout_idx][qside_idx][qlong_idx] += 1.0;
                }
            }
        }
        if (event%10 == 0) std::cout << "Calculations for Event Number:  " << event << "  done." << std::endl ;
    }

    // Write output files
    for (int kTbin = 0; kTbin < kT_num_bins; ++kTbin) { // CHANGE THE FILENAME FOR EVERY PARTICLE PAIR
        std::string filename = "Correlation_Function_Dataset_Antiproton_etafall_0.5_" + std::to_string(kTbin) + ".data";
        std::ofstream file(filename);
        for (int ii = 0; ii < Q_num_bins; ++ii) {
            for (int jj = 0; jj < Q_num_bins; ++jj) {
                for (int kk = 0; kk < Q_num_bins; ++kk) {
                    if (denominator[kTbin][ii][jj][kk] == 0) continue;

                    double correlfunc = numerator[kTbin][ii][jj][kk] / denominator[kTbin][ii][jj][kk];
                    double delta_correlfunc = std::sqrt((1.0 / numerator[kTbin][ii][jj][kk]) + (1.0 / denominator[kTbin][ii][jj][kk]));

                    file << Q_min + ii * Q_bin_width << " " << Q_min + jj * Q_bin_width << " " << Q_min + kk * Q_bin_width << " "
                         << correlfunc << " " << delta_correlfunc << "\n";
                }
            }
        }
        file.close();
    }
}

int main() {
    std::string filename = "particle_list_0_10_etafall_0.5_200GeV.bin";
    auto particles = readParticleData(filename, TotalEvents);
    calculateCorrelation(particles, TotalEvents, MixedEvents);
    std::cout << "Data files generated successfully." << std::endl;
    return 0;
}
