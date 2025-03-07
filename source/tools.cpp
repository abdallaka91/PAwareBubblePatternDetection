#include "tools.h"
#include "struct.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
// #include "debugging_tools.h"

using std::vector;

void PoAwN::tools::circshift(vector<uint16_t> &CCSK_bin_seq, int shift)
{
    const uint16_t q = CCSK_bin_seq.size();
    const auto shift_masked = shift & (q - 1);
    rotate(CCSK_bin_seq.begin(), CCSK_bin_seq.begin() + (q - shift_masked), CCSK_bin_seq.end());
}

void PoAwN::tools::create_ccsk_rotated_table(const uint16_t *CCSK_bin_seq,
                                             const uint16_t q,
                                             vector<vector<uint16_t>> &ccsk_rotated_table)
{
    std::vector<uint16_t> eta0(CCSK_bin_seq, CCSK_bin_seq + q);
    ccsk_rotated_table[0] = eta0;

    for (int i = 1; i < q; ++i)
    {
        circshift(eta0, -1);
        ccsk_rotated_table[i] = eta0;
    }
}
void PoAwN::tools::RandomBinaryGenerator(const base_code_t code_param,
                                         const vector<vector<uint16_t>> &bin_table,
                                         const bool repeatable,
                                         const int SEED,
                                         vector<vector<uint16_t>> &KBIN,
                                         vector<uint16_t> &KSYMB)
{
    KSYMB.resize(code_param.K);
    KBIN.resize(code_param.K);
    for (uint16_t k = 0; k < code_param.K; k++)
    {
        static std::mt19937 gen(repeatable ? SEED : std::random_device{}());
        static std::uniform_int_distribution<int> dist(0, code_param.q - 1);
        uint16_t randv = (uint16_t)dist(gen);
        KSYMB[k] = randv;
        KBIN[k] = bin_table[randv];
    }
}

void PoAwN::tools::AWGN_gen(double MEAN,
                            double STD,
                            bool repeatable,
                            double SEED,
                            vector<vector<softdata_t>> &noise_table)
{
    std::mt19937 gen(repeatable ? std::mt19937(SEED) : std::mt19937(std::random_device{}()));
    int N = noise_table.size();
    for (int i = 0; i < N; i++)
    {
        int q1 = noise_table[i].size();
        for (int j = 0; j < q1; j++)
        {
            std::normal_distribution<double> dist(MEAN, STD);
            noise_table[i][j] = (softdata_t)dist(gen);
        }
    }
}

void PoAwN::tools::awgn_channel_noise(const vector<vector<uint16_t>> &NBIN,
                                      const double sigma,
                                      const bool repeatable,
                                      const double SEED,
                                      vector<vector<softdata_t>> &noisy_sig)
{
    uint16_t N = NBIN.size();
    uint16_t q1 = NBIN[0].size();
    vector<vector<softdata_t>> noise_table;
    noise_table.resize(N, vector<softdata_t>(q1, 0));
    AWGN_gen(SEED, sigma, repeatable, SEED, noise_table);
        for (int i = 0; i < N; i++)
        for (int j = 0; j < q1; j++)
        {
            noisy_sig[i][j] = (NBIN[i][j] == 0) ? 1 : -1;
            noisy_sig[i][j] += noise_table[i][j];
        }
}

void PoAwN::tools::Encoder(const vector<vector<uint16_t>> &ADDDEC, const vector<vector<uint16_t>> &MULDEC,
                           const vector<vector<uint16_t>> &polar_coeff,
                           const vector<uint16_t> u_symb, vector<uint16_t> &NSYMB)
{
    uint16_t N = u_symb.size();
    uint16_t n = log2(N);
    vector<vector<uint16_t>> temp_symb(N, vector<uint16_t>(n + 1, 0));
    for (uint16_t i = 0; i < N; ++i)
        temp_symb[i][0] = u_symb[i];

    uint16_t a, b;
    uint16_t tmp_add;
    uint16_t tmp_mul;
    for (uint16_t l = 1; l <= n; l++)
    {
        for (uint16_t k = 0; k < N / 2; k++)
        {
            uint16_t pw1 = pow(2, l - 1);
            a = 2 * k - k % pw1;
            b = pw1 + 2 * k - k % pw1;
            tmp_add = ADDDEC[temp_symb[a][l - 1]][temp_symb[b][l - 1]];
            temp_symb[a][l] = tmp_add;
            tmp_mul = MULDEC[temp_symb[b][l - 1]][polar_coeff[l - 1][k]];
            temp_symb[b][l] = tmp_mul;
        }
    }
    for (uint16_t i = 0; i < N; ++i)
        NSYMB[i] = temp_symb[i][n];
}