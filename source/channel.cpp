#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "Decoder_functions.h"
#include "GF_tools.h"
#include "init.h"
#include "struct.h"
#include "tools.h"
#include "channel.h"
#define PI 3.14159265358979323846

using namespace PoAwN::structures;
using namespace PoAwN::tools;
using namespace PoAwN::init;
using namespace PoAwN::decoding;
using std::array;
using std::cout;
using std::endl;
using std::stod;
using std::stoi;
using std::string;
using std::vector;

void PoAwN::channel::EncodeChanBPSK_BinCCSK(const decoder_parameters &dec_param,
                                       const table_GF &table,
                                       const float SNR,
                                       const vector<vector<uint16_t>> &bin_table,
                                       vector<decoder_t> &chan_LLR_sorted,
                                       vector<uint16_t> &KSYMB)
{
    uint16_t N = dec_param.N, K = dec_param.K, q = dec_param.q;
    uint16_t nm = dec_param.nm;
    float sigma = sqrt(1.0 / (pow(10, SNR / 10.0)));// N0/2 or N0?
    vector<uint16_t> NSYMB(N);
    vector<vector<uint16_t>> KBIN(K), NBIN;
    NBIN.resize(N, vector<uint16_t>());
    vector<vector<softdata_t>> noisy_sig(N, vector<softdata_t>(bin_table[0].size(), (softdata_t)0.0));
    vector<uint16_t> u_symb(N, 0);
    RandomBinaryGenerator(K, q, bin_table, 0, 0, KBIN, KSYMB);
    for (int i = 0; i < K; i++)
        u_symb[dec_param.reliab_sequence[i]] = KSYMB[i];
    Encoder(table.ADDDEC, table.MULDEC, dec_param.polar_coeff, u_symb, NSYMB);
    // vector<uint16_t> temp=NSYMB;
    // inv_Encoder(table.ADDDEC, table.DIVDEC, dec_param.polar_coeff, NSYMB, temp);
    for (int i = 0; i < int(NSYMB.size()); i++)
        NBIN[i] = bin_table[NSYMB[i]];
    awgn_channel_noise(NBIN, sigma, 0, 0, noisy_sig);
    vector<vector<softdata_t>> chan_LLR(N, vector<softdata_t>(q, 0));
    Channel_LLR(noisy_sig, bin_table, q, sigma, chan_LLR);
    LLR_sort(chan_LLR, nm, chan_LLR_sorted);
}

void PoAwN::channel::EncodeChanGF_CCSK(const decoder_parameters &dec_param,
                                       const table_GF &table,
                                       const float SNR,
                                       const vector<vector<uint16_t>> &CCSK_rotated_codes,
                                       vector<decoder_t> &chan_LLR_sorted,
                                       vector<uint16_t> &KSYMB)
{
    uint16_t N = dec_param.N, K = dec_param.K, q = dec_param.q, csk_sz = CCSK_rotated_codes[0].size();
    uint16_t nm = dec_param.nm;
    float sigma = sqrt(1.0 / ((pow(10, SNR / 10.0))));// N0/2 or N0?
    vector<uint16_t> NSYMB(N);
    vector<vector<uint16_t>> KBIN(K);

    vector<uint16_t> u_symb(N, 0);
    RandomBinaryGenerator(K, q, CCSK_rotated_codes, 0, 0, KBIN, KSYMB);
    for (int i = 0; i < K; i++)
        u_symb[dec_param.reliab_sequence[i]] = KSYMB[i];
    Encoder(table.ADDDEC, table.MULDEC, dec_param.polar_coeff, u_symb, NSYMB);

    softdata_t modulation_I[q];
    softdata_t modulation_Q[q];
    for (int i = 0; i < q; i++)
    {
        modulation_I[i] = sin(i * 2 * PI / 64);
        modulation_Q[i] = cos(i * 2 * PI / 64);
    }
    vector<vector<vector<softdata_t>>> NBIN(N, vector<vector<softdata_t>>(csk_sz, vector<softdata_t>(2)));
    vector<vector<vector<softdata_t>>> rnd_noise(N, vector<vector<softdata_t>>(csk_sz, vector<softdata_t>(2)));
    vector<vector<softdata_t>> noisy_sym(csk_sz, vector<softdata_t>(2));
    vector<vector<softdata_t>> chan_LLR(N, vector<softdata_t>(q, 0));

    softdata_t min1;
    for (int n = 0; n < N; n++)
    {
        min1 = std::numeric_limits<softdata_t>::max();
        AWGN_gen(0, sigma, 0, 0, rnd_noise[n]);
        for (int i = 0; i < csk_sz; i++)
        {
            NBIN[n][i][0] = modulation_I[CCSK_rotated_codes[NSYMB[n]][i]];
            NBIN[n][i][1] = modulation_Q[CCSK_rotated_codes[NSYMB[n]][i]];
            noisy_sym[i][0] = NBIN[n][i][0] + rnd_noise[n][i][0];
            noisy_sym[i][1] = NBIN[n][i][1] + rnd_noise[n][i][1];
        }

        for (int k = 0; k < q; k++)
        {
            for (int i = 0; i < csk_sz; i++)
                chan_LLR[n][k] += pow(noisy_sym[i][0] - modulation_I[CCSK_rotated_codes[k][i]], 2) + pow(noisy_sym[i][1] - modulation_Q[CCSK_rotated_codes[k][i]], 2);
            chan_LLR[n][k] /= 2 * pow(sigma, 2);
        }
        for (int k = 0; k < q; k++)
        {
            chan_LLR[n][k] = chan_LLR[n][k] + 10 * log10(2 / pow(sigma, 2));
            if (chan_LLR[n][k] < min1)
                min1 = chan_LLR[n][k];
        }
        for (int k = 0; k < q; k++)
            chan_LLR[n][k] -= min1;

        LLR_sort(chan_LLR, nm, chan_LLR_sorted);
    }
}
