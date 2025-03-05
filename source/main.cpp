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

int main(int argc, char *argv[])
{

    if (argc != 13 && argc != 1)
    {
        cout << "validate: NbMonteCarlo, SNR, is_bpsk, q, N, K, nL, nH, nm, nb, Zc, offset" << std::endl;
        return 1;
    }

    int NbMonteCarlo = stoi(argv[1]);
    float EbN0 = stod(argv[2]);
    bool is_bpsk = (bool)stod(argv[3]);

    uint16_t q, N, K, n, nL, nH, nm, nb, Zc, nopM, p, frozen_val = 0;
    softdata_t offset;
    q = stoi(argv[4]);
    p = log2(q);
    N = stoi(argv[5]);
    K = stoi(argv[6]);
    n = log2(N);
    nL = stoi(argv[7]);
    nH = stoi(argv[8]);
    nm = stoi(argv[9]);
    nb = stoi(argv[10]);
    Zc = stoi(argv[11]);
    offset = stod(argv[12]);
    nopM = nm + 4;

     base_code_t code_param(N, K, n, q, p, frozen_val);
    code_param.sig_mod = is_bpsk ? "bpsk" : "ccsk";
    
    int gf_rand_SEED = 0;
    float nse_rand_SEED = 1.2544;
    bool repeatable_randgen = 0;

    table_GF table;

    cout << "Loading code_param..." << endl;
    LoadCode(code_param, EbN0);

    vector<vector<uint16_t>> CCSK_rotated_codes(q, vector<uint16_t>());
    create_ccsk_rotated_table(CCSK_bin_seq.data(), CCSK_bin_seq.size(), CCSK_rotated_codes);
    cout << "OK!, " << "Loading tables..." << endl;
    // void LoadTables(base_code_t & code, table_GF & table,  const uint16_t *GF_polynom_primitive)
    LoadTables(code_param, table, GF_polynom_primitive.data());
    cout << "Done!" << endl;
    cout << "Simulation starts..." << endl;

     decoder_parameters dec_param(code_param, offset, nm, nL, nH, nb, Zc, nopM);


    float sigma = sqrt(1.0 / (pow(10, EbN0 / 10.0)));

    dec_param.Roots_V.resize(n + 1);
    dec_param.Roots_indices.resize(n);
    dec_param.clusts_CNs.resize(n);
    dec_param.clusts_VNs.resize(n);
    dec_param.coefs_id.resize(n);

    dec_param.Roots_V[n].resize(1U << n, false);
    for (uint16_t l = 0; l < n; l++)
    {
        dec_param.Roots_V[l].resize(1U << l, false);
        dec_param.Roots_indices[l].resize(pow(2, l));
        dec_param.clusts_CNs[l].resize(pow(2, l));
        dec_param.clusts_VNs[l].resize(pow(2, l));
        dec_param.coefs_id[l].resize(pow(2, l));
        for (uint16_t s = 0; s < dec_param.Roots_V[l].size(); s++)
        {
            uint16_t sz1 = N >> (l + 1U), sz2 = sz1 << 1U;
            dec_param.clusts_CNs[l][s].resize(sz1);
            dec_param.clusts_VNs[l][s].resize(sz1);
            dec_param.coefs_id[l][s].resize(sz1);
            for (uint16_t t = 0; t < sz1; ++t)
            {
                dec_param.clusts_CNs[l][s][t] = s * sz2 + (t << 1U) - (t % sz1);
                dec_param.clusts_VNs[l][s][t] = dec_param.clusts_CNs[l][s][t] + (N >> (l + 1U));
                dec_param.coefs_id[l][s][t] = dec_param.clusts_VNs[l][s][t] - (s + 1) * sz1;
            }
        }
        for (uint16_t s = 0; s < dec_param.Roots_indices[l].size(); s++)
        {

            uint16_t sz1 = N >> l;
            dec_param.Roots_indices[l][s].resize(sz1);
            for (uint16_t t = 0; t < sz1; ++t)
                dec_param.Roots_indices[l][s][t] = s * sz1 + t;
        }
    }

    vector<uint16_t> NSYMB(N);
    vector<vector<uint16_t>> KBIN, NBIN;
    unsigned int FER = 0;
    vector<vector<decoder_t>> L;
    NBIN.resize(code_param.N, vector<uint16_t>());

    decoder_t temp_dec;

    temp_dec.intrinsic_LLR.reserve(dec_param.nm);
    temp_dec.intrinsic_GF.reserve(dec_param.nm);

    vector<vector<softdata_t>> chan_LLR(N, vector<softdata_t>(code_param.q, 0));
    vector<uint16_t> info_sec_rec(K, dec_param.MxUS);
    vector<vector<softdata_t>> noisy_sig(code_param.N,
                                         vector<softdata_t>(is_bpsk ? code_param.p : code_param.q, (softdata_t)0.0));

    vector<vector<uint16_t>> bin_table;
    bin_table = code_param.sig_mod == "bpsk" ? table.BINDEC : CCSK_rotated_codes;
    vector<uint16_t> u_symb(code_param.N, 0);

    for (int i0 = 1; i0 <= NbMonteCarlo; ++i0)
    {
        L.resize(code_param.n + 1, vector<decoder_t>(code_param.N, temp_dec));
        vector<uint16_t> KSYMB;
        vector<vector<uint16_t>> KBIN;

        RandomBinaryGenerator(code_param, bin_table,
                              repeatable_randgen, gf_rand_SEED, KBIN, KSYMB);
       
        for (int i = 0; i < code_param.K; i++)
            u_symb[code_param.reliab_sequence[i]] = KSYMB[i];
        Encoder(table.ADDDEC, table.MULDEC, code_param.polar_coeff, u_symb, NSYMB);

        if (is_bpsk)
            for (int i = 0; i < int(NSYMB.size()); i++)
                NBIN[i] = table.BINDEC[NSYMB[i]];
        else
            for (int i = 0; i < int(NSYMB.size()); i++)
                NBIN[i] = CCSK_rotated_codes[NSYMB[i]];

        awgn_channel_noise(NBIN, sigma, repeatable_randgen, nse_rand_SEED, noisy_sig);
        Channel_LLR(noisy_sig, is_bpsk ? table.BINDEC : CCSK_rotated_codes, code_param.q,
                    sigma, chan_LLR);
        LLR_sort(chan_LLR, dec_param.nm, L[0]);

        decode_SC(dec_param, table.ADDDEC, table.MULDEC, table.DIVDEC, L, info_sec_rec);

        for (uint16_t i = 0; i < dec_param.K; i++)
        {
            if (KSYMB[i] != info_sec_rec[i])
            {
                FER++;
                break;
            }
        }
        if ((i0  % 200 == 0 && i0 > 0))
            cout << "\rSNR: " << EbN0 << " dB, FER = " << FER << "/" << (float)i0 << " = " << (float)FER / (float)i0  << std::flush;
    }
    cout << endl;
}