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
#include <iomanip>
#include <algorithm>

using namespace PoAwN::structures;
using namespace PoAwN::tools;
using namespace PoAwN::init;
using namespace PoAwN::decoding;
using namespace PoAwN::channel;
using std::array;
using std::cout;
using std::endl;
using std::stod;
using std::stoi;
using std::string;
using std::vector;

int main(int argc, char *argv[])
{

    if (argc != 13)
    {
        cout << "validate: NbMonteCarlo, SNR, sig_mod(BPSK, CCSK_bin, CCSK_NB), q, N, K, nL, nH, nm, Zc, offset" << std::endl;
        return 1;
    }
    uint16_t q, N, K, n, nL, nH, nm, nb, Zc, nopM, p, frozen_val = 0;
    softdata_t offset;
    int NbMonteCarlo = stoi(argv[1]);
    float Pt, EbN0 = stod(argv[2]);
    string sig_mod = argv[3];
    std::transform(sig_mod.begin(), sig_mod.end(), sig_mod.begin(), ::toupper);
    q = stoi(argv[4]);
    p = log2(q);
    N = stoi(argv[5]);
    K = stoi(argv[6]);
    n = log2(N);
    nL = stoi(argv[7]);
    nH = stoi(argv[8]);
    nm = stoi(argv[9]);
    Zc = stoi(argv[10]);
    offset = stod(argv[11]);
    Pt = stod(argv[12]);
    nopM = nm + 4;
    nb = 4;

    base_code_t code_param(N, K, n, q, p, frozen_val);
    code_param.sig_mod = sig_mod;

    int gf_rand_SEED = 0;
    float nse_rand_SEED = 1.2544;
    bool repeatable_randgen = 0;

    table_GF table;

    cout << "Loading code_param..." << endl;
    LoadCode(code_param, EbN0);

    cout << "OK!, " << "Loading tables..." << endl;
    // void LoadTables(base_code_t & code, table_GF & table,  const uint16_t *GF_polynom_primitive)
    LoadTables(code_param, table, GF_polynom_primitive.data());
    cout << "Done!" << endl;
    cout << "Simulation starts..." << endl;

    decoder_parameters dec_param(code_param, offset, nm, nL, nH, nb, Zc, nopM);
    LoadBubblesIndcatorlists(dec_param, EbN0, Pt);

    vector<vector<std::array<int, 2>>> ns0(n);
    vector<vector<uint16_t>> ns(n);


    for (int i = 0; i < n; i++)
    {
        ns0[i].resize(1 << i);
        ns[i].resize(1<<i);
        for (int j = 0; j < 1 << i; j++)
        {
            ns0[i][j].fill(0);
            for (int k = 0; k < dec_param.Bubb_Indicator[i][j][0].size(); k++)
            {
                if (dec_param.Bubb_Indicator[i][j][0][k]+1 > ns0[i][j][0])
                    ns0[i][j][0] = dec_param.Bubb_Indicator[i][j][0][k] + 1;
                if (dec_param.Bubb_Indicator[i][j][1][k]+1 > ns0[i][j][1])
                    ns0[i][j][1] = dec_param.Bubb_Indicator[i][j][1][k] + 1;
            }
            ns[i][j]=std::max(ns0[i][j][0], ns0[i][j][1]);
        }
    }

    

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
    frozen_lay_pos(dec_param, dec_param.ucap);

    CCSK_seq ccsk_seq;
    vector<vector<uint16_t>> CCSK_rotated_codes(q, vector<uint16_t>());
    if (code_param.sig_mod == "CCSK_BIN")
        create_ccsk_rotated_table(ccsk_seq.CCSK_bin_seq[code_param.p - 1], ccsk_seq.CCSK_bin_seq[code_param.p - 1].size(), CCSK_rotated_codes);
    else if (code_param.sig_mod == "CCSK_NB")
        create_ccsk_rotated_table(ccsk_seq.CCSK_GF_seq[code_param.p - 1], ccsk_seq.CCSK_GF_seq[code_param.p - 1].size(), CCSK_rotated_codes);

    decoder_t temp_dec;
    temp_dec.intrinsic_LLR.reserve(dec_param.nm);
    temp_dec.intrinsic_GF.reserve(dec_param.nm);

    vector<vector<decoder_t>> L;
    vector<uint16_t> info_sec_rec(K, dec_param.MxUS);
    unsigned int FER = 0;
    vector<uint16_t> KSYMB(K);
    for (int i0 = 1; i0 <= NbMonteCarlo; ++i0)
    {
        L.assign(code_param.n + 1, vector<decoder_t>(code_param.N, temp_dec));

        vector<vector<uint16_t>> KBIN;
        if (code_param.sig_mod == "CCSK_BIN")
            EncodeChanBinCCSK(dec_param, table, EbN0, CCSK_rotated_codes, L[0], KSYMB);
        else if (code_param.sig_mod == "CCSK_NB")
            EncodeChanGF_CCSK(dec_param, table, EbN0, CCSK_rotated_codes, L[0], KSYMB);

        decode_SC(dec_param, table.ADDDEC, table.MULDEC, table.DIVDEC, L, info_sec_rec);

        for (uint16_t i = 0; i < dec_param.K; i++)
        {
            if (KSYMB[i] != info_sec_rec[i])
            {
                FER++;
                break;
            }
        }
        if ((i0 % 200 == 0 && i0 > 0))
            cout << "\rSNR: " << EbN0 << " dB, FER = " << FER << "/" << (float)i0 << " = " << (float)FER / (float)i0 << std::flush;
    }
    cout << endl;
}