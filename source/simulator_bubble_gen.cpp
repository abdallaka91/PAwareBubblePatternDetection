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
#include "HelperFunc.h"
#include <fstream>
#include <iomanip>
#include <cstring>
#include <string>
#include <iomanip>
#include <algorithm>
#include "channel.h"

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
        cout << "validate: NbMonteCarlo, SNR, sig_mod(BPSK, CCSK_bin, CCSK_NB), q, N, K, nL, nH, nm, Zc, offset, Pt" << std::endl;
        return 1;
    }
    uint16_t q, N, K, n, nL, nH, nm, nb, Zc, nopM, p, frozen_val = 0;
    softdata_t offset;
    int NbMonteCarlo = stoi(argv[1]);
    float Pt, EbN0 = stod(argv[2]);
    string sig_mod = argv[3];
    ;
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

    float sigma = sqrt(1.0 / (pow(10, EbN0 / 10.0)));

    dec_param.Roots_V.resize(n + 1);
    dec_param.Roots_indices.resize(n);
    dec_param.clusts_CNs.resize(n);
    dec_param.clusts_VNs.resize(n);
    dec_param.coefs_id.resize(n);

    dec_param.Roots_V[n].resize(1U << n, false);

    vector<vector<vector<vector<uint16_t>>>> Bt;
    vector<vector<vector<vector<float>>>> Cs;
    vector<vector<vector<vector<uint64_t>>>> Rs;
    Cs.resize(n);
    Bt.resize(n);
    Rs.resize(n);

    for (uint16_t l = 0; l < n; l++)
    {
        dec_param.Roots_V[l].resize(1U << l, false);
        dec_param.Roots_indices[l].resize(pow(2, l));
        dec_param.clusts_CNs[l].resize(pow(2, l));
        dec_param.clusts_VNs[l].resize(pow(2, l));
        dec_param.coefs_id[l].resize(pow(2, l));
        Cs[l].resize(pow(2, l));
        Bt[l].resize(pow(2, l));
        for (uint16_t s = 0; s < dec_param.Roots_V[l].size(); s++)
        {
            Cs[l][s].assign(nH, vector<float>(nL, 0));
            Bt[l][s].assign(nH, vector<uint16_t>(nL, 0));
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

    CCSK_seq ccsk_seq;
    vector<vector<uint16_t>> CCSK_rotated_codes(q, vector<uint16_t>());
    if (code_param.sig_mod == "CCSK_BIN")
        create_ccsk_rotated_table(ccsk_seq.CCSK_bin_seq[code_param.p - 1], ccsk_seq.CCSK_bin_seq[code_param.p - 1].size(), CCSK_rotated_codes);
    else if (code_param.sig_mod == "CCSK_NB")
        create_ccsk_rotated_table(ccsk_seq.CCSK_GF_seq[code_param.p - 1], ccsk_seq.CCSK_GF_seq[code_param.p - 1].size(), CCSK_rotated_codes);

    decoder_t temp_dec;
    temp_dec.intrinsic_LLR.resize(dec_param.nm, 0);
    temp_dec.intrinsic_GF.resize(dec_param.nm, 0);

    vector<vector<decoder_t>> L;
    vector<uint16_t> info_sec_rec(K, dec_param.MxUS);
    unsigned int FER = 0;
    vector<uint16_t> KSYMB(K);
    bool succ_dec, succ_writing, newsim;
    int succ_dec_frame = 0, i0 = 0;
    while (succ_dec_frame < NbMonteCarlo)
    {
        i0++;
        succ_dec = 1;
        L.assign(code_param.n + 1, vector<decoder_t>(code_param.N, temp_dec));

        vector<vector<uint16_t>> KBIN;
        if (code_param.sig_mod == "CCSK_BIN")
            EncodeChanBPSK_BinCCSK(dec_param, table, EbN0, CCSK_rotated_codes, L[0], KSYMB);
        else if (code_param.sig_mod == "CCSK_NB")
            EncodeChanGF_CCSK(dec_param, table, EbN0, CCSK_rotated_codes, L[0], KSYMB);
        else
            EncodeChanBPSK_BinCCSK(dec_param, table, EbN0, table.BINDEC, L[0], KSYMB);


            decode_SC_bubble_gen(dec_param, table.ADDDEC, table.MULDEC, table.DIVDEC, L, info_sec_rec, Bt);

        for (uint16_t i = 0; i < dec_param.K; i++)
        {
            if (KSYMB[i] != info_sec_rec[i])
            {
                succ_dec = false;
                FER++;
                break;
            }
        }
        if (succ_dec)
        {
            succ_dec_frame++;
            for (uint16_t l = 0; l < n; l++)
                for (uint16_t s = 0; s < N >> (n - l); s++)
                    for (int j0 = 0; j0 < nH; j0++)
                        for (int j1 = 0; j1 < nL; j1++)
                            Cs[l][s][j0][j1] += (float)Bt[l][s][j0][j1];
        }
        for (uint16_t l = 0; l < n; l++)
            for (uint16_t s = 0; s < N >> (n - l); s++)
                for (auto &rw : Bt[l][s])
                    for (auto &elem : rw)
                        elem = 0;

        if ((i0 % 100 == 0 && i0 > 0))
            cout << "\rSNR: " << EbN0 << " dB, FER = " << FER << "/" << (float)i0 << " = " << (float)FER / (float)i0 << std::flush;
    }

    for (uint16_t l = 0; l < n; l++)
        for (uint16_t s = 0; s < N >> (n - l); s++)
            for (int j0 = 0; j0 < nH; j0++)
                for (int j1 = 0; j1 < nL; j1++)
                {
                    Cs[l][s][j0][j1] /= (float)(1 << (n - (l + 1))); // normalize 2^(n-l+1)
                    Cs[l][s][j0][j1] /= (float)NbMonteCarlo;
                    if (Cs[l][s][j0][j1] > Pt)
                        Bt[l][s][j0][j1] = true;
                }

    cout << "\rSNR: " << EbN0 << " dB, FER = " << FER << "/" << (float)i0 << " = " << (float)FER / (float)i0 << std::flush;
    cout << endl;
    newsim = true;
    std::ostringstream fname;
    string bubble_direct;
    if (code_param.sig_mod == "BPSK")
        bubble_direct = "./BubblesPattern/bpsk/N";
    else if (code_param.sig_mod == "CCSK_BIN")
        bubble_direct = "./BubblesPattern/ccsk_bin/N";
    else
        bubble_direct = "./BubblesPattern/ccsk_nb/N";

    // fname << "/mnt/c/Users/Abdallah Abdallah/Desktop/BubblePattern/"
    //       << "N" << code_param.N << "_GF" << code_param.q
    //       << "_SNR" << std::fixed << std::setprecision(2) << EbN0 << ".txt";
    fname << bubble_direct << code_param.N << "/CotributionMatrices/" << "bubbles_N" << code_param.N << "_GF" << code_param.q
          << "_SNR" << std::fixed << std::setprecision(2) << EbN0 << "_" << dec_param.nH << "x" << dec_param.nL << "_Pt"
          << std::fixed << std::setprecision(2) << Pt << "_Cs_mat.txt";

    std::string filename = fname.str();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 1u << i; j++)
        {
            succ_writing = appendToFile(fname.str(), Cs[i][j], newsim);
            newsim = false;
        }
    }
    if (succ_writing)
        std::cout << "Cs Matrices written to: " << filename << std::endl;

    newsim = true;
    fname.str("");
    fname.clear();

    fname << bubble_direct << code_param.N << "/CotributionMatrices/" << "bubbles_N" << code_param.N << "_GF" << code_param.q
          << "_SNR" << std::fixed << std::setprecision(2) << EbN0 << "_" << dec_param.nH << "x" << dec_param.nL << "_Pt"
          << std::fixed << std::setprecision(2) << Pt << "_Bt_mat.txt";

    filename = fname.str();

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 1u << i; j++)
        {
            succ_writing = appendToFile(fname.str(), Bt[i][j], newsim);
            newsim = false;
        }
    }
    if (succ_writing)
        std::cout << "Bt Matrices written to: " << filename << std::endl;

    newsim = true;
    fname.str("");
    fname.clear();
    fname << bubble_direct << code_param.N << "/bubbles_N" << code_param.N << "_GF" << code_param.q
          << "_SNR" << std::fixed << std::setprecision(2) << EbN0 << "_" << dec_param.nH << "x" << dec_param.nL << "_Pt"
          << std::fixed << std::setprecision(2) << Pt << "_Bt_lsts.txt";

    filename = fname.str();

    std::filesystem::create_directories(std::filesystem::path(filename).parent_path());

    FILE *file = fopen(filename.c_str(), newsim ? "w" : "a");
    if (!file)
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    else
    {
        for (int l = 0; l < n; l++)
        {
            for (int s = 0; s < (1u << l); s++)
            {
                std::ostringstream lne;
                for (int j0 = 0; j0 < nH; j0++)
                    for (int j1 = 0; j1 < nL; j1++)
                        if (Bt[l][s][j0][j1])
                            lne << j0 << "  " << j1 << std::setw(8);
                if (!(l == n - 1 && s == (1u << l) - 1))
                    lne << "\n";
                fprintf(file, "%s", lne.str().c_str());
                newsim = false;
            }
            fprintf(file, "%s", "\n");
        }
        fclose(file);
        std::cout << "Bt Matrices written to: " << filename << std::endl;
    }
}