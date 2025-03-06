
#include "Decoder_functions.h"
#include "struct.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <vector>
#include <iostream>
#include <queue>

void PoAwN::decoding::Channel_LLR(const vector<vector<softdata_t>> &chan_observ,
                                  const vector<vector<uint16_t>> &bin_symb_seq,
                                  uint16_t q,
                                  softdata_t sigma,
                                  vector<vector<softdata_t>> &chan_LLR)
{
    const int N = chan_observ.size();
    const int q1 = chan_observ[0].size();

    constexpr const softdata_t two_pow16 = softdata_t(1 << 16);

    vector<softdata_t> hard_decison(q1, two_pow16 - 1);

    // softdata_t mn_llr = two_pow16 - 1; // Never used outside of the loop
    const softdata_t fct = softdata_t(2) / (sigma * sigma);
    for (int i = 0; i < N; i++)
    {
        softdata_t mn_llr = std::numeric_limits<softdata_t>::max();
        for (int j = 0; j < q; j++)
        {
            softdata_t temp = 0;

            for (int k = 0; k < q1; k++)
            {
                // true == 1 and false == 0
                hard_decison[k] = softdata_t(uint8_t(chan_observ[i][k] <= 0));
                temp += chan_observ[i][k] * ((softdata_t)bin_symb_seq[j][k] - hard_decison[k]);
            }
            temp *= fct;
            chan_LLR[i][j] = temp;
            if (temp < mn_llr)
                mn_llr = temp;
        }
        if (mn_llr > std::numeric_limits<softdata_t>::min())
            for (int j = 0; j < q; j++)
            {
                chan_LLR[i][j] -= mn_llr;
            }
    }
}

void PoAwN::decoding::tuples_sorter(const decoder_t &tuples,
                                    const uint16_t nm,
                                    const uint16_t q,
                                    decoder_t &sorted_tuples)
{

    sorted_tuples.intrinsic_LLR.reserve(nm);
    sorted_tuples.intrinsic_GF.reserve(nm);
    size_t N1 = tuples.intrinsic_GF.size();
    vector<softdata_t> temp_llr(q, std::numeric_limits<softdata_t>::max());
    vector<uint16_t> temp_GF(q);
    iota(temp_GF.begin(), temp_GF.end(), 0);

    for (size_t i = 0; i < N1; i++)
        if (tuples.intrinsic_LLR[i] < temp_llr[tuples.intrinsic_GF[i]])
            temp_llr[tuples.intrinsic_GF[i]] = tuples.intrinsic_LLR[i];

    vector<int> indices(temp_llr.size());
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&](int i, int j)
         { return temp_llr[i] < temp_llr[j]; });

    size_t count = 0;

    for (auto idx = indices.begin(); idx < indices.end(); idx++)
    {
        sorted_tuples.intrinsic_GF.push_back(temp_GF[*idx]);
        sorted_tuples.intrinsic_LLR.push_back(temp_llr[*idx]);
        if (++count == nm)
            break;
    }
}

void PoAwN::decoding::ECN_AEMS_bubble(const decoder_t &theta_1,
                                      const decoder_t &phi_1,
                                      const decoder_parameters &dec_param,
                                      uint16_t coef,
                                      const vector<vector<uint16_t>> &ADDDEC,
                                      const vector<vector<uint16_t>> &DIVDEC,
                                      decoder_t &theta)
{
    theta.intrinsic_LLR.reserve(dec_param.nm);
    theta.intrinsic_GF.reserve(dec_param.nm);
    bool rel_theta = (theta_1.intrinsic_LLR[dec_param.Zc] > phi_1.intrinsic_LLR[dec_param.Zc]);
    decoder_t phi_1_p = phi_1;
    for (size_t i = 0; i < phi_1.intrinsic_GF.size(); i++)
        phi_1_p.intrinsic_GF[i] = DIVDEC[phi_1_p.intrinsic_GF[i]][coef];

    vector<softdata_t> srtr(dec_param.nb);
    vector<uint16_t> srtrg(dec_param.nb);
    vector<vector<uint16_t>> sij(2, vector<uint16_t>(dec_param.nb));
    vector<uint16_t> a_gf = rel_theta ? theta_1.intrinsic_GF : phi_1_p.intrinsic_GF;
    vector<softdata_t> a = rel_theta ? theta_1.intrinsic_LLR : phi_1_p.intrinsic_LLR;
    vector<uint16_t> b_gf = rel_theta ? phi_1_p.intrinsic_GF : theta_1.intrinsic_GF;
    vector<softdata_t> b = rel_theta ? phi_1_p.intrinsic_LLR : theta_1.intrinsic_LLR;

    uint16_t nH = a_gf.size() < dec_param.nH ? a_gf.size() : dec_param.nH;
    uint16_t nL = b_gf.size() < dec_param.nL ? b_gf.size() : dec_param.nL;
    vector<vector<bool>> Ti(nH, vector<bool>(nL, false));

    for (int i = 0; i < dec_param.nb; ++i)
    {
        srtr[i] = a[i];
        srtrg[i] = ADDDEC[a_gf[i]][b_gf[0]];
        sij[0][i] = i;
        sij[1][i] = 0;
        Ti[i][0] = true;
    }

    int cnt = 0, nop = 0;
    while (nop < dec_param.nopM)
    {
        ++nop;
        auto min_it = min_element(srtr.begin(), srtr.end());
        int n = distance(srtr.begin(), min_it);
        int i = sij[0][n];
        int j = sij[1][n];

        if (find(theta.intrinsic_GF.begin(), theta.intrinsic_GF.end(), srtrg[n]) == theta.intrinsic_GF.end())
        {
            theta.intrinsic_LLR.push_back(*min_it);
            theta.intrinsic_GF.push_back(srtrg[n]);
            ++cnt;
        }
        if (i == nH - 1 || j == nL - 1 || cnt == dec_param.nm)
            break;
        int H = (i == 0) ? 1 : (i == dec_param.nb - 1) ? 0
                                                       : 1;
        int Hb = 1 - H;
        int i1 = (!Ti[i + Hb][j + H]) ? i + Hb : i + H;
        int j1 = (!Ti[i + Hb][j + H]) ? j + H : j + Hb;
        Ti[i1][j1] = true;
        srtr[n] = a[i1] + b[j1];
        srtrg[n] = ADDDEC[a_gf[i1]][b_gf[j1]];
        sij[0][n] = i1;
        sij[1][n] = j1;
    }
}

void PoAwN::decoding::ECN_EMS(const decoder_t &theta_1,
                              const decoder_t &phi_1,
                              const vector<vector<uint16_t>> &ADDDEC,
                              const vector<vector<uint16_t>> &DIVDEC,
                              const uint16_t nm,
                              uint16_t q,
                              const uint16_t coef,
                              decoder_t &theta)
{
    decoder_t phi_1_p = phi_1;
    for (size_t i = 0; i < phi_1.intrinsic_GF.size(); i++)
        phi_1_p.intrinsic_GF[i] = DIVDEC[phi_1_p.intrinsic_GF[i]][coef];

    size_t theta_1_n = theta_1.intrinsic_LLR.size();
    size_t Phi_1_n = theta_1.intrinsic_LLR.size();
    decoder_t bubble;
    bubble.intrinsic_GF.resize(theta_1_n * Phi_1_n);
    bubble.intrinsic_LLR.resize(theta_1_n * Phi_1_n);

    for (size_t i = 0; i < theta_1_n; i++)
    {
        for (size_t j = 0; j < Phi_1_n; ++j)
        {
            bubble.intrinsic_GF[i * Phi_1_n + j] = ADDDEC[theta_1.intrinsic_GF[i]][phi_1_p.intrinsic_GF[j]];
            bubble.intrinsic_LLR[i * Phi_1_n + j] = theta_1.intrinsic_LLR[i] + phi_1_p.intrinsic_LLR[j];
        }
    }
    tuples_sorter(bubble, nm, q, theta);
}

void PoAwN::decoding::ECN_EMS_L(const decoder_t &theta_1,
                                const decoder_t &phi_1,
                                const vector<vector<uint16_t>> &ADDDEC,
                                const vector<vector<uint16_t>> &DIVDEC,
                                const decoder_parameters &dec_param,
                                const uint16_t coef,
                                decoder_t &theta,
                                vector<vector<uint16_t>> &Bt1)
{
    decoder_t phi_1_p = phi_1;
    bool rel_theta = (theta_1.intrinsic_LLR[dec_param.Zc] > phi_1.intrinsic_LLR[dec_param.Zc]);
    for (size_t i = 0; i < phi_1.intrinsic_GF.size(); i++)
        phi_1_p.intrinsic_GF[i] = DIVDEC[phi_1_p.intrinsic_GF[i]][coef];

    vector<uint16_t> a_gf = rel_theta ? theta_1.intrinsic_GF : phi_1_p.intrinsic_GF;
    vector<softdata_t> a = rel_theta ? theta_1.intrinsic_LLR : phi_1_p.intrinsic_LLR;
    vector<uint16_t> b_gf = rel_theta ? phi_1_p.intrinsic_GF : theta_1.intrinsic_GF;
    vector<softdata_t> b = rel_theta ? phi_1_p.intrinsic_LLR : theta_1.intrinsic_LLR;

    uint16_t nH = a_gf.size() < dec_param.nH ? a_gf.size() : dec_param.nH;
    uint16_t nL = b_gf.size() < dec_param.nL ? b_gf.size() : dec_param.nL;

    decoder_t tuples;
    tuples.intrinsic_GF.resize(nH * nL);
    tuples.intrinsic_LLR.resize(nH * nL);
    vector<vector<uint16_t>> idxs(nH * nL, vector<uint16_t>(2, 0));

    for (size_t i = 0; i < nH; i++)
        for (size_t j = 0; j < nL; ++j)
        {
            idxs[i * nL + j][0] = i;
            idxs[i * nL + j][1] = j;
            tuples.intrinsic_GF[i * nL + j] = ADDDEC[a_gf[i]][b_gf[j]];
            tuples.intrinsic_LLR[i * nL + j] = a[i] + b[j];
        }

    uint16_t nm = dec_param.nm;
    uint16_t q = dec_param.q;

    theta.intrinsic_LLR.reserve(nm);
    theta.intrinsic_GF.reserve(nm);
    size_t N1 = tuples.intrinsic_GF.size();
    vector<softdata_t> temp_llr(q, std::numeric_limits<softdata_t>::max());
    vector<uint16_t> temp_GF(q);
    iota(temp_GF.begin(), temp_GF.end(), 0);

    vector<vector<uint16_t>> sort_bub_idxs;
    sort_bub_idxs.reserve(nm);
    vector<vector<uint16_t>> temp_idxs(q, vector<uint16_t>(2, 0));
    ;

    for (size_t i = 0; i < N1; i++)
        if (tuples.intrinsic_LLR[i] < temp_llr[tuples.intrinsic_GF[i]])
        {
            temp_llr[tuples.intrinsic_GF[i]] = tuples.intrinsic_LLR[i];
            temp_idxs[tuples.intrinsic_GF[i]] = idxs[i];
        }

    vector<int> indices(temp_llr.size());
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&](int i, int j)
         { return temp_llr[i] < temp_llr[j]; });

    size_t count = 0;

    for (auto idx = indices.begin(); idx < indices.end(); idx++)
    {
        theta.intrinsic_GF.push_back(temp_GF[*idx]);
        theta.intrinsic_LLR.push_back(temp_llr[*idx]);
        sort_bub_idxs.push_back({temp_idxs[*idx][0], temp_idxs[*idx][1]});
        if (++count == nm)
            break;
    }
    for (int i = 0; i < dec_param.nm; i++)
        Bt1[sort_bub_idxs[i][0]][sort_bub_idxs[i][1]] += 1;
}

void PoAwN::decoding::LLR_sort(const vector<vector<softdata_t>> &chan_LLR,
                               const uint16_t nm,
                               vector<decoder_t> &chan_LLR_sorted)
{
    int N = chan_LLR.size();
    int q = chan_LLR[0].size();
    chan_LLR_sorted = vector<decoder_t>(N);
    vector<size_t> indices(q);
    iota(indices.begin(), indices.end(), 0);
    vector<softdata_t> temp_llr_arr(q);
    vector<uint16_t> temp_gf_arr(q);
    for (int k = 0; k < N; ++k)
    {
#include <limits> // Make sure this is included

        chan_LLR_sorted[k].intrinsic_LLR.resize(nm,
                                                std::numeric_limits<typename decltype(chan_LLR_sorted[k].intrinsic_LLR)::value_type>::max());
        chan_LLR_sorted[k].intrinsic_GF.resize(nm,
                                               std::numeric_limits<typename decltype(chan_LLR_sorted[k].intrinsic_GF)::value_type>::max());
        for (int i = 0; i < q; ++i)
            temp_llr_arr[i] = chan_LLR[k][i];
        partial_sort(indices.begin(), indices.begin() + nm, indices.end(),
                     [&temp_llr_arr](int i1, int i2)
                     { return temp_llr_arr[i1] < temp_llr_arr[i2]; });
        for (int i = 0; i < nm; ++i)
        {
            chan_LLR_sorted[k].intrinsic_LLR[i] = temp_llr_arr[indices[i]];
            chan_LLR_sorted[k].intrinsic_GF[i] = indices[i];
        }
    }
}

void PoAwN::decoding::VN_update(const decoder_t &theta_1,
                                const decoder_t &phi_1,
                                const vector<vector<uint16_t>> &ADDDEC,
                                const vector<vector<uint16_t>> &DIVDEC,
                                uint16_t nm,
                                uint16_t coef,
                                uint16_t hard_decision,
                                uint16_t q,
                                softdata_t ofst,
                                decoder_t &phi)
{
    vector<softdata_t> VN_in2_llr;
    vector<softdata_t> VN_in1_llr;

    phi.intrinsic_GF.reserve(nm);
    phi.intrinsic_LLR.reserve(nm);

    auto mx_theta_llr = theta_1.intrinsic_LLR.back();
    auto mx_phi_llr = phi_1.intrinsic_LLR.back();

    VN_in2_llr.resize(q, (softdata_t)mx_theta_llr + ofst);
    VN_in1_llr.resize(q, (softdata_t)mx_phi_llr + ofst);

    uint16_t temp_gf1;
    for (uint16_t i = 0; i < theta_1.intrinsic_LLR.size(); i++)
    {
        temp_gf1 = ADDDEC[hard_decision][theta_1.intrinsic_GF[i]];
        VN_in2_llr[temp_gf1] = theta_1.intrinsic_LLR[i];
    }

    uint16_t temp_gf2;
    softdata_t mn_tmp = std::numeric_limits<softdata_t>::max();
    for (uint16_t i = 0; i < phi_1.intrinsic_LLR.size(); i++)
    {
        temp_gf2 = DIVDEC[phi_1.intrinsic_GF[i]][coef];
        VN_in1_llr[temp_gf2] = phi_1.intrinsic_LLR[i];
    }

    vector<softdata_t> temp_LLR(q);
    for (uint16_t i = 0; i < q; i++)
    {
        temp_LLR[i] = VN_in1_llr[i] + VN_in2_llr[i];
        if (temp_LLR[i] < mn_tmp)
            mn_tmp = temp_LLR[i];
    }
    if (mn_tmp > 1e-5)
        for (uint16_t i = 0; i < q; i++)
            temp_LLR[i] -= mn_tmp;

    std::priority_queue<std::pair<softdata_t, uint16_t>> maxHeap;

    for (uint16_t i = 0; i < temp_LLR.size(); ++i)
    {
        maxHeap.push({temp_LLR[i], i});
        if (maxHeap.size() > nm)
            maxHeap.pop();
    }
    vector<std::pair<softdata_t, uint16_t>> kSmallest;
    kSmallest.reserve(nm + 1);
    while (!maxHeap.empty())
    {
        kSmallest.push_back(maxHeap.top());
        maxHeap.pop();
    }
    for (int i = nm - 1; i >= 0; i--)
    {
        phi.intrinsic_GF.push_back((uint16_t)kSmallest[i].second);
        phi.intrinsic_LLR.push_back((softdata_t)kSmallest[i].first);
    }
}

void PoAwN::decoding::decode_SC(const decoder_parameters &dec_param,
                                const vector<vector<uint16_t>> &ADDDEC,
                                const vector<vector<uint16_t>> &DIVDEC,
                                const vector<vector<uint16_t>> &MULDEC,
                                vector<vector<decoder_t>> &L,
                                vector<uint16_t> &info_sec_rec,
                                vector<vector<vector<vector<uint16_t>>>> &Bt)
{
    vector<vector<bool>> Roots = dec_param.Roots_V;
    uint16_t MxUS = dec_param.MxUS, n = dec_param.n, N = dec_param.N;
    int l = 0, s = 0;
    uint16_t hard_decsion, temp_coef, i1, i2, i3, SZc, SZc1, l1;
    vector<uint16_t> Root;

    vector<vector<uint16_t>> V(n + 1, vector<uint16_t>(dec_param.N, dec_param.MxUS));
    V[n] = (vector<uint16_t>(dec_param.N, dec_param.frozen_val));

    for (uint16_t i = 0; i < dec_param.K; i++)
        V[n][dec_param.reliab_sequence[i]] = dec_param.MxUS;

    while (l > -1)
    {
        if (Roots[l][s])
        {
            if (s % 2 == 1)
            {
                l -= 1;
                s = (s - 1) / 2;
                Root = dec_param.Roots_indices[l][s];
                SZc = Root.size();
                SZc1 = SZc >> 1;

                for (uint16_t t = 0; t < SZc1; t++)
                {
                    l1 = l + 1;
                    i1 = Root[t], i2 = Root[t + SZc1], i3 = dec_param.coefs_id[l][s][t];
                    temp_coef = dec_param.polar_coeff[l][i3];
                    V[l][i1] = ADDDEC[V[l1][i1]][V[l1][i2]];
                    V[l][i2] = MULDEC[V[l1][i2]][temp_coef];
                }
                Roots[l][s] = true;
            }
            else
            {
                l -= 1;
                s /= 2;
            }
        }
        else if (Roots[l + 1][2 * s])
        {
            Root = dec_param.Roots_indices[l][s];
            SZc = Root.size();
            SZc1 = SZc >> 1;
            decoder_t theta_1, phi_1;
            for (uint16_t t = 0; t < SZc1; t++)
            {
                theta_1 = L[l][Root[t]];
                phi_1 = L[l][Root[t + SZc1]];
                i3 = dec_param.coefs_id[l][s][t];
                temp_coef = dec_param.polar_coeff[l][i3];
                decoder_t phi;
                hard_decsion = V[l + 1][Root[t]];
                VN_update(theta_1, phi_1, ADDDEC, DIVDEC, dec_param.nm,
                          temp_coef, hard_decsion, dec_param.q, dec_param.offset, phi);
                L[l + 1][Root[t + SZc1]] = phi;
            }
            l += 1;
            s = 2 * s + 1;
            if (l == n)
            {
                Roots[n][s] = true;
                if (V[n][s] == MxUS)
                    V[n][s] = L[n][s].intrinsic_GF[0];
                if (s == N - 1)
                    break;
            }
        }
        else
        {
            Root = dec_param.Roots_indices[l][s];
            SZc = Root.size();
            SZc1 = SZc >> 1;
            decoder_t phi_1, theta_1;
            for (uint16_t t = 0; t < SZc1; t++)
            {
                i3 = dec_param.coefs_id[l][s][t];
                theta_1 = L[l][Root[t]];
                phi_1 = L[l][Root[t + SZc1]];
                temp_coef = dec_param.polar_coeff[l][i3];
                decoder_t theta;
                // ECN_EMS(theta_1, phi_1, ADDDEC, DIVDEC, dec_param.nm, dec_param.q, temp_coef, theta);
                ECN_EMS_L(theta_1, phi_1, ADDDEC, DIVDEC, dec_param, temp_coef, theta, Bt[l][s]);
                // ECN_AEMS_bubble(theta_1, phi_1, dec_param, temp_coef, ADDDEC, DIVDEC, theta);
                L[l + 1][Root[t]] = theta;
            }
            l = l + 1;
            s = 2 * s;
            if (l == n)
            {
                Roots[n][s] = true;
                if (V[n][s] == MxUS)
                    V[n][s] = L[n][s].intrinsic_GF[0];
            }
        }
    }
    info_sec_rec.resize(dec_param.K, dec_param.MxUS);
    for (uint16_t i = 0; i < dec_param.K; i++)
    {
        info_sec_rec[i] = V[n][dec_param.reliab_sequence[i]];
    }
}
