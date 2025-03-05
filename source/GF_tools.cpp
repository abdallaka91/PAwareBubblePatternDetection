#include <cmath>
#include <vector>
#include "GF_tools.h"

using std::vector;

void PoAwN::GFtools::GF_bin_seq_gen(uint16_t q, uint16_t primitivePoly,
                                    vector<vector<uint16_t>> &GF_binSeq, vector<vector<uint16_t>> &Dec_binSeq)
{
    uint16_t p = static_cast<uint16_t>(log2(q));
    Dec_binSeq.resize(q, vector<uint16_t>(p, 0));
    vector<uint16_t> fieldElements(q,100);
    fieldElements[0] = 0;
    uint16_t current = 1;
    for (uint16_t i = 1; i < q ; ++i)
    {
        fieldElements[i] = current;
        current <<= 1;
        if (current & (1 << p))
            current ^= primitivePoly;
    }

    vector<uint16_t> binRep(p, 0);
    GF_binSeq.resize(q);
    GF_binSeq[0] = binRep;
    for (int j = 1; j<q;j++)
    {
        uint16_t elem= fieldElements[j];
        vector<uint16_t> binRep(p, 0);
        for (int i = p - 1; i >= 0; --i)
        {
            binRep[(uint16_t)i] = elem & 1;
            elem >>= 1;
        }
        GF_binSeq[j] = binRep;
    }
    for (int i = 0; i < q; i++)
    {
        Dec_binSeq[fieldElements[i]] = GF_binSeq[i];
    }
}

void PoAwN::GFtools::GF_bin2GF(const vector<vector<uint16_t>> &GF_bin, vector<uint16_t> &GF_dec, vector<uint16_t> &dec_GF)
{
    uint16_t p = GF_bin[0].size();
    uint16_t q = GF_bin.size();
    dec_GF.resize(q, 5000);
    GF_dec.resize(q, 5000);
    for (uint16_t i = 0; i < q; i++)
    {
        GF_dec[i]=0;
        for (uint16_t k = 0; k < p; k++)
        {
            GF_dec[i] += GF_bin[i][k] * pow(2, p - k - 1);
        }
        dec_GF[GF_dec[i]] = i;
    }
}

void PoAwN::GFtools::GF_add_mat_gen(const vector<uint16_t> &GF_dec, vector<vector<uint16_t>> &GF_add_mat, vector<vector<uint16_t>> &DEC_add_mat)
{
    uint16_t q = GF_dec.size();
    DEC_add_mat.resize(q, std::vector<uint16_t>(q, 65535));
    GF_add_mat.resize(q);
    for (uint16_t i = 0; i < q; i++)
    {
        vector<uint16_t> Raw(q, 0);
        for (uint16_t j = 0; j < q; j++)
        {
            Raw[j] = GF_dec[i] ^ GF_dec[j];
        }
        GF_add_mat[i] = Raw;
        for (int k = 0; k < q; k++)
            DEC_add_mat[GF_dec[i]][GF_dec[k]] = Raw[k];
    }
}

void PoAwN::GFtools::GF_mul_mat_gen(const vector<uint16_t> &GFDEC, vector<vector<uint16_t>> &GF_mul_mat, vector<vector<uint16_t>> &DEC_mul_mat)
{
    uint16_t q = GFDEC.size();
    DEC_mul_mat.resize(q, std::vector<uint16_t>(q, 65535));
    // vector<vector<uint16_t>> debug_mul;
    // debug_mul.resize(q, std::vector<uint16_t>(q, 65535));
    GF_mul_mat.resize(q);
    for (uint16_t i = 0; i < q; i++)
    {
        vector<uint16_t> Raw(q, 0);
        for (uint16_t j = 0; j < q; j++)
        {
            if (i * j == 0)
                Raw[j] = 0;
            else
            {
                uint16_t i1 = i - 1;
                uint16_t j1 = j - 1;
                uint16_t pw1 = i1 + j1;
                uint16_t i2 = pw1 + 1;
                uint16_t i3 = i2;
                if (i3 >= q)
                    i3 = i2 % (q - 1);
                Raw[j] = GFDEC[i3];
            }
        }
        GF_mul_mat[i] = Raw;
        for (int k = 0; k < q; k++)
            DEC_mul_mat[GFDEC[i]][GFDEC[k]] = Raw[k];
    }
}

void PoAwN::GFtools::GF_div_mat_gen(const vector<uint16_t> &GFDEC, vector<vector<uint16_t>> &GF_div_mat, vector<vector<uint16_t>> &DEC_div_mat)
{
    uint16_t q = GFDEC.size();
    DEC_div_mat.resize(q, std::vector<uint16_t>(q, 65535));
    GF_div_mat.resize(q);
    for (uint16_t i = 0; i < q; i++)
    {
        vector<uint16_t> Raw(q, 0);
        for (uint16_t j = 0; j < q; j++)
        {
            if (i == 0)
                Raw[j] = 0;
            else if (j == 0)
                Raw[j] = 65535;
            else
            {
                int i1 = (int)i - 1;
                int j1 = (int)j - 1;
                int pw1 = i1 - j1;
                int pw = pw1;
                if (pw1 < 0)
                    pw = q - 2 + pw1 + 1;
                int i2 = pw + 1;
                Raw[j] = GFDEC[(uint16_t)i2];
            }
        }
        GF_div_mat[i] = Raw;
        for (int k = 0; k < q; k++)
            DEC_div_mat[GFDEC[i]][GFDEC[k]] = Raw[k];
    }
}
