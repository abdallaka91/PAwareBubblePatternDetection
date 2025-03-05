#ifndef STRUCT_H_INCLUDED
#define STRUCT_H_INCLUDED

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace PoAwN {
namespace structures {
    using std::vector;

struct csk_t {
    uint16_t PNsize;
    vector<uint16_t> PN;
    vector<vector<uint16_t>> CSK_arr;
};

using softdata_t = float;

struct decoder_t {

    vector<softdata_t> intrinsic_LLR;
    vector<uint16_t> intrinsic_GF;

    decoder_t(vector<softdata_t> v1 = {}, vector<uint16_t> v2 = {})
        : intrinsic_LLR(v1), intrinsic_GF(v2) {}
};

struct table_GF {
    vector<vector<uint16_t>> BINGF;  
    vector<vector<uint16_t>> BINDEC;
    vector<vector<uint16_t>> ADDGF;
    vector<vector<uint16_t>> ADDDEC;
    vector<vector<uint16_t>> MULGF;
    vector<vector<uint16_t>> DIVGF;
    vector<uint16_t> DECGF;
    vector<uint16_t> GFDEC;
    vector<vector<uint16_t>> MULDEC;
    vector<vector<uint16_t>> DIVDEC;
};

struct base_code_t {
    uint16_t N, K, n, q, p, frozen_val;
    vector<uint16_t> reliab_sequence, frozen_ind;
    float Rate;
    vector<vector<uint16_t>> polar_coeff;
    std::string sig_mod;
    uint16_t MxUS = 1 << 15;
    base_code_t(uint16_t N = 0, uint16_t K=0, uint16_t n=0, uint16_t q=0, uint16_t p=0, uint16_t frozen_val = 0)
    : N(N), K(K), n(n), q(q), p(p), frozen_val(frozen_val){}
};

struct decoder_parameters : public base_code_t {
    vector<vector<bool>> Roots_V;
    vector<vector<vector<uint16_t>>> Roots_indices, clusts_CNs, clusts_VNs, coefs_id;
    softdata_t offset;
    uint16_t nm, nL, nH, nb, Zc, nopM;
    decoder_parameters(const base_code_t& base, softdata_t offset = 0, uint16_t nm = 0, 
        uint16_t nL = 0, uint16_t nH = 0, uint16_t nb = 0, uint16_t Zc = 0, uint16_t nopM = 0)
    : base_code_t(base),  
      offset(offset), nm(nm), nL(nL), nH(nH), nb(nb), Zc(Zc), nopM(nopM) {}
};

// clang-format off
constexpr const std::array<uint16_t, 9> GF_polynom_primitive = {
    7, 11, 19, 37, 67, 137, 285, 259, 1033
};

constexpr const std::array<uint16_t, 64> CCSK_bin_seq = {
    0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1
};
// clang-format on

} // namespace structures
} // namespace PoAwN

#endif
