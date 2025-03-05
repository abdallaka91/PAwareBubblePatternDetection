#ifndef TOOLS
#define TOOLS

#include "struct.h"
#include <cstdint>
#include <vector>

namespace PoAwN {
namespace tools {

using std::vector;
using namespace PoAwN::structures;

void circshift(std::vector<uint16_t> & CCSK_bin_seq, int shift);

void create_ccsk_rotated_table(const uint16_t * CCSK_bin_seq, const uint16_t q, vector<vector<uint16_t>> & ccsk_rotated_table);

void RandomBinaryGenerator(const base_code_t code_param, const vector<vector<uint16_t>> & bin_table,
                           const bool repeatable, const int SEED,
                           vector<vector<uint16_t>> & KBIN, vector<uint16_t> & KSYMB);
void AWGN_gen(double                       MEAN,
              double                       STD,
              bool                         repeatable,
              double                       SEED,
              vector<vector<softdata_t>> & noise_table);

void awgn_channel_noise(const vector<vector<uint16_t>> & NBIN,
                        const double                           sigma,
                        const bool                             repeatable,
                        const double                           SEED,
                        vector<vector<softdata_t>> &           noisy_sig);
void Encoder(const vector<vector<uint16_t>> & ADDDEC, const vector<vector<uint16_t>> & MULDEC,
             const vector<vector<uint16_t>> & polar_coeff,
             const vector<uint16_t> u_symb, vector<uint16_t> & NSYMB);
} // namespace tools
} // namespace PoAwN
#endif