#ifndef DECODER_FUNCTIONS
#define DECODER_FUNCTIONS
#include "struct.h"

#include <cstdint>
#include <vector>

namespace PoAwN
{

    namespace decoding
    {

        using structures::decoder_parameters;
        using structures::decoder_t;
        using structures::softdata_t;
        using structures::vector;

        void Channel_LLR(const vector<vector<softdata_t>> &chan_observ,
                         const vector<vector<uint16_t>> &bin_symb_seq,
                         uint16_t q,
                         softdata_t sigma,
                         vector<vector<softdata_t>> &chan_LLR);

        void LLR_sort(const vector<vector<softdata_t>> &chan_LLR,
                      const uint16_t nm,
                      vector<decoder_t> &chan_LLR_sorted);

        void ECN_EMS_bubble(const decoder_t &theta_1,
                            const decoder_t &phi_1,
                            const decoder_parameters &dec_param,
                            uint16_t coef,
                            const vector<vector<uint16_t>> &ADDDEC,
                            const vector<vector<uint16_t>> &DIVDECC,
                            decoder_t &theta);
        void ECN_PA(const decoder_t &theta_1,
                    const decoder_t &phi_1,
                    const vector<vector<uint16_t>> &ADDDEC,
                    const vector<vector<uint16_t>> &DIVDEC,
                    const decoder_parameters &dec_param,
                    const uint16_t coef,
                    const vector<vector<uint16_t>> &Bubb_Indicator,
                    decoder_t &theta1);

        void ECN_EMS(const decoder_t &theta_1,
                     const decoder_t &phi_1,
                     const vector<vector<uint16_t>> &ADDDEC,
                     const vector<vector<uint16_t>> &DIVDEC,
                     const uint16_t nm,
                     uint16_t q,
                     const uint16_t coef,
                     decoder_t &theta);

        void tuples_sorter(const decoder_t &tuples,
                           const uint16_t nm,
                           const uint16_t q,
                           decoder_t &sorted_tuples);

        void VN_update(const decoder_t &theta_1,
                       const decoder_t &phi_1,
                       const vector<vector<uint16_t>> &ADDDEC,
                       const vector<vector<uint16_t>> &DIVDEC,
                       uint16_t nm,
                       uint16_t coef,
                       uint16_t hard_decision,
                       uint16_t q,
                       softdata_t ofst,
                       decoder_t &phi);
        void ECN_AEMS_bubble(const decoder_t &theta_1,
                             const decoder_t &phi_1,
                             const decoder_parameters &dec_param,
                             uint16_t coef,
                             const vector<vector<uint16_t>> &ADDDEC,
                             const vector<vector<uint16_t>> &DIVDECC,
                             decoder_t &theta);
        void decode_SC(const decoder_parameters &dec_param,
                       const vector<vector<uint16_t>> &ADDDEC,
                       const vector<vector<uint16_t>> &MULDEC,
                       const vector<vector<uint16_t>> &DIVDEC,
                       vector<vector<decoder_t>> &L,
                       vector<uint16_t> &info_sec_rec);
        void ECN_EMS_L(const decoder_t &theta_1,
                       const decoder_t &phi_1,
                       const vector<vector<uint16_t>> &ADDDEC,
                       const vector<vector<uint16_t>> &DIVDEC,
                       const decoder_parameters &dec_param,
                       const uint16_t coef,
                       decoder_t &theta1,
                       vector<vector<uint16_t>> &Cs1);

        void decode_SC(const decoder_parameters &dec_param,
                       const vector<vector<uint16_t>> &ADDDEC,
                       const vector<vector<uint16_t>> &MULDEC,
                       const vector<vector<uint16_t>> &DIVDEC,
                       vector<vector<decoder_t>> &L,
                       vector<uint16_t> &info_sec_rec,
                       vector<vector<vector<vector<uint16_t>>>> &Cs);

    } // namespace decoding

} // namespace PoAwN

#endif