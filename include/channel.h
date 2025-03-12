#ifndef CHANNEL
#define CHANNEL
#include "struct.h"

#include <cstdint>
#include <vector>
#include"struct.h"

namespace PoAwN
{
    namespace channel
    {
        using structures::decoder_parameters;
        using structures::decoder_t;
        using structures::vector;
        using structures::table_GF;
        void EncodeChanBPSK_BinCCSK(const decoder_parameters &dec_param,
            const table_GF &table,
            const float SNR,
            const vector<vector<uint16_t>> &bin_table,
            vector<decoder_t> &chan_LLR_sorted,
            vector<uint16_t> &KSYMB);

            void EncodeChanGF_CCSK(const decoder_parameters &dec_param,
                const table_GF &table,
                const float SNR,
                const vector<vector<uint16_t>> &CCSK_rotated_codes,
                vector<decoder_t> &chan_LLR_sorted,
                vector<uint16_t> &KSYMB);
    }
}

#endif