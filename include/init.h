#ifndef INIT_H
#define INIT_H

#include "struct.h"
#include <vector>
namespace PoAwN
{
    namespace init
    {
        void LoadCode(PoAwN::structures::base_code_t &code, float SNR);
        void LoadTables(PoAwN::structures::base_code_t &code,
                        PoAwN::structures::table_GF &table,
                        const uint16_t *GF_polynom_primitive);
    } // namespace init
} // namespace PoAwN
#endif