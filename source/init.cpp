#include "GF_tools.h"
#include "struct.h"
#include "tools.h"
#include <fstream>
#include <iomanip>
#include <cstring>
#include <iostream>
#include "init.h"

void PoAwN::init::LoadCode(PoAwN::structures::base_code_t &code, float SNR)
{
    code.Rate = (float)code.K / (float)code.N;
    code.reliab_sequence.resize(code.q, 0);

    // Construct file name using std::ostringstream
    std::ostringstream fname;
    fname << "./matrices/N" << code.N << "/mat_N" << code.N << "_GF"
          << code.q << "_SNR" << std::fixed << std::setprecision(2) << SNR
          << ".txt";

    std::cout << fname.str() << std::endl;
    std::ifstream opfile(fname.str());
    if (!opfile)
    {
        std::cerr << "Sequence Unavailable!!" << std::endl;
        exit(-1010);
    }

    int tmp;
    code.reliab_sequence.resize(code.N);
    for (int i = 0; i < code.N; i++)
    {
        opfile >> tmp;
        code.reliab_sequence[i] = tmp;
    }

    if (code.sig_mod == "bpsk")
    {
        code.polar_coeff.resize(code.n, std::vector<uint16_t>(code.N / 2));
        // Read polar coefficients
        for (int i = 0; i < code.n; i++)
        {
            for (int j = 0; j < code.N / 2; j++)
            {
                opfile >> tmp;
                code.polar_coeff[i][j] = tmp;
            }
        }
    }
    else
    {
        code.polar_coeff.resize(code.n, std::vector<uint16_t>(code.N / 2));
        for (int i = 0; i < code.n; i++)
        {
            for (int j = 0; j < code.N / 2; j++)
            {
                code.polar_coeff[i][j] = 1;
            }
        }
    }
}

void PoAwN::init::LoadTables(PoAwN::structures::base_code_t &code,
                PoAwN::structures::table_GF &table,
                const uint16_t *GF_polynom_primitive)
{
    uint16_t prim_pol = GF_polynom_primitive[code.p - 2];
    PoAwN::GFtools::GF_bin_seq_gen(code.q, prim_pol, table.BINGF, table.BINDEC);
    PoAwN::GFtools::GF_bin2GF(table.BINGF, table.GFDEC, table.DECGF);
    PoAwN::GFtools::GF_add_mat_gen(table.DECGF, table.ADDGF, table.ADDDEC);
    PoAwN::GFtools::GF_mul_mat_gen(table.GFDEC, table.MULGF, table.MULDEC);
    PoAwN::GFtools::GF_div_mat_gen(table.GFDEC, table.DIVGF, table.DIVDEC);
}