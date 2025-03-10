#include "GF_tools.h"
#include "struct.h"
#include "tools.h"
#include <fstream>
#include <iomanip>
#include <cstring>
#include <iostream>
#include "init.h"
#include <sstream>
void PoAwN::init::LoadCode(PoAwN::structures::base_code_t &code, float SNR)
{
    code.Rate = (float)code.K / (float)code.N;
    code.reliab_sequence.resize(code.q, 0);

    std::ostringstream fname;
    std::string mat_direct;
    if (code.sig_mod == "bpsk")
        mat_direct = "./matrices/bpsk/N";
    else if (code.sig_mod == "CCSK_BIN")
        mat_direct = "./matrices/ccsk_bin/N";
    else
        mat_direct = "./matrices/ccsk_nb/N";

    fname << mat_direct << code.N << "/mat_N" << code.N << "_GF"
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

void PoAwN::init::LoadBubblesIndcatorlists(PoAwN::structures::decoder_parameters &dec,
     const float SNR,
      const float Pt)
{
    uint16_t n = dec.n, nH = dec.nH, nL = dec.nL;
    std::ostringstream fname;
    std::string mat_direct;
    if (dec.sig_mod == "bpsk")
        mat_direct = "./BubblesPattern/bpsk/N";
    else if (dec.sig_mod == "CCSK_BIN")
        mat_direct = "./BubblesPattern/ccsk_bin/N";
    else
        mat_direct = "./BubblesPattern/ccsk_nb/N";

    fname << mat_direct << dec.N << "/bubbles_N" << dec.N << "_GF" << dec.q
          << "_SNR" << std::fixed << std::setprecision(2) << SNR << "_" << dec.nH << "x" <<  std::fixed << std::setprecision(2) << dec.nL << "_Pt"
          << std::fixed << std::setprecision(2) <<Pt << "_Bt_lsts.txt";
    std::string filename = fname.str();

    std::ifstream file(filename);
    std::vector<std::string> lines;
    std::string line;
    int linecount = 0;
    while (std::getline(file, line))
    {
        if (!line.empty())
        {
            lines.push_back(line);
            linecount++;
        }
    }
    if (linecount != (1u << n) - 1)
    {
        std::cerr << "File data is not enough, it should  contain " << 1u << n - 1 << "lines, each contains the list of (i,j) coordinates of cluster s at layer l!" << std::endl;
        return;
    }

    int line_cnt0,  line_cnt= 0;
    dec.Bubb_Indicator.resize(n);
    for (int l = 0; l < n; l++)
    {
        dec.Bubb_Indicator[l].resize(1 << l);
        line_cnt0 = (1u << l) - 1;
        for (int s = 0; s < (1u << l); s++)
        {
            line_cnt = line_cnt0+s;
            dec.Bubb_Indicator[l][s].resize(2);
            {
                std::istringstream iss(lines[line_cnt]);
                std::vector<int> numbers;
                std::string token;
                while (iss >> token)
                {
                    std::string clean;
                    for (char c : token)
                        if (std::isdigit(c))
                            clean += c;
                        else if (!clean.empty())
                        { 
                            numbers.push_back(std::stoi(clean));
                            clean.clear();
                        }
                    if (!clean.empty())
                        numbers.push_back(std::stoi(clean));
                }

                for (size_t i = 0; i < numbers.size(); i += 2)
                {
                    dec.Bubb_Indicator[l][s][0].push_back(numbers[i]);
                    dec.Bubb_Indicator[l][s][1].push_back(numbers[i + 1]);
                }
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