#ifndef HELPERFUNC_H
#define HELPERFUNC_H

#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem> 


template <typename T>
bool appendToFile(const std::string &filename, const std::vector<std::vector<T>> &table, const bool newsim, const bool newclust)
{
    std::filesystem::path filepath(filename);
    std::filesystem::create_directories(filepath.parent_path());

    std::ofstream file(filename, newsim ? std::ios::out : std::ios::app);
    if (newclust)
    {
        std::string dashes(400, '-');
        file << "\n"
             << dashes << "\n\n";
    }
    else if (!newsim)
        file << "\n";

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }

    for (const auto &row : table)
    {
        for (const auto &value : row)
        {
            file << value << " ";
        }
        file << "\n";
    }
    file.close();
    return true;
}


inline void save_intrinsic_LLR(const std::vector<std::vector<PoAwN::structures::softdata_t>> &chan_LLR, const std::string &filename)
{
    std::ofstream file(filename, std::ios::app);
    for (const auto &symbol_LLR : chan_LLR)
    {
        for (const auto &llr : symbol_LLR)
            file << llr << " ";
        file << "\n";
    }
}

template <typename T>
bool appendSequenceToFile(const std::string &filename, const std::vector<T> &sequence)
{
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) return false;

    for (const auto &symbol : sequence)
        file << symbol << " ";
    file << "\n";

    return true;
}

#endif // HELPERFUNC_H
