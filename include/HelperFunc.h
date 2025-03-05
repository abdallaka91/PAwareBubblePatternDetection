#ifndef HELPERFUNC_H
#define HELPERFUNC_H

#include <vector>
#include <fstream>
#include <iostream>

// Function Template (Must Be Defined in Header)
template <typename T>
void appendToFile(const std::string& filename, const std::vector<std::vector<T>>& table) {
    std::ofstream file(filename, std::ios::app);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    file << "\n";
    for (const auto& row : table) {
        for (const auto& elem : row) {
            file << elem << " ";
        }
        file << "\n";
    }
    file.close();
}

#endif // HELPERFUNC_H
