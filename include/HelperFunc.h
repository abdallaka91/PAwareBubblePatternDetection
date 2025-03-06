#ifndef HELPERFUNC_H
#define HELPERFUNC_H

#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>  // For creating directories

// Function to append data to file, creating directories if they don't exist
template <typename T>
void appendToFile(const std::string& filename, const std::vector<std::vector<T>>& table) {
    // Ensure the directory exists before trying to open the file
    std::filesystem::path filepath(filename);
    std::filesystem::create_directories(filepath.parent_path());  // Create the directory if it doesn't exist

    // Open file in append mode
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write table data to the file
    for (const auto& row : table) {
        for (const auto& value : row) {
            file << value << " ";
        }
        file << "\n";  // Newline after each row
    }

    file.close();
    std::cout << "Data appended to " << filename << std::endl;
}

#endif // HELPERFUNC_H
