#ifndef HELPERFUNC_H
#define HELPERFUNC_H

#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem> // For creating directories

// Function to append data to file, creating directories if they don't exist
template <typename T>
bool appendToFile(const std::string &filename, const std::vector<std::vector<T>> &table,const bool newsim, const bool newclust)
{
    // Ensure the directory exists before trying to open the file
    std::filesystem::path filepath(filename);
    std::filesystem::create_directories(filepath.parent_path()); // Create the directory if it doesn't exist

    // Open file in append mode
    std::ofstream file(filename, newsim ? std::ios::out : std::ios::app);
    if (newclust)
    {
        std::string dashes(400, '-');
        file << "\n"<< dashes << "\n\n";
    }
    else if(!newsim)
        file << "\n";

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return 0;
    }

    // Write table data to the file
    for (const auto &row : table)
    {
        for (const auto &value : row)
        {
            file << value << " ";
        }
        file << "\n"; // Newline after each row
    }
    file.close();
    return 1;
}

#endif // HELPERFUNC_H
