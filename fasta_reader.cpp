#include <iostream>
#include <fstream>
#include <string>

// Extract the sequence at the given index from a FASTA file
std::string getFastaSequenceByIndex(const std::string& filename, size_t targetIndex) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << "\n";
        return "";
    }

    std::string line;
    std::string currentSequence;
    size_t currentIndex = -1;

    while (std::getline(inFile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            currentIndex++;
            if (currentIndex > targetIndex) break;
            currentSequence.clear();
        } else if (currentIndex == targetIndex) {
            currentSequence += line;
        }
    }

    if (currentIndex < targetIndex) {
        std::cerr << "Error: Index " << targetIndex << " out of range.\n";
        return "";
    }

    return currentSequence;
}
