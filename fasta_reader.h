#ifndef FASTA_READER_H
#define FASTA_READER_H
#include <string>

std::string getFastaSequenceByIndex(const std::string& filename, size_t targetIndex);

#endif