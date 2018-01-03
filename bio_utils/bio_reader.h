#include "sequence.h"
#include <vector>
#include <tuple>

#ifndef BIO_READER_H
#define BIO_READER_H

enum FileType {
    FASTA, FASTQ, SAM, BAM
};

std::vector<Sequence> read_from_file(const char *file_name, const FileType &file_type);

std::tuple<std::vector<char *>, std::vector<char *>, std::vector<int>> read_fasta(const char *file_name);

std::tuple<std::vector<char *>, std::vector<char *>, std::vector<int>, std::vector<char *>>
read_fastq(const char *file_name);

#endif //BIO_READER_H
