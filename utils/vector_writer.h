#include <vector>

#ifndef C_LONG_READ_ALIGNER_VECTOR_WRITER_H
#define C_LONG_READ_ALIGNER_VECTOR_WRITER_H

void write_to_file(char const *file_name, const std::vector<int *> &data, int row_size);

std::vector<int *> read_from_file(char const *file_name);

#endif //C_LONG_READ_ALIGNER_VECTOR_WRITER_H
