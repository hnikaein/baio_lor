#include "vector_writer.h"
#include <cstdio>
#include <cstring>

using namespace std;

void write_to_file(const char *const file_name, const vector<int *> &data, const int row_size) {
    auto file = fopen(file_name, "wb");
    int buffer[10 * BUFSIZ], write_size = 0;
    buffer[write_size++] = static_cast<int>(data.size());
    buffer[write_size++] = row_size;
    for (int *vec: data) {
        memcpy(buffer + write_size, vec, row_size * sizeof(int));
        write_size += row_size;
        if (write_size > BUFSIZ) {
            fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
            write_size = 0;
        }
    }
    fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
    fclose(file);
}

vector<int *> read_from_file(char const *file_name) {
    auto file = fopen(file_name, "rb");
    if (!file)
        return vector<int *>();
    int buffer[2];
    fread(buffer, sizeof(int), 2, file);
    auto all_file = new int[buffer[0] * buffer[1]];
    fread(all_file, sizeof(int), static_cast<size_t>(buffer[0] * buffer[1]), file);
    fclose(file);
    vector<int *> result(static_cast<unsigned long>(buffer[0]));
    for (int i = 0; i < buffer[0]; ++i)
        result[i] = all_file + i * buffer[1];
    return result;
}
