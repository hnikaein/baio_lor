#include "vector_writer.h"
#include <cstdio>
#include <cstring>
#include <sys/stat.h>

using namespace std;

void write_to_file(const char *const file_name, const vector<int> *data, const int size) {
    auto file = fopen(file_name, "wb");
    const int MYBUFSIZE = 10 * BUFSIZ;
    int buffer[MYBUFSIZE], write_size = 0, last_data_size = -1, additional_zeros = 0;
    buffer[write_size++] = size;
    for (int i = 0; i < size; ++i) {
        auto data_size = static_cast<int>(data[i].size());
        if (last_data_size == 0)
            if (last_data_size == data_size) {
                additional_zeros++;
                continue;
            } else {
                memcpy(buffer + (write_size++), &additional_zeros, sizeof(int));
                additional_zeros = 0;
            }
        last_data_size = data_size;
        if (write_size + data_size >= MYBUFSIZE) {
            fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
            write_size = 0;
            if (data_size >= MYBUFSIZE) {
                fwrite(&data_size, static_cast<size_t>(1), sizeof(int), file);
                fwrite(&data[i][0], static_cast<size_t>(data_size), sizeof(int), file);
            }
        }
        if (data_size < MYBUFSIZE) {
            memcpy(buffer + (write_size++), &data_size, sizeof(int));
            if (data_size) {
                memcpy(buffer + write_size, &data[i][0], data_size * sizeof(int));
                write_size += data_size;
            }
        }
    }
    fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
    fclose(file);
}

vector<int> *read_from_file(char const *file_name) {
    struct stat st{};
    auto file = fopen(file_name, "rb");
    if (!file || stat(file_name, &st) != 0)
        throw "file is not present or not readable";
    auto all_file = new int[st.st_size / 4];
    fread(all_file, sizeof(int), static_cast<size_t>(st.st_size / 4), file);
    fclose(file);
    auto buffer_p = all_file;
    int size = *(buffer_p++);
    auto result = new vector<int>[size];
    for (int i = 0; i < size; ++i) {
        int vec_size = *(buffer_p++);
        if (vec_size == 0) {
            int additional_zeros = *(buffer_p++);
            i += additional_zeros;
        } else {
            result[i].assign(buffer_p, buffer_p + vec_size);
            buffer_p += vec_size;
        }
    }
    delete[] all_file;
    return result;
}
