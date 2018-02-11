#include "sequence.h"
#include "../utils/logger.h"
#include <cstring>
#include <sys/stat.h>

using namespace std;

Sequence::Sequence(const char *const seq_str, const unsigned long seq_str_len, const char *name,
                   const char *quality_str) : seq_str(seq_str), size(seq_str_len), name(name),
                                              quality_str(quality_str) {}

string Sequence::get_name() {
    return string(name);
}

Sequence Sequence::get_reversed() {
    auto *new_seq_str = new char[size + 1];
    new_seq_str -= 1;
    for (int i = 0; i < size; ++i) {
        if (seq_str[i] == 'A')
            new_seq_str[size - i] = 'T';
        if (seq_str[i] == 'T')
            new_seq_str[size - i] = 'A';
        if (seq_str[i] == 'C')
            new_seq_str[size - i] = 'G';
        if (seq_str[i] == 'G')
            new_seq_str[size - i] = 'C';
    }
    new_seq_str += 1;
    new_seq_str[size] = '\0';
    Sequence sequence(new_seq_str, size, name, quality_str);
    sequence.delete_flag = 0b010;
    return move(sequence);
}

tuple<vector<Sequence>, vector<Sequence>, vector<int>>
Sequence::chunkenize_big_sequence(const vector<Sequence> &seqs, unsigned int chunk_size, const bool double_needed) {
    auto chunk_diff = static_cast<int>(chunk_size / 2);
    vector<Sequence> chunks, double_chunks;
    vector<int> starts;
    for (auto &seq : seqs) {
        starts.push_back(static_cast<int>(chunks.size()));
        for (int e1 = 0, chunk_i = 0;
             e1 + chunk_diff < seq.size || (chunk_diff >= seq.size && e1 == 0); e1 += chunk_diff, chunk_i++) {
            auto new_name = Logger::formatString("%08d_%lu_%s", chunk_i, chunk_size, seq.name);
            Sequence s1(seq.seq_str + e1, min(static_cast<unsigned>(chunk_size), static_cast<unsigned>(seq.size - e1)),
                        strdup(new_name.c_str()));
            s1.delete_flag = 0b100;
            chunks.push_back(move(s1));
            if (double_needed) {
                const int sequence_real_start = max(e1 - chunk_diff, 0);
                Sequence s2(seq.seq_str + sequence_real_start,
                            min(static_cast<unsigned long>(e1) + chunk_size + chunk_diff, seq.size) -
                            sequence_real_start, strdup(new_name.c_str()));
                s2.delete_flag = 0b100;
                double_chunks.push_back(move(s2));
            }
        }
    }
    starts.push_back(static_cast<int>(chunks.size()));
    return {move(chunks), move(double_chunks), move(starts)};
}

const char *Sequence::get_name_c() {
    return name;
}

int Sequence::write_to_file(const char *file_name, const bool append, const bool force_write) {
    if (!force_write) {
        struct stat st{};
        if (stat(file_name, &st) == 0)
            return 1;
    }
    auto file = fopen(file_name, append ? "a" : "w");
    char tempch = '>', endlch = '\n';
    fwrite(&tempch, static_cast<size_t>(1), sizeof(char), file);
    fwrite(name, static_cast<size_t>(strlen(name)), sizeof(char), file);
    fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    fwrite(seq_str, static_cast<size_t>(size), sizeof(char), file);
    fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    fclose(file);
    return 0;
}

Sequence::~Sequence() {
    if (delete_flag & 0b100)
        free((void *) name);
    if (delete_flag & 0b010)
        delete[] seq_str;
    if (delete_flag & 0b001)
        delete[] quality_str;
}

Sequence::Sequence(Sequence &&sequence) noexcept : name(sequence.name), quality_str(sequence.quality_str),
                                                   seq_str(sequence.seq_str), size(sequence.size),
                                                   delete_flag(sequence.delete_flag) {
    sequence.name = sequence.quality_str = sequence.seq_str = nullptr;
}
