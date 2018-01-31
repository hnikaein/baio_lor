#include "sequence.h"
#include "../utils/logger.h"
#include <algorithm>
#include <boost/algorithm/string/replace.hpp>
#include <boost/format.hpp>

using namespace std;

Sequence::Sequence(const char *const seq_str, const unsigned long seq_str_len, const char *name,
                   const char *quality_str) : seq_str(seq_str), size(seq_str_len), name(name),
                                              quality_str(quality_str) {}

Sequence::Sequence(const std::string &seq_str) : seq_str(seq_str.c_str()), size(seq_str.size()), name(""),
                                                 quality_str("") {}

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
    return Sequence(new_seq_str, size, name, quality_str);
}

vector<Sequence> Sequence::chunkenize_big_sequence(const vector<Sequence> &seqs, unsigned int chunk_size) {
    auto chunk_diff = static_cast<int>(chunk_size / 2);
    vector<Sequence> result1, result2;
    for (auto seq : seqs) {
        for (int e1 = 0, chunk_i = 0;
             e1 + chunk_diff < seq.size || (chunk_diff >= seq.size && e1 == 0); e1 += chunk_diff, chunk_i++) {
            string new_name = Logger::formatString("%08d_%lu_%s", chunk_i, chunk_size, seq.name);
            auto *new_name_str = new char[new_name.size() + 1];
            strncpy(new_name_str, new_name.c_str(), new_name.size() + 1);
            result1.emplace_back(seq.seq_str + e1, min(static_cast<int>(chunk_size), static_cast<int>(seq.size - e1)),
                                 new_name_str);
        }
    }
    return move(result1);
}

const char *Sequence::get_name_c() {
    return name;
}

void Sequence::write_to_file(const char *file_name, const bool append) {
    auto file = fopen(file_name, append ? "a" : "w");
    char tempch = '>', endlch = '\n';
    fwrite(&tempch, static_cast<size_t>(1), sizeof(char), file);
    fwrite(name, static_cast<size_t>(strlen(name)), sizeof(char), file);
    fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    fwrite(seq_str, static_cast<size_t>(size), sizeof(char), file);
    fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    fclose(file);
}



