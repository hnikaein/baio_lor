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
    auto s = string(seq_str);
    reverse(s.begin(), s.end());
    boost::replace_all(s, "A", "$");
    boost::replace_all(s, "T", "A");
    boost::replace_all(s, "$", "T");
    boost::replace_all(s, "C", "$");
    boost::replace_all(s, "G", "C");
    boost::replace_all(s, "$", "G");
    auto *seq_str = new char[size + 1];
    strncpy(seq_str, s.c_str(), size + 1);
    return Sequence(seq_str, size, name, quality_str);
}

vector<Sequence> Sequence::chunkenize_big_sequence(const vector<Sequence> &seqs, unsigned int chunk_size) {
    auto chunk_diff = static_cast<int>(chunk_size / 2);
    vector<Sequence> result1, result2;
    for (auto seq : seqs)
        for (int e1 = 0, chunk_i = 0; e1 + chunk_diff < seq.size; e1 += chunk_diff, chunk_i++) {
            string new_name = Logger::formatString("%08d_%lu_%s", chunk_i, chunk_size, seq.name);
            auto *new_name_str = new char[new_name.size() + 1];
            strncpy(new_name_str, new_name.c_str(), new_name.size() + 1);
            result1.emplace_back(seq.seq_str + e1, min(static_cast<int>(chunk_size), static_cast<int>(seq.size - e1)),
                                 new_name_str);
            const int sequence_real_start = max(e1 - chunk_diff, 0);
//            result2.emplace_back(seq.seq_str + sequence_real_start, e1 + chunk_size + chunk_diff - sequence_real_start,
//                                 new_name_str);
        }
    return result1;
}

const char *Sequence::get_name_c() {
    return name;
}



