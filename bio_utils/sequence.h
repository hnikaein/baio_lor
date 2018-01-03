#include <vector>
#include <string>


#ifndef SEQUENCE_H
#define SEQUENCE_H

class Sequence {
public:
    static std::tuple<std::vector<Sequence>, std::vector<Sequence>>
    chunkenize_big_sequence(const std::vector<Sequence> &seqs, unsigned int chunk_size);

    Sequence(const char *seq_str, unsigned long seq_str_len, const char *name = "", const char *quality_str = "");

    Sequence(const std::string &seq_str);

    std::string get_name();

    const char *get_name_c();

    Sequence get_reversed();

    const char *const seq_str;
    unsigned long size;
private:
    const char *name, *quality_str;
};

#endif //SEQUENCE_H
