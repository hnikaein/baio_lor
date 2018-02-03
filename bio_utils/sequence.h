#include <vector>
#include <string>


#ifndef SEQUENCE_H
#define SEQUENCE_H

class Sequence {
public:
    static std::tuple<std::vector<Sequence>, std::vector<Sequence>, std::vector<int>>
    chunkenize_big_sequence(const std::vector<Sequence> &seqs, unsigned int chunk_size, bool double_needed = false);

    Sequence(const char *seq_str, unsigned long seq_str_len, const char *name = "", const char *quality_str = "");

    Sequence(Sequence &&sequence) noexcept;

    std::string get_name();

    const char *get_name_c();

    Sequence get_reversed();

    int write_to_file(const char *file_name, bool append = false, bool force_write = true);

    ~Sequence();

    const char *seq_str, *quality_str;
    unsigned long size;
    int delete_flag = 0b000;
private:
    const char *name;
};

#endif //SEQUENCE_H
