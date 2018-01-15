#include "bio_reader.h"
#include <sys/stat.h>

using namespace std;

vector<Sequence> read_from_file(const char *const file_name, const FileType &file_type) {
    vector<char *> names, seqs, quality;
    vector<int> lens;
    if (file_type == FASTA)
        tie(names, seqs, lens) = read_fasta(file_name);
    else if (file_type == FASTQ)
        tie(names, seqs, lens, quality) = read_fastq(file_name);
    auto l = vector<Sequence>();
    if (quality.empty())
        for (int i = 0; i < names.size(); i++)
            l.emplace_back(seqs[i], lens[i], names[i], "");
    else
        for (int i = 0; i < names.size(); i++)
            l.emplace_back(seqs[i], lens[i], names[i], quality[i]);
    return l;
}


std::tuple<std::vector<char *>, std::vector<char *>, std::vector<int>> read_fasta(const char *const file_name) {
    vector<char *> seqs, names;
    vector<int> lens;
    struct stat st{};
    auto file = fopen(file_name, "r");
    if (!file || stat(file_name, &st) != 0)
        throw "fasta file not found or not readable";

    auto result = new char[st.st_size];
    bool name_state = true;
    int result_i = 0, last_result_i = 0;
    char buffer[BUFSIZ], last_ch, ch = '\0';
    size_t read_len;

    fread(buffer, 1, 1, file);
    names.push_back(result);

    while ((read_len = fread(buffer, 1, BUFSIZ, file))) {
        for (int i = 0; i < read_len; i++) {
            last_ch = ch;
            ch = buffer[i];
            if (ch == '>') {
                result[result_i++] = '\0';
                name_state = true;
                names.push_back(result + result_i);
                lens.push_back(result_i - last_result_i - 1);
                continue;
            } else if (ch == '\n' || ch == '\r') {
                if (last_ch != '\n' && last_ch != '\r' && name_state) {
                    result[result_i++] = '\0';
                    name_state = false;
                    seqs.push_back(result + result_i);
                    last_result_i = result_i;
                }
                continue;
            }
            result[result_i++] = name_state ? ch : static_cast<char>(toupper(ch));
        }
    }
    result[result_i++] = '\0';
    lens.push_back(result_i - last_result_i);
    fclose(file);
    return make_tuple(move(names), move(seqs), move(lens));
}

std::tuple<std::vector<char *>, std::vector<char *>, std::vector<int>, std::vector<char *>>
read_fastq(const char *const file_name) {
    vector<char *> seqs, names, quality;
    vector<int> lens;
    auto file = fopen(file_name, "r");
    if (!file)
        return make_tuple(move(names), move(seqs), move(lens), move(quality));
    struct stat st{};
    if (stat(file_name, &st) != 0)
        throw "fastq file not found or not readable";

    auto result = new char[st.st_size];

    int state = 0;
    int result_i = 0, last_result_i = 0;
    char buffer[BUFSIZ], last_ch, ch = '\0';
    size_t read_len;
    names.push_back(result + 1);


    while ((read_len = fread(buffer, 1, BUFSIZ, file))) {
        for (int i = 0; i < read_len; i++) {
            last_ch = ch;
            ch = buffer[i];
            if (ch == '\n' || ch == '\r') {
                if (last_ch == '\n' || last_ch == '\r')
                    continue;
                if (state == 0 || state == 1 || state == 3)
                    result[result_i++] = '\0';
                state = (state + 1) % 4;
                if (state == 0)
                    names.push_back(result + result_i + 1);
                if (state == 1) {
                    seqs.push_back(result + result_i);
                    last_result_i = result_i;
                }
                if (state == 2)
                    lens.push_back(result_i - last_result_i - 1);
                if (state == 3)
                    quality.push_back(result + result_i);
            } else if (state != 2)
                result[result_i++] = state == 1 ? static_cast<char>(toupper(ch)) : ch;
        }
    }
    fclose(file);
    result[result_i] = '\0';
    if (names.size() > lens.size()) // in case of a trailing newline
        names.pop_back();
    return make_tuple(move(names), move(seqs), move(lens), move(quality));
}
