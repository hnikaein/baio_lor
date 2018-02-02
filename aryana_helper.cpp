#include <cmath>
#include <cstring>
#include <unistd.h>
#include <sstream>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include "aryana_helper.h"
#include "configs.h"
#include "utils/logger.h"

long long total_candidates = 0, best_factor_candidates = 0;

using namespace std;
const char *aryana_helper_ref_file_base_name;
map<int, vector<pair<int, int>>> aryana_helper_results;
vector<Sequence> *aryana_helper_reads;
unordered_set<string> headers;
vector<string> result_lines;
string first_header;


int run_aryana_for_ref(const int ref_num) {
    auto log_chunk = static_cast<int>(log2(CHUNK_SIZES[CHUNK_SIZES_LEN - 1]));
    map<string, int> reads_ref_id;

    auto ref_file_name = Logger::formatString("%s/%s_%d_%d.fasta", CHUNKS_FOLDER_NAME, aryana_helper_ref_file_base_name,
                                              log_chunk, ref_num);

    string reads_string;
    for (auto read_num: aryana_helper_results[ref_num]) {
        Sequence *read = &(*aryana_helper_reads)[read_num.first];
        reads_string.append("@");
        reads_string.append(read->get_name_c());
        reads_string.append("\n");
        reads_string.append(read->seq_str);
        reads_string.append("\n+\n");
        reads_string.append(read->quality_str);
        reads_string.append("\n");
        if (read_num.second > 0)
            reads_ref_id[read->get_name()] = read_num.second;
    }

    int p[2];
    pipe(p);
    write(p[1], (void *) reads_string.c_str(), reads_string.size());
    close(p[1]);
    stdin = fdopen(p[0], "rb");

    unsigned long buf_size = (reads_string.size() + 500 * aryana_helper_results[ref_num].size()) * 3 / 2;
    char buffer[buf_size];
    FILE *tmp_stdout = stdout;
    stdout = fmemopen(buffer, buf_size, "wb");

    aryana_args args{};
    args.discordant = 1;
    args.threads = 1;
    args.potents = 10;
    args.debug = 0;
    args.seed_length = 8; // XXX changed
    args.best_factor = 0.6;
    args.bisulfite = 0;
    args.order = 0;
    args.exactmatch_num = 50;
    args.report_multi = 0;
    args.mismatch_limit = -1;
    args.mismatch_penalty = 5;
    args.gap_open_penalty = 5;
    args.gap_ext_penalty = 3;
    args.out_buffer_factor = 100000;
    args.ignore = ignore_none;
    args.orientation = orien_all;
    args.min_dis = 0;
    args.max_dis = 10000;
    args.reference = strdup(ref_file_name.c_str());
    args.read_file = const_cast<char *>("-");
    args.single = 1;
    args.paired = 0;
    bwa_aln_core2(&args);

    stdout = tmp_stdout;
    istringstream res(buffer);
    getline(res, first_header);
    string line;
    getline(res, line);
    line = line.substr(7, line.find("LN:") - 8);
    size_t sz;
    int offset = max((stoi(line) - 1) * stoi(line.substr(9), &sz) / 2, 0);
    line = "@SQ\tSN:" + line.substr(10 + sz);
    headers.insert(line);
    while (getline(res, line)) {
        size_t pos = 0, prev_pos = 0;
        int j = 0;
        string token, new_line, name;
        while ((pos = line.find('\t', pos + 1)) != string::npos) {
            token = line.substr(prev_pos, pos - prev_pos);
            if (j == 0)
                name = token;
            if (j == 1 && reads_ref_id.count(name))
                token = to_string(stoi(token) + 256);
            if (j == 2)
                token = token.substr(10 + sz);
            if (j == 3)
                token = to_string(stoi(token) + offset);
            new_line.append(token);
            new_line.append("\t");
            prev_pos = pos + 1;
            j++;
        }
        token = line.substr(prev_pos);
        new_line.append(token);
        result_lines.push_back(new_line);
    }
}

void run_aryana(const char *ref_file_base_name, const char *reads_file_name, vector<Sequence> &reads,
                map<int, vector<int>> &results) {
    aryana_helper_ref_file_base_name = ref_file_base_name;
    for (int i = 0; i < reads.size(); ++i)
        for (int j = 0; j < results[i].size(); j++) {
            auto num_ref = results[i][j];
            if (!aryana_helper_results.count(num_ref))
                aryana_helper_results[num_ref] = vector<pair<int, int>>();
            aryana_helper_results[num_ref].emplace_back(i, j);
        }
    aryana_helper_reads = &reads;
    for (const auto &p: aryana_helper_results) {
        printf("%d\n", p.first);
        run_aryana_for_ref(p.first);
    }
//    multiproc(THREADS_COUNT, run_aryana_for_ref, static_cast<int>(reads.size()));
    ofstream result_file((string(reads_file_name) + ".sam").c_str());
    result_file << first_header << endl;
    for (const auto &header :headers)
        result_file << header << endl;
    sort(result_lines.begin(), result_lines.end());
    for (const auto &result_line: result_lines)
        result_file << result_line << endl;
    result_file.close();
}
