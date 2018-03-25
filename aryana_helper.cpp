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
#include <fcntl.h>

long long total_candidates = 0, best_factor_candidates = 0;

extern Logger *logger;
using namespace std;
const char *aryana_helper_ref_file_base_name;
map<int, map<int, int>> aryana_helper_results;
vector<Sequence> *aryana_helper_reads;
unordered_set<string> headers;
vector<string> result_lines;
string first_header;
int total_unassigned = 0;


int run_aryana_for_ref(const int ref_num, int threads_count) {
    logger->debug("begin of run_aryana for ref_num: %d", ref_num);
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
        reads_string.append(read->seq_str, 0, 40000);
        reads_string.append("\n+\n");
        reads_string.append(read->quality_str, 0, 40000);
        reads_string.append("\n");
        if (read_num.second > 0)
            reads_ref_id[read->get_name()] = read_num.second;
    }

    int p[2];
    pipe(p);
    fcntl(p[1], F_SETPIPE_SZ, reads_string.size() + 100);
    write(p[1], (void *) reads_string.c_str(), reads_string.size());
    close(p[1]);
    stdin = fdopen(p[0], "rb");

    auto reads_size = static_cast<int>(aryana_helper_results[ref_num].size());
    unsigned long buf_size = (reads_string.size() + 500 * reads_size) * 3 / 2;
    char buffer[buf_size] = {};
    FILE *tmp_stdout = stdout;

    aryana_args args{};
    args.discordant = 1;
    args.threads = min(threads_count, max(reads_size / 7, 1));
    args.potents = 100; // XXX changed
    args.debug = logger->log_level - 4;
    args.seed_length = 10; // XXX changed
    args.best_factor = 0.6;
    args.bisulfite = 0;
    args.order = 0;
    args.exactmatch_num = 50;
    args.report_multi = 0;
    args.mismatch_limit = -1;
    args.mismatch_penalty = 3; // XXX changed
    args.gap_open_penalty = 3; // XXX changed
    args.gap_ext_penalty = 2; // XXX changed
    args.out_buffer_factor = 100000;
    args.ignore = ignore_none;
    args.orientation = orien_all;
    args.min_dis = 0;
    args.max_dis = 10000;
    args.reference = strdup(ref_file_name.c_str());
    args.read_file = const_cast<char *>("-");
    args.single = 1;
    args.paired = 0;
    args.indel_ratio_between_seeds = 2; // XXX: Added by me
    args.platform = pacbio; // XXX changed
    logger->debug("begin of aryana for ref_num and count: %d -> %d", ref_num, reads_size);
    stdout = fmemopen(buffer, buf_size, "wb");
    bwa_aln_core2(&args);
    fflush(stdout);
    close(p[0]);
    stdout = tmp_stdout;
    free(args.reference);
    logger->debug("end of aryana for ref_num: %d", ref_num);

    if (!buffer[0]) {
        logger->info("aryana have no result for ref_num: %d and result count: %d", ref_num, reads_size);
        if (reads_size == 1)
            logger->debug("which is %d", aryana_helper_results[ref_num][0]);
        return 0;
    }
    istringstream res(buffer);
    getline(res, first_header);
    string line;
    getline(res, line);
    line = line.substr(7, line.find("LN:") - 8);
    size_t sz;
    int offset = max((stoi(line) - 1) * stoi(line.substr(9), &sz) / 2, 0);
    line = "@SQ\tSN:" + line.substr(10 + sz) + "\tLN:0";
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
            if (j == 2) {
                if (token[0] == '*')
                    break;
                token = token.substr(10 + sz);
            }
            if (j == 3)
                token = to_string(stoi(token) + offset);
            new_line.append(token);
            new_line.append("\t");
            prev_pos = pos + 1;
            j++;
        }
        if (token.size() == 1) {
            total_unassigned++;
            continue;
        }
        token = line.substr(prev_pos);
        new_line.append(token);
        result_lines.push_back(new_line);
    }
    return 0;
}

void run_aryana(const char *ref_file_base_name, const char *output_file_name, vector<Sequence> &reads,
                map<int, vector<int>> &results, int threads_count) {
    aryana_helper_ref_file_base_name = ref_file_base_name;
    for (int i = 0; i < reads.size(); ++i)
        for (int j = 0; j < results[i].size(); j++) {
            auto num_ref = results[i][j];
            if (!aryana_helper_results.count(num_ref))
                aryana_helper_results[num_ref] = map<int, int>();
            if (!aryana_helper_results[num_ref].count(i))
                aryana_helper_results[num_ref][i] = j;
            else
                logger->debug("multiple candidate for read number: %d, ref number: %d", i, num_ref);
        }
    aryana_helper_reads = &reads;
    for (const auto &p: aryana_helper_results)
        run_aryana_for_ref(p.first, threads_count);
    logger->info("total unassigned: %d", total_unassigned);
    ofstream result_file(output_file_name); // TODO convert to fwrite
    result_file << first_header << endl;
    for (const auto &header :headers)
        result_file << header << endl;
    sort(result_lines.begin(), result_lines.end());
    for (const auto &result_line: result_lines)
        result_file << result_line << endl;
    result_file.close();
}

