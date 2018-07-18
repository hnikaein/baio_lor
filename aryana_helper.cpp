#include "aryana_helper.h"
#include <vector>
#include <map>
#include <cstring>
#include <unistd.h>
#include <sstream>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <fcntl.h>
#include <sys/stat.h>
#include "configs.h"
#include "utils/logger.h"
#include "utils/multiproc.h"
#include "aryana/main.h"
#include "bio_utils/sequence.h"

using namespace std;

long long total_candidates = 0, best_factor_candidates = 0;
extern Logger *logger;
extern int threads_count;
extern char *ref_file_name;
extern char *output_file_name;
extern vector<Sequence> reads;
extern vector<Sequence> genome_double_chunks;
extern map<int, vector<int>> lor_results;

char *index_suffixes[4] = {const_cast<char *>(".bwt"), const_cast<char *>(".sa"), const_cast<char *>(".bin"),
                           const_cast<char *>(".ann")};
unsigned int index_buffer_size = CHUNK_SIZES[CHUNK_SIZES_LEN - 1] * 4;
map<int, map<int, int>> reads_in_ref;
pthread_mutex_t output_lock, index_read_lock[4];
string first_header;
unordered_set<string> headers;
vector<string> result_lines;
int total_unassigned = 0;
int index_size;
long *poses;
FILE *index_files[4];


int run_aryana_for_ref(const int ref_num) {
    logger->debug("begin of run_aryana for ref_num: %d", ref_num);
    unordered_set<string> alternative_results;

    string reads_string;
    for (auto read_num: reads_in_ref[ref_num]) {
        Sequence *read = &(reads[read_num.first]);
        reads_string.append("@");
        reads_string.append(read->get_name_c());
        reads_string.append("\n");
        reads_string.append(read->seq_str, 0, 40000);
        reads_string.append("\n+\n");
        reads_string.append(read->quality_str, 0, 40000);
        reads_string.append("\n");
        if (read_num.second > 0)
            alternative_results.insert(read->get_name());
    }

    int p[2];
    pipe(p);
    fcntl(p[1], F_SETPIPE_SZ, reads_string.size() + 100);
    write(p[1], (void *) reads_string.c_str(), reads_string.size());
    close(p[1]);
    auto stdin_file = fdopen(p[0], "rb");

    auto reads_size = static_cast<int>(reads_in_ref[ref_num].size());
    unsigned long buf_size = (reads_string.size() + 500 * reads_size) * 3 / 2;
    char buffer[buf_size] = {};
    auto stdout_file = fmemopen(buffer, buf_size, "wb");

    char index_buffer[4][index_buffer_size];
    FILE *simulated_index_file[4];
    size_t index_files_size[4];
    for (int i = 0; i < 4; ++i) {
        pthread_mutex_lock(index_read_lock + i);
        fseek(index_files[i], poses[ref_num * 4 + i], SEEK_SET);
        index_files_size[i] = fread(index_buffer[i], 1,
                                    ref_num == index_size - 1 ? index_buffer_size : static_cast<size_t>(
                                            poses[(ref_num + 1) * 4 + i] - poses[ref_num * 4 + i]),
                                    index_files[i]);
        pthread_mutex_unlock(index_read_lock + i);
        simulated_index_file[i] = fmemopen(index_buffer[i], index_files_size[i], i < 3 ? "rb" : "r");
    }

    aryana_args args{};
    args.discordant = 1;
    args.threads = 1; //min(threads_count, max(reads_size / 7, 1));
    args.potents = 100; // XXX changed
    args.debug = logger->log_level - 5;
    args.seed_length = 10; // XXX changed
    args.best_factor = 0.6;
    args.bisulfite = 0;
    args.order = 0;
    args.exactmatch_num = 50;
    args.report_multi = 0;
    args.mismatch_limit = -1;
    args.mismatch_penalty = mismath_penalty; // XXX changed
    args.gap_open_penalty = gap_open_penalty; // XXX changed
    args.gap_ext_penalty = gap_extend_penalty; // XXX changed
    args.out_buffer_factor = 100000;
    args.ignore = ignore_none;
    args.orientation = orien_all;
    args.min_dis = 0;
    args.max_dis = 10000;
    args.reference = const_cast<char *>("");
    args.index_files = simulated_index_file; //TODO
    args.index_files_size = index_files_size;
    args.read_file = new char[10];
    args.read_file[0] = '+';
    for (int i = 1; i < 10; i++)
        args.read_file[i] = '\0';
    memcpy(args.read_file + 2, &stdin_file, sizeof(stdin_file));
    args.stdout_file = stdout_file;
    args.single = 1;
    args.paired = 0;
    args.indel_ratio_between_seeds = 2; // XXX: Added by me
    args.platform = pacbio; // XXX changed
    logger->debug("begin of aryana for ref_num and count: %d -> %d", ref_num, reads_size);
    bwa_aln_core2(&args);
    fflush(stdout_file);
    close(p[0]);
    delete[] args.read_file;
    logger->debug("end of aryana for ref_num: %d", ref_num);

    if (!buffer[0]) {
        logger->info("aryana have no result for ref_num: %d and result count: %d", ref_num, reads_size);
        if (reads_size == 1)
            logger->debug("which is %d", reads_in_ref[ref_num][0]);
        return 0;
    }
    istringstream res(buffer);
    pthread_mutex_lock(&output_lock);
    getline(res, first_header);
    string line;
    getline(res, line);
    line = line.substr(7, line.find("LN:") - 8);
    size_t sz;
    int offset = max((stoi(line) - 1) * stoi(line.substr(9), &sz) / 2, 0);
    line = "@SQ\tSN:" + line.substr(10 + sz) + "\tLN:0";
    headers.insert(line);
    pthread_mutex_unlock(&output_lock);
    while (getline(res, line)) {
        size_t pos = 0, prev_pos = 0;
        int j = 0;
        string token, new_line, name;
        while ((pos = line.find('\t', pos + 1)) != string::npos) {
            token = line.substr(prev_pos, pos - prev_pos);
            if (j == 0)
                name = token;
            if (j == 1 && alternative_results.count(name)) {
                int token_int = stoi(token);
                if (token_int < 256)
                    token_int += 256;
                token = to_string(token_int);
            }
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
            pthread_mutex_lock(&output_lock);
            total_unassigned++;
            pthread_mutex_unlock(&output_lock);
            continue;
        }
        token = line.substr(prev_pos);
        new_line.append(token);
        pthread_mutex_lock(&output_lock);
        result_lines.push_back(new_line);
        pthread_mutex_unlock(&output_lock);
    }
    return 0;
}

void run_aryana() {
    for (int read_i = 0; read_i < reads.size(); ++read_i)
        for (int j = 0; j < lor_results[read_i].size(); j++) {
            auto num_ref = lor_results[read_i][j];
            if (!reads_in_ref.count(num_ref))
                reads_in_ref[num_ref] = map<int, int>();
            if (!reads_in_ref[num_ref].count(read_i))
                reads_in_ref[num_ref][read_i] = j;
            else
                logger->debug("multiple candidate for read number: %d, ref number: %d", read_i, num_ref);
        }
    vector<int> active_refs;
    for (const auto &p: reads_in_ref)
        active_refs.push_back(p.first);

    auto index_of_index = fopen(Logger::formatString("%s.ind", ref_file_name).c_str(), "rb");
    fread(&index_size, sizeof(int), 1, index_of_index);
    poses = new long[index_size * 4];
    fread(poses, sizeof(long), static_cast<size_t>(index_size * 4), index_of_index);
    fclose(index_of_index);
    for (int i = 0; i < 4; ++i)
        index_files[i] = fopen(Logger::formatString("%s%s", ref_file_name, index_suffixes[i]).c_str(), "rb");

    multiproc(threads_count, run_aryana_for_ref, active_refs);

    for (int i = 0; i < 4; ++i)
        fclose(index_files[i]);
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

int aryana_index_part(const int i) {
    char src_file_base_name[strlen(ref_file_name) + 10];
    sprintf(src_file_base_name, "%s.%d", ref_file_name, i);
    genome_double_chunks[i].write_to_file(src_file_base_name, false, false);
    char *argv[] = {const_cast<char *>("index"), src_file_base_name};
    bwa_index(2, argv);
    argv[0] = const_cast<char *>("fa2bin");
    fa2bin(2, argv);
}

void create_aryana_index() {
    auto bwt_file_name = Logger::formatString("%s.bwt", ref_file_name);
    struct stat st{};
    if (stat(bwt_file_name.c_str(), &st) == 0)
        return;
    char *dumb_suffixes[3] = {const_cast<char *>(""), const_cast<char *>(".pac"), const_cast<char *>(".amb")};
    FILE *dst_files[4];
    for (int i = 0; i < 4; ++i)
        dst_files[i] = fopen(Logger::formatString("%s%s", ref_file_name, index_suffixes[i]).c_str(), "wb");
    auto index_of_index = fopen(Logger::formatString("%s.ind", ref_file_name).c_str(), "wb");

    char buffer[index_buffer_size];
    char *src_file_base_name = (char *) calloc(strlen(ref_file_name) + 10, 1);
    int genome_chunk_length = static_cast<int>(genome_double_chunks.size());
    fwrite(&genome_chunk_length, sizeof(int), 1, index_of_index);
    for (int i = 0; i < genome_chunk_length; i += 100) {
        int part_end = min(i + 100, genome_chunk_length);
        multiproc(1, aryana_index_part, part_end, i);
        for (int ii = i; ii < part_end; ++ii) {
            sprintf(src_file_base_name, "%s.%d", ref_file_name, ii);
            long poses[4];
            for (int j = 0; j < 4; ++j) {
                poses[j] = ftell(dst_files[j]);
                auto src_file_name = Logger::formatString("%s%s", src_file_base_name, index_suffixes[j]);
                auto src = fopen(src_file_name.c_str(), "rb");
                size_t size = static_cast<int>(fread(buffer, 1, index_buffer_size, src));
                fwrite(buffer, 1, size, dst_files[j]);
                fclose(src);
                remove(src_file_name.c_str());
            }
            for (int j = 0; j < 3; ++j)
                remove(Logger::formatString("%s%s", src_file_base_name, dumb_suffixes[j]).c_str());
            fwrite(poses, sizeof(long), 4, index_of_index);
        }
    }
    fclose(index_of_index);
    for (int i = 0; i < 4; ++i)
        fclose(dst_files[i]);
    free(src_file_base_name);
}

