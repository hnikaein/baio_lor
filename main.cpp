#include "configs.h"
#include "hash_utils.h"
#include "basket_min_hash.h"
#include "bio_utils/bio_reader.h"
#include "utils/time_profile.h"
#include "utils/logger.h"
#include "utils/vector_writer.h"
#include "utils/heap.h"
#include "utils/multiproc.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <sys/stat.h>
#include "aryana/main.h"
#include "aryana_helper.h"

using namespace std;

//vector<Sequence> ref_chunk_sequences[CHUNK_SIZES_LEN];
vector<Sequence> reads, ref_genome;
vector<int> *ref_hash_sketchs[CHUNK_SIZES_LEN];
double alt_matchs_ratio = ALT_MATCHS_RATIO_DEFAULT;
vector<int> genome_parts_starts[CHUNK_SIZES_LEN];
BasketMinHash *similarity_claz;
map<int, vector<int>> total_results;
char *ref_file_base_name;
vector<Sequence> genome_double_chunks;
unsigned int log_chunk;
extern Logger *logger;


int create_aryana_indexes(const int part) {
    auto file_name_str = Logger::formatString("%s/%s_%d_%d.fasta", CHUNKS_FOLDER_NAME, ref_file_base_name,
                                              log_chunk, part);
    char *file_name = strdup(file_name_str.c_str());
    if (!genome_double_chunks[part].write_to_file(file_name, false, false)) {
        char *argv[] = {const_cast<char *>("index"), file_name};
        bwa_index(2, argv);
        argv[0] = const_cast<char *>("fa2bin");
        fa2bin(2, argv);
    }
    free(file_name);
    return 0;
}

string vector_str(int *vec, unsigned long size) {
    string s;
    for (int i = 0; i < size; ++i) {
        s.append(",");
        s.append(to_string(vec[i]));
    }
    return s;
}

void write_results(const char *ref_file_base_name, const char *reads_file_name) {
    char endlch = '\n';
    auto log_chunk = static_cast<int>(log2(CHUNK_SIZES[CHUNK_SIZES_LEN - 1]));

    auto chunks_address = Logger::formatString("%s/%s_%d_%%d.fasta", CHUNKS_FOLDER_NAME, ref_file_base_name, log_chunk);

    auto file = fopen((string(reads_file_name) + ".lra").c_str(), "w");

    fwrite(chunks_address.c_str(), static_cast<size_t>(chunks_address.size()), sizeof(char), file);
    fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    fwrite(reads_file_name, static_cast<size_t>(strlen(reads_file_name)), sizeof(char), file);
    fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    for (int i = 0; i < reads.size(); ++i) {
        string line = vector_str(total_results[i].data(), total_results[i].size());
        fwrite(line.c_str() + 1, static_cast<size_t>(line.size() - 1), sizeof(char), file);
        fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    }
    fclose(file);
}

vector<pair<int, int>>
align_read__find_res(const int chunk_i, const int *sketch_read, const int *sketch_read_reverse) {
    auto temp_result = vector<pair<int, int>>();
    int max_score = 0;
    int max_alt_matchs = 0;
    for (auto sketch:{sketch_read, sketch_read_reverse}) {
        max_alt_matchs += 2 * MAX_ALT_MATCHS;
        auto scores = Heap(genome_parts_starts[chunk_i].back());
        for (int i = 0; i < SKETCH_SIZE; ++i)
            for (int j:ref_hash_sketchs[chunk_i][sketch[i]])
                scores.inc_element(j);
        while (scores.has_element() && temp_result.size() < max_alt_matchs) {
            auto p = scores.pop();
            if (p.first < alt_matchs_ratio * max_score)
                break;
            temp_result.push_back(p);
            if (p.first > max_score)
                max_score = p.first;
        }
    }
    sort(temp_result.begin(), temp_result.end(),
         [](const pair<int, int> &a, const pair<int, int> &b) { return a.first < b.first; });
    unsigned long temp_result_size = temp_result.size();
    if (temp_result_size > 0) {
        int first_score = temp_result[temp_result_size - 1].first;
        int i;
        for (i = static_cast<int>(temp_result_size - 2); i >= 0; --i)
            if (temp_result[i].first < alt_matchs_ratio * first_score)
                break;
        auto result = vector<pair<int, int>>(temp_result_size - i - 1);
        for (int k = i + 1; k < temp_result_size; ++k)
            result[k - i - 1] = move(temp_result[k]);
        return result;
    } else
        return vector<pair<int, int>>();
}

int tot_res = 0;

int align_read(const int read_i) {
    Sequence &read = reads[read_i];
    unsigned int chunk_i = 0;
    while (chunk_i < CHUNK_SIZES_LEN and CHUNK_SIZES[chunk_i] < read.size)
        chunk_i += 1;
    if (chunk_i == CHUNK_SIZES_LEN)
        chunk_i -= 1;
    const unsigned int chunk_size = CHUNK_SIZES[chunk_i];

    auto sketch_read = similarity_claz->get_sketch(read, chunk_i, GINGLE_LENGTH, GAP_LENGTH);
    if (sketch_read[0] == -1) {
        logger->info("removing read %d", read_i);
        return 0;
    }
    auto sketch_read_reverse = similarity_claz->get_sketch(read.get_reversed(), chunk_i, GINGLE_LENGTH, GAP_LENGTH);
    if (sketch_read_reverse[0] == -1) {
        logger->info("removing read %d", read_i);
        return 0;
    }
    auto scores = align_read__find_res(chunk_i, sketch_read, sketch_read_reverse);
    delete[] sketch_read;
    delete[] sketch_read_reverse;

    int read_begin = stoi(strstr(read.get_name_c(), "startpos=") + 9);
    int read_len = stoi(strstr(read.get_name_c(), "length=") + 7);
    int correct_chunk_index = -1, correct_score = -1, correct_rank = -1;
    auto results = unordered_set<int>();
    auto results_score = vector<int>();
    int chunk_ratio = CHUNK_SIZES[CHUNK_SIZES_LEN - 1] / CHUNK_SIZES[chunk_i];
    while (!scores.empty() && results.size() < MAX_ALT_MATCHS) {
        auto te = scores.back();
        scores.pop_back();
        int max_simm = te.first, max_simm_i = te.second;
        if (results.count(max_simm_i + 1) || results.count(max_simm_i - 1))
            continue;

        results.insert(max_simm_i);
        results_score.push_back(max_simm);
        logger->debugl2("our result: read:%04d num:%d score:%d index:%d", read_i, results.size(), max_simm, max_simm_i);
        for (int i = 0; i < genome_parts_starts[chunk_i].size() - 1; ++i)
            if (max_simm_i >= genome_parts_starts[chunk_i][i] && max_simm_i < genome_parts_starts[chunk_i][i + 1]) {
                max_simm_i -= genome_parts_starts[chunk_i][i];
                total_results[read_i].push_back(max_simm_i / chunk_ratio + genome_parts_starts[CHUNK_SIZES_LEN - 1][i]);
                break;
            }

        // check is correct or not
        if (correct_score == -1) {
            if ((max_simm_i - 1) * int(chunk_size) / 2 <= read_begin &&
                read_begin + read_len <= (max_simm_i + 3) * chunk_size / 2) {
                correct_rank = static_cast<int>(results.size());
                correct_chunk_index = max_simm_i;
                correct_score = max_simm;
            }
        }
    }
    logger->debugl2("Read name: %s", read.get_name_c());
    int alt_score = correct_rank == 1 ? (results_score.size() > 1 ? results_score[1] : -1) : results_score[0];
    logger->debug("our result: read:%04d\t%04d\t%d\t%03d\t%03d\t%03d\t%d\t%d", read_i, correct_chunk_index,
                  correct_rank, correct_score, alt_score, read_len, int(log2(CHUNK_SIZES[chunk_i])),
                  CHUNK_SIZES[chunk_i]);
    tot_res += results.size();
    return correct_chunk_index > -1;
}


unsigned int sketchize_chunk_i;
vector<Sequence> genome_chunks;
pthread_mutex_t sketchize_mutex[BIG_PRIME_NUMBER + 1];

int make_ref_sketch__skethize(const int i) {
    auto chunk_i = sketchize_chunk_i;
    vector<int> *ref_hash_sketch = ref_hash_sketchs[chunk_i];
    int *sketch = similarity_claz->get_sketch(genome_chunks[i], chunk_i, GINGLE_LENGTH, 0);
    if (sketch[0] == -1)
        return 0;

    // if (i == 90864) { printf(sketch_str(sketch, SKETCH_SIZE).c_str()); printf("\n"); char buffer[32001];
    // strncpy(buffer, genome_chunks[i].seq_str, 32000); buffer[32000] = 0; printf(buffer); printf("\n"); }

    pthread_mutex_lock(&sketchize_mutex[sketch[0]]);
    ref_hash_sketch[sketch[0]].push_back(i);
    pthread_mutex_unlock(&sketchize_mutex[sketch[0]]);
    for (int j = 1; j < SKETCH_SIZE; ++j)
        if (sketch[j] != sketch[j - 1]) {
            pthread_mutex_lock(&sketchize_mutex[sketch[j]]);
            ref_hash_sketch[sketch[j]].push_back(i);
            pthread_mutex_unlock(&sketchize_mutex[sketch[j]]);
        }

    delete[] sketch;
    if (i % 1000 == 0)
        logger->debug("loading reference in chunk size 2^%d: %d records", log_chunk, i);
    return 0;
}


auto make_ref_sketch(const char *const ref_file_name, const BasketMinHash &similarity_class, const unsigned int chunk_i,
                     const int gingle_length, const int gap_length, const char *const index_base_file_name = nullptr,
                     const bool write_index = WRITE_INDEX, bool read_index = READ_INDEX) {
    auto chunk_size = CHUNK_SIZES[chunk_i];

    sketchize_chunk_i = chunk_i;
    log_chunk = static_cast<unsigned int>(log2(chunk_size));

    auto index_file_name =
            Logger::formatString("%s_%d_%d_%d_%d_%d.gin",
                                 index_base_file_name == nullptr ? ref_file_name : index_base_file_name, SKETCH_SIZE,
                                 gingle_length, gap_length, LOG_MAX_BASENUMBER, log_chunk);
    add_time();
    if (read_index)
        try {
            ref_hash_sketchs[chunk_i] = read_from_file(index_file_name.c_str());
            logger->info("chunk %d read completed", chunk_i);
        } catch (...) {
            read_index = false;
        }
    if (!read_index) {
        if (ref_genome.empty()) {
            ref_genome = read_from_file(ref_file_name, FASTA);
            add_time();
            logger->info("loaded reference: %d ms", last_time());

            // char buffer[1001]; strncpy(buffer, ref_genome[343].seq_str + 7106518, 1000); buffer[1001] = 0;
            // printf(buffer);printf("\n");
        }
        vector<int> *ref_hash_sketch = ref_hash_sketchs[chunk_i] = new vector<int>[BIG_PRIME_NUMBER + 2];

        tie(genome_chunks, genome_double_chunks, ref_hash_sketch[BIG_PRIME_NUMBER + 1]) =
                Sequence::chunkenize_big_sequence(ref_genome, chunk_size, chunk_i == CHUNK_SIZES_LEN - 1);
        if (chunk_i == CHUNK_SIZES_LEN - 1) {
            mkdir(CHUNKS_FOLDER_NAME, 0777);
            multiproc(THREADS_COUNT, create_aryana_indexes, static_cast<int>(genome_double_chunks.size()));
        }
        add_time();
        multiproc(THREADS_COUNT, make_ref_sketch__skethize, static_cast<int>(genome_chunks.size()));
    }
    // BIG_PRIME_NUMBER + 1 for genome_parts_starts and last one for ref_chunk_size[chunk_i]
    genome_parts_starts[chunk_i] = ref_hash_sketchs[chunk_i][BIG_PRIME_NUMBER + 1];
    add_time();

    if (!read_index && write_index)
        write_to_file(index_file_name.c_str(), ref_hash_sketchs[chunk_i], BIG_PRIME_NUMBER + 2);
    add_time();
    logger->info("loaded reference sketch/names in chunk size 2^%d: %sms: %d records", log_chunk, get_times_str(true),
                 genome_parts_starts[chunk_i].back());
    genome_chunks.clear();
    if (chunk_i == CHUNK_SIZES_LEN - 1 && !read_index)
        delete[] ref_genome[0].get_name_c();
}


int main(int argsc, char *argv[]) {
    auto *ref_file_name = const_cast<char *>(REF_FILE_NAME_DEFAULT);
    auto *reads_file_name = const_cast<char *>(READS_FILE_NAME_DEFAULT);
    int log_level = Logger::DEBUG;
    bool read_index = READ_INDEX, write_index = WRITE_INDEX;
    for (int i = 0; i < argsc; ++i) {
        char *key = argv[i];
        char *value = nullptr;
        for (int j = 0; argv[i][j] > 0; ++j)
            if (argv[i][j] == '=') {
                argv[i][j] = 0;
                value = argv[i] + j + 1;
                break;
            }
        if (!strcmp(key, "--ref"))
            ref_file_name = value;
        if (!strcmp(key, "--reads"))
            reads_file_name = value;
        if (!strcmp(key, "--alt-match-ratio"))
            alt_matchs_ratio = stod(value);
        if (!strcmp(key, "--log"))
            log_level = stoi(value);
        if (!strcmp(key, "--read-index"))
            read_index = strcmp(value, "false") != 0;
        if (!strcmp(key, "--write-index"))
            write_index = strcmp(value, "false") != 0;
    }
    ref_file_base_name = strstr(ref_file_name, "/") + 1;
    logger = new Logger(log_level);

    auto config_str = Logger::formatString(
            "shingle length:%d, gap_length:%d, base_number:%d, chunk begin:%d, chunk end:%d", GINGLE_LENGTH, GAP_LENGTH,
            LOG_MAX_BASENUMBER, int(log2(CHUNK_SIZES[0])), int(log2(CHUNK_SIZES[CHUNK_SIZES_LEN - 1])));
    const char *config = config_str.c_str();
    logger->info("begin with config: %s", config);
    similarity_claz = new BasketMinHash(1, zigma_hash, BIG_PRIME_NUMBER);
    add_time();

    for (unsigned int chunk_i = 0; chunk_i < CHUNK_SIZES_LEN; chunk_i++)
        make_ref_sketch(ref_file_name, *similarity_claz, chunk_i, GINGLE_LENGTH, 0, nullptr, write_index, read_index);
    add_time();

    reads = read_from_file(reads_file_name, FASTQ);
    add_time();
    logger->info("load reads time: %d ms", last_time());

    for (int i = 0; i < reads.size(); ++i)
        total_results[i] = vector<int>();
    int corrects = multiproc(THREADS_COUNT, align_read, static_cast<int>(reads.size()));
    add_time();
    logger->info("total results:%d", tot_res);

    write_results(ref_file_base_name, reads_file_name);
    logger->info("correct reads for config(%s): %d\nmaking sam files", config, corrects);
    run_aryana(ref_file_base_name, reads_file_name, reads, total_results);
    logger->info("total times %s", get_times_str(true));
    delete[] (reads[0].get_name_c() - 1);
    return 0;
}