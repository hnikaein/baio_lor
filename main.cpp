#include "configs.h"
#include "hash_utils.h"
#include "basket_min_hash.h"
#include "utils/time_profile.h"
#include "utils/logger.h"
#include "bio_utils/bio_reader.h"
#include "utils/vector_writer.h"
#include "utils/heap.h"
#include <iostream>
#include <cmath>
#include <boost/format.hpp>
#include <cstring>
#include <omp.h>

using namespace std;

//vector<Sequence> ref_chunk_sequences[CHUNK_SIZES_LEN];
vector<Sequence> reads, ref_genome;
Logger *logger;
vector<int> *ref_hash_sketchs[CHUNK_SIZES_LEN];
double alt_matchs_ratio = ALT_MATCHS_RATIO_DEFAULT;
vector<int> genome_parts_starts[CHUNK_SIZES_LEN];

vector<pair<int, int>>
align_read__find_res2(const int chunk_i, const int *sketch_read, const int *sketch_read_reverse) {
    auto temp_result = vector<pair<int, int>>();
    int max_score = 0;
    for (auto sketch:{sketch_read, sketch_read_reverse}) {
        auto scores = Heap(genome_parts_starts[chunk_i].back());
        for (int i = 0; i < SKETCH_SIZE; ++i)
            for (int j:ref_hash_sketchs[chunk_i][sketch[i]])
                scores.inc_element(j);
        while (scores.has_element()) {
            auto p = scores.pop();
            if (p.first < alt_matchs_ratio * max_score)
                break;
            temp_result.push_back(p);
            if (p.first > max_score)
                max_score = p.first;
        }
    }
//    sort(temp_result.begin(), temp_result.end(),
//         [](const pair<int, int> &a, const pair<int, int> &b) { return a.first > b.first; });
//    int i;
//    for (i = 1; i < temp_result.size(); ++i)
//        if (temp_result[i].first < alt_matchs_ratio * temp_result[0].first)
//            break;
//    temp_result.erase(temp_result.begin() + i, temp_result.end());
//    return temp_result;
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

auto align_read(const int read_i, const BasketMinHash &similarity_claz, const int gingle_length, const int gap_length) {
    Sequence &read = reads[read_i];
    unsigned int chunk_i = 0;
    while (chunk_i < CHUNK_SIZES_LEN and CHUNK_SIZES[chunk_i] < read.size)
        chunk_i += 1;
    if (chunk_i == CHUNK_SIZES_LEN)
        chunk_i -= 1;
    const unsigned int chunk_size = CHUNK_SIZES[chunk_i];
    const string &read_name = read.get_name();
    int read_begin = stoi(read_name.substr(read_name.find("startpos=") + 9, read_name.find("_number")));
    int read_len = stoi(read_name.substr(read_name.find("length=") + 7, read_name.find("bp_")));
    int correct_index = -1, correct_score = -1;
    // TODO   read.write_to_file("tmp/%04dr.fastq" % read_i);

    add_time();
    auto sketch_read = similarity_claz.get_sketch(read, chunk_i, gingle_length, gap_length);
    auto sketch_read_reverse = similarity_claz.get_sketch(read.get_reversed(), chunk_i, gingle_length, gap_length);

    add_time();
//    auto scores = align_read__find_res(chunk_i, sketch_read, sketch_read_reverse);
    // XXX
    auto scores = align_read__find_res2(chunk_i, sketch_read, sketch_read_reverse);
    delete sketch_read;
    delete sketch_read_reverse;
    add_time();

    auto results = unordered_set<int>();
    while (!scores.empty()) {
        auto te = scores.back();
        scores.pop_back();
        int max_simm = te.first, max_simm_i = te.second;
        // if (results.count(max_simm_i + 1) || results.count(max_simm_i - 1))
        //    continue;
        // TODO ref_chunk_sequences[chunk_i][max_simm_i].write_to_file("tmp/%04dg.fasta" % read_i,append=True if i else False)

        results.insert(max_simm_i);

        logger->debugl2("our result: read:%04d num:%d score:%d index:%d", read_i, results.size(), max_simm, max_simm_i);

        for (int i = 0; i < genome_parts_starts[chunk_i].size() - 1; ++i)
            if (max_simm_i >= genome_parts_starts[chunk_i][i] && max_simm_i < genome_parts_starts[chunk_i][i + 1]) {
                max_simm_i -= genome_parts_starts[chunk_i][i];
                break;
            }
        if ((max_simm_i - 1) * int(chunk_size) / 2 <= read_begin &&
            read_begin + read_len <= (max_simm_i + 3) * chunk_size / 2) {
            correct_index = max_simm_i;
            correct_score = max_simm;
        }
    }
    add_time();
    // TODO run_aryana(read_i);

    add_time();
    logger->debugl2("Read name: %s", read.get_name_c());
    logger->debug("our result: read:%04d\t%04d\t%d\t%d\t%d\t%d\t%s", read_i, correct_index, correct_score, read_len,
                  int(log2(CHUNK_SIZES[chunk_i])), CHUNK_SIZES[chunk_i], get_times_str(true));
    tot_res += results.size();
    return correct_index > -1; // TODO return times

}

auto make_ref_sketch(const char *const ref_file_name, const BasketMinHash &similarity_class, const unsigned int chunk_i,
                     const int gingle_length, const int gap_length, const char *const index_base_file_name = nullptr,
                     const bool write_index = WRITE_INDEX, bool read_index = READ_INDEX) {
    auto chunk_size = CHUNK_SIZES[chunk_i];
    auto log_chunk = (int) log2(chunk_size);
    auto index_file_name = str(boost::format("%s_%d_%d_%d_%d_%d.gin") %
                               (index_base_file_name == nullptr ? ref_file_name : index_base_file_name) % SKETCH_SIZE %
                               gingle_length % gap_length % LOG_MAX_BASENUMBER % log_chunk);
    add_time();
    if (read_index)
        try {
            ref_hash_sketchs[chunk_i] = read_from_file(index_file_name.c_str());
            logger->info("chunk %d read completed", chunk_i);
        } catch (...) {
            read_index = false;
        }
    if (!read_index) {
        if (ref_genome.empty())
            ref_genome = read_from_file(ref_file_name, FASTA);
        add_time();
        logger->info("loaded reference: %d ms", last_time());

        vector<Sequence> genome_chunks = Sequence::chunkenize_big_sequence(ref_genome, chunk_size);
        vector<int> *ref_hash_sketch = ref_hash_sketchs[chunk_i] = new vector<int>[BIG_PRIME_NUMBER + 2];
        auto genome_chunks_size = static_cast<int>(genome_chunks.size());

        for (int i = 0; i < genome_chunks_size; ++i)
            if (!strncmp(genome_chunks[i].get_name().c_str(), "00000000", 8))
                ref_hash_sketch[BIG_PRIME_NUMBER + 1].push_back(i);
#pragma omp parallel for
        for (int i = 0; i < genome_chunks_size; ++i) {
            int *sketch = similarity_class.get_sketch(genome_chunks[i], chunk_i, gingle_length, gap_length);
#pragma omp critical
            {
                ref_hash_sketch[sketch[0]].push_back(i);
                for (int j = 1; j < SKETCH_SIZE; ++j)
                    if (sketch[j] != sketch[j - 1])
                        ref_hash_sketch[sketch[j]].push_back(i);
            }
            delete sketch;
            if ((i / THREADS_COUNT) % (1000 / THREADS_COUNT) == 0)
                logger->debug("loading reference in chunk size 2^%d: %d records", log_chunk, i);
        }
        ref_hash_sketch[BIG_PRIME_NUMBER + 1].push_back(static_cast<int>(genome_chunks_size));
    }
    // BIG_PRIME_NUMBER + 1 for genome_parts_starts and last one for ref_chunk_size[chunk_i]
    genome_parts_starts[chunk_i] = ref_hash_sketchs[chunk_i][BIG_PRIME_NUMBER + 1];
    add_time();

    if (!read_index && write_index)
        write_to_file(index_file_name.c_str(), ref_hash_sketchs[chunk_i], BIG_PRIME_NUMBER + 2);
    add_time();
    logger->info("loaded reference sketch/names in chunk size 2^%d: %sms: %d records", log_chunk, get_times_str(true),
                 genome_parts_starts[chunk_i].back());
}


int main(int argsc, char *argv[]) {
    auto *ref_file_name = const_cast<char *>(REF_FILE_NAME_DEFAULT);
    auto *reads_file_name = const_cast<char *>(READS_FILE_NAME_DEFAULT);
    int log_level = Logger::DEBUG;
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
    }
    logger = new Logger(log_level);
    auto config_str = str(boost::format(
            "shingle length:%d, gap_length:%d, base_number:%d, chunk begin:%d, chunk end:%d") %
                          GINGLE_LENGTH % GAP_LENGTH % LOG_MAX_BASENUMBER % log2(CHUNK_SIZES[0]) %
                          log2(CHUNK_SIZES[CHUNK_SIZES_LEN - 1]));
    const char *config = config_str.c_str();
    logger->info("begin with config: %s", config);
    auto basket_min_hash = BasketMinHash(1, zigma_hash, BIG_PRIME_NUMBER);
    add_time();
    for (unsigned int chunk_i = 0; chunk_i < CHUNK_SIZES_LEN; chunk_i++)
        make_ref_sketch(ref_file_name, basket_min_hash, chunk_i, GINGLE_LENGTH, 0);
    add_time();
    reads = read_from_file(reads_file_name, FASTQ);
    add_time();
    logger->info("load reads time: %d ms", last_time());
    int corrects = 0;
    auto reads_size = static_cast<int>(reads.size());
    for (int i = 0; i < reads_size; i++)
        corrects += align_read(i, basket_min_hash, GINGLE_LENGTH, GAP_LENGTH);
    add_time();
    logger->info("total results:%d", tot_res);
    logger->info("correct reads for config(%s): %d\nmaking sam files", config, corrects);
    logger->info("total times %s", get_times_str(true));
    return 0;
}