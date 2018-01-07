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

using namespace std;

vector<int *> ref_chunk_sketchs[CHUNK_SIZES_LEN];
vector<Sequence> ref_chunk_sequences[CHUNK_SIZES_LEN], reads, ref_genome;
Logger logger(Logger::DEBUG);
vector<int> ref_hash_sketchs[CHUNK_SIZES_LEN][BIG_PRIME_NUMBER + 1];

int vector_equal_counts(const int *a, const int *b, int size) {
    int equality = 0;
    for (int i = 0; i < size; ++i)
        if (a[i] == b[i])
            equality += 1;
    return equality;
}

vector<pair<int, int>>
align_read__find_res2(const int chunk_i, const int *sketch_read, const int *sketch_read_reverse) {
    auto temp_result = vector<pair<int, int>>();
    int max_score = 0;
    for (auto sketch:{sketch_read, sketch_read_reverse}) {
        auto scores = Heap(static_cast<int>(ref_chunk_sketchs[chunk_i].size()));
        for (int i = 0; i < SKETCH_SIZE; ++i)
            for (int j:ref_hash_sketchs[chunk_i][sketch[i]])
                scores.inc_element(j);
        while (scores.has_element()) {
            auto p = scores.pop();
            if (p.first < ALT_MATCHS_RATIO * max_score)
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
//        if (temp_result[i].first < ALT_MATCHS_RATIO * temp_result[0].first)
//            break;
//    temp_result.erase(temp_result.begin() + i, temp_result.end());
//    return temp_result;
    sort(temp_result.begin(), temp_result.end(),
         [](const pair<int, int> &a, const pair<int, int> &b) { return a.first < b.first; });
    unsigned long temp_result_size = temp_result.size();
    int first_score = temp_result[temp_result_size - 1].first;
    int i;
    for (i = static_cast<int>(temp_result_size - 2); i >= 0; --i)
        if (temp_result[i].first < ALT_MATCHS_RATIO * first_score)
            break;
    auto result = vector<pair<int, int>>(temp_result_size - i - 1);
    for (int k = i + 1; k < temp_result_size; ++k)
        result[k - i - 1] = move(temp_result[k]);
    return result;

}

auto align_read__find_res(const int chunk_i, const int *sketch_read, const int *sketch_read_reverse) {
    vector<int *> &ref_sketchs = ref_chunk_sketchs[chunk_i];
    int max_sim_i[MAX_ALT_MATCHS] = {}, max_sim[MAX_ALT_MATCHS] = {};
    int last_match_i = 1;
    int count = 0;
    for (int i = 0; i < ref_sketchs.size(); i++) {
        for (int sketch_sim: {vector_equal_counts(ref_sketchs[i], sketch_read, SKETCH_SIZE),
                              vector_equal_counts(ref_sketchs[i], sketch_read_reverse, SKETCH_SIZE)}) {
            count += sketch_sim;
            if (sketch_sim <= max_sim[MAX_ALT_MATCHS - 1] || sketch_sim <= ALT_MATCHS_RATIO * max_sim[0])
                continue;
            for (int j = 0; j < last_match_i; j++) {
                if (sketch_sim >= max_sim[j]) {
                    if (i > max_sim_i[j] + 1) { // TODO should check over all of array
                        if (last_match_i < MAX_ALT_MATCHS)
                            last_match_i++;
                        for (int k = last_match_i - 1; k > j; k -= 1) {
                            max_sim_i[k] = max_sim_i[k - 1];
                            max_sim[k] = max_sim[k - 1];
                        }
                    }
                    max_sim_i[j] = i;
                    max_sim[j] = sketch_sim;
                    break;
                }
                if (i == max_sim_i[j] + 1)
                    break;
            }
        }
    }
    logger.info("%d", count);
    auto result = vector<pair<int, int>>();
    for (int i = last_match_i - 1; i >= 0; i--)
        if (max_sim[i] >= ALT_MATCHS_RATIO * max_sim[0])
            result.emplace_back(max_sim[i], max_sim_i[i]);

    return result;
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
        if (results.count(max_simm_i + 1) || results.count(max_simm_i - 1))
            continue;
        // TODO ref_chunk_sequences[chunk_i][max_simm_i].write_to_file("tmp/%04dg.fasta" % read_i,append=True if i else False)

        results.insert(max_simm_i);

        logger.debugl2("our result: read:%04d num:%d score:%d index:%d name:%s", read_i, results.size(), max_simm,
                       max_simm_i, ref_chunk_sequences[chunk_i][max_simm_i].get_name_c());

        if ((max_simm_i - 1) * int(chunk_size) / 2 <= read_begin &&
            read_begin + read_len <= (max_simm_i + 3) * chunk_size / 2) {
            correct_index = max_simm_i;
            correct_score = max_simm;
        }
    }
    add_time();
    // TODO run_aryana(read_i);

    add_time();
    logger.debugl2("Read name: %s", read.get_name_c());
    logger.debug("our result: read:%04d\t%04d\t%d\t%d\t%d\t%d\t%s", read_i, correct_index, correct_score, read_len,
                 int(log2(CHUNK_SIZES[chunk_i])), CHUNK_SIZES[chunk_i], get_times_str(true));
    tot_res += results.size();
    return correct_index > -1; // TODO return times

}

auto make_ref_sketch(const char *const ref_file_name, const BasketMinHash &similarity_class, const unsigned int chunk_i,
                     const int gingle_length, const int gap_length, const char *const index_base_file_name = nullptr,
                     const bool write_index = true, bool read_index = true) {
    auto chunk_size = CHUNK_SIZES[chunk_i];
    auto log_chunk = (int) log2(chunk_size);
    auto index_file_name = str(boost::format("%s_%d_%d_%d_%d_%d.gin") %
                               (index_base_file_name == nullptr ? ref_file_name : index_base_file_name) %
                               similarity_class.sketch_ratio % gingle_length % gap_length % LOG_MAX_BASENUMBER %
                               log_chunk);
    vector<Sequence> genome_chunks;
    add_time();
    tie(genome_chunks, ref_chunk_sequences[chunk_i]) = Sequence::chunkenize_big_sequence(ref_genome, chunk_size);
    add_time();
    if (read_index) {
        ref_chunk_sketchs[chunk_i] = read_from_file(index_file_name.c_str());
        if (ref_chunk_sketchs[chunk_i].empty())
            read_index = false;
    }
    vector<int *> &ref_sketchs = ref_chunk_sketchs[chunk_i];
    if (!read_index) {
        // TODO    p = Pool(THREADS_COUNT);
        for (auto seq:genome_chunks) {
            ref_sketchs.push_back(similarity_class.get_sketch(seq, chunk_i, gingle_length, gap_length));
            if ((ref_sketchs.size() / THREADS_COUNT) % (100 / THREADS_COUNT) == 0)
                logger.debug("loading reference in chunk size 2^%d: %d records", log_chunk, ref_sketchs.size());
        }
    }
    add_time();
    vector<int> *ref_hash_sketch = ref_hash_sketchs[chunk_i];
    for (int i = 0; i < ref_sketchs.size(); ++i)
        for (int j = 0; j < SKETCH_SIZE; ++j)
//            if (ref_sketchs[i][j] > BIG_PRIME_NUMBER)
//                cout << "salam";
//            else
            ref_hash_sketch[ref_sketchs[i][j]].push_back(i);
    add_time();

    if (!read_index && write_index)
        write_to_file(index_file_name.c_str(), ref_chunk_sketchs[chunk_i], SKETCH_SIZE);
    add_time();
    logger.info("loaded reference sketch/names in chunk size 2^%d: %sms: %d records",
                log_chunk, get_times_str(true), ref_chunk_sequences[chunk_i].size());
}


int main() {
    auto config_str = str(boost::format(
            "sketch ratio:%d, shingle length:%d, gap_length:%d, base_number:%d, chunk begin:%d, chunk end:%d") %
                          SKETCH_RATIO % GINGLE_LENGTH % GAP_LENGTH % LOG_MAX_BASENUMBER % log2(CHUNK_SIZES[0]) %
                          log2(CHUNK_SIZES[CHUNK_SIZES_LEN - 1]));
    const char *config = config_str.c_str();
    logger.info("begin with config: %s", config);
    auto basket_min_hash = BasketMinHash(SKETCH_RATIO, 1, zigma_hash, BIG_PRIME_NUMBER);
    add_time();
    ref_genome = read_from_file(REF_FILE_NAME, FASTA);
    add_time();
    logger.info("loaded reference: %d ms", last_time());
    for (unsigned int chunk_i = 0; chunk_i < CHUNK_SIZES_LEN; chunk_i++)
        make_ref_sketch(REF_FILE_NAME, basket_min_hash, chunk_i, GINGLE_LENGTH, 0);
    add_time();
    reads = read_from_file(READS_FILE_NAME, FASTQ);
    add_time();
    logger.info("load reads time: %d ms", last_time());
    int corrects = 0;
    auto reads_size = static_cast<int>(reads.size());
    for (int i = 0; i < reads_size; i++)
        corrects += align_read(i, basket_min_hash, GINGLE_LENGTH, GAP_LENGTH);
    add_time();
    logger.info("total results:%d", tot_res);
    logger.info("correct reads for config(%s): %d\nmaking sam files", config, corrects);
    logger.info("total times %s", get_times_str(true));
    return 0;
}