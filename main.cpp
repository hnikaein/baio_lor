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
#include <sys/stat.h>
#include <getopt.h>
#include "aryana/main.h"
#include "aryana_helper.h"

using namespace std;

vector<Sequence> reads, ref_genome, genome_double_chunks;
vector<int> *ref_hash_sketchs[CHUNK_SIZES_LEN];
vector<int> genome_parts_starts[CHUNK_SIZES_LEN];
map<int, vector<int>> total_results;

BasketMinHash *similarity_claz;
unsigned int log_chunk;


auto *ref_file_name = const_cast<char *>(REF_FILE_NAME);
auto *reads_file_name = const_cast<char *>(READS_FILE_NAME);
auto *output_file_name = const_cast<char *>(OUTPUT_FILE_NAME);
char *ref_file_base_name;
int log_level = Logger::INFO, threads_count = THREADS_COUNT;
bool read_index = READ_INDEX, write_index = WRITE_INDEX, simlord_reads = SIMLORD_READS;
extern Logger *logger;
double alt_matchs_ratio = ALT_MATCHS_RATIO;

void read_args(int argc, char *argv[]) {
    static struct option long_options[] =
            {
                    {"threads",            required_argument, nullptr, 't'},
                    {"ref",                required_argument, nullptr, 'r'},
                    {"query",              required_argument, nullptr, 'q'},
                    {"output",             required_argument, nullptr, 'o'},
                    {"log",                required_argument, nullptr, 'l'},
                    {"alt-match-ratio",    required_argument, nullptr, 'm'},
                    {"no-read-index",      no_argument,       nullptr, 'n'},
                    {"no-write-index",     no_argument,       nullptr, 'w'},
                    {"simlord-reads",      no_argument,       nullptr, 's'},
                    {"mismatch-penalty",   no_argument,       nullptr, 'M'},
                    {"gap-open-penalty",   no_argument,       nullptr, 'O'},
                    {"gap-extend-penalty", no_argument,       nullptr, 'E'},
            };

    int option_index = 0, c;
    while ((c = getopt_long(argc, argv, "t:r:q:o:m:l:nwsM:O:E:", long_options, &option_index)) >= 0)
        switch (c) {
            case 't':
                threads_count = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'r':
                ref_file_name = strdup(optarg);
                break;
            case 'q':
                reads_file_name = strdup(optarg);
                break;
            case 'o':
                output_file_name = strdup(optarg);
                break;
            case 'm':
                alt_matchs_ratio = strtod(optarg, nullptr);
                break;
            case 'l':
                log_level = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'n':
                read_index = false;
                break;
            case 'w':
                write_index = false;
                break;
            case 's':
                simlord_reads = true;
                break;
            case 'M':
                MISMATH_PENALTY = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'O':
                GAP_OPEN_PENALTY = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'E':
                GAP_EXTEND_PENALTY = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            default:
                break;
        }
    ref_file_base_name = strstr(ref_file_name, "/");
    if (ref_file_base_name == nullptr)
        ref_file_base_name = ref_file_name;
    else
        ref_file_base_name++;
    logger = new Logger(log_level);
    optind = 1;
}

int create_aryana_indexes(const int part) {
    if (part % 100 == 0)
        logger->debugl2("create_aryana_indexes part\t\t%d", part);
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

    int read_begin = simlord_reads ? stoi(strstr(read.get_name_c(), "startpos=") + 9) : 0;
    int read_len = simlord_reads ? stoi(strstr(read.get_name_c(), "length=") + 7) : 0;
    int correct_chunk_index = -1, correct_score = -1, correct_rank = -1;
    auto results = unordered_set<int>();
    auto results_score = vector<int>();
    int max_chunk_size = CHUNK_SIZES[CHUNK_SIZES_LEN - 1];
    int chunk_ratio = max_chunk_size / CHUNK_SIZES[chunk_i];
    while (!scores.empty() && results.size() < MAX_ALT_MATCHS) {
        auto te = scores.back();
        scores.pop_back();
        int max_simm = te.first, max_simm_i = te.second;
        if (results.count(max_simm_i + 1) || results.count(max_simm_i - 1))
            continue;

        results.insert(max_simm_i);
        results_score.push_back(max_simm);
        logger->debugl2("our result: read:%04d num:%d score:%d index:%d", read_i, results.size(), max_simm, max_simm_i);
        int aryana_chunk = 0, aryana_chunk_local = 0;
        for (int i = 0; i < genome_parts_starts[chunk_i].size() - 1; ++i)
            if (max_simm_i >= genome_parts_starts[chunk_i][i] && max_simm_i < genome_parts_starts[chunk_i][i + 1]) {
                max_simm_i -= genome_parts_starts[chunk_i][i];
                aryana_chunk_local = max_simm_i / chunk_ratio;
                if (aryana_chunk_local ==
                    genome_parts_starts[CHUNK_SIZES_LEN - 1][i + 1] - genome_parts_starts[CHUNK_SIZES_LEN - 1][i])
                    aryana_chunk_local--;
                aryana_chunk = aryana_chunk_local + genome_parts_starts[CHUNK_SIZES_LEN - 1][i];
                total_results[read_i].push_back(aryana_chunk);
                break;
            }

        // check is correct or not
        if (simlord_reads)
            if (correct_score == -1) {
                if ((aryana_chunk_local - 1) * max_chunk_size / 2 <= read_begin &&
                    read_begin + read_len <= (aryana_chunk_local + 3) * max_chunk_size / 2) {
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
int make_ref_sketch_gingle_length, make_ref_sketch_gap_length;

int make_ref_sketch__skethize(const int i) {
    auto chunk_i = sketchize_chunk_i;
    vector<int> *ref_hash_sketch = ref_hash_sketchs[chunk_i];
    int *sketch = similarity_claz->get_sketch(genome_chunks[i], chunk_i, make_ref_sketch_gingle_length,
                                              make_ref_sketch_gap_length);
    if (sketch[0] == -1)
        return 0;

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
        }
        vector<int> *ref_hash_sketch = ref_hash_sketchs[chunk_i] = new vector<int>[BIG_PRIME_NUMBER + 2];

        logger->debugl2("before chunkenize");
        tie(genome_chunks, genome_double_chunks, ref_hash_sketch[BIG_PRIME_NUMBER + 1]) =
                Sequence::chunkenize_big_sequence(ref_genome, chunk_size, chunk_i == CHUNK_SIZES_LEN - 1);
        logger->debugl2("after chunkenize");
        if (chunk_i == CHUNK_SIZES_LEN - 1) {
            mkdir(CHUNKS_FOLDER_NAME, 0777);
            multiproc(threads_count, create_aryana_indexes, static_cast<int>(genome_double_chunks.size()));
        }
        add_time();
        make_ref_sketch_gingle_length = gingle_length;
        make_ref_sketch_gap_length = gap_length;
        multiproc(threads_count, make_ref_sketch__skethize, static_cast<int>(genome_chunks.size()));
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
        delete[] ref_genome[0].get_name_c(); // because we allocate all other names and seqs in bahave of this
}


int main(int argc, char *argv[]) {
    read_args(argc, argv);

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
    int corrects = multiproc(threads_count, align_read, static_cast<int>(reads.size()));
    add_time();
    logger->info("total results:%d", tot_res);
    logger->debug("correct reads for config(%s): %d\nmaking sam files", config, corrects);

    run_aryana(ref_file_base_name, output_file_name, reads, total_results, threads_count);
    delete[] (reads[0].get_name_c() - 1);
    add_time();
    logger->info("total times %s", get_times_str(true));
    return 0;
}
