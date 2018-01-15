#ifndef C_LONG_READ_ALIGNER_CONFIGS_H
#define C_LONG_READ_ALIGNER_CONFIGS_H

const unsigned int CHUNK_SIZES[] = {1 << 10, 1 << 11, 1 << 12, 1 << 13, 1 << 14, 1 << 15};
const unsigned int SKETCH_SIZES[] = {1 << 8, 1 << 8, 1 << 8, 1 << 9, 1 << 9, 1 << 9};

#define BIG_PRIME_NUMBER            14'868'587  // for 10
//#define BIG_PRIME_NUMBER            12'744'527  // for 9
#define SKETCH_SIZE                 SKETCH_SIZES[chunk_i]
#define LOG_MAX_BASENUMBER          5
#define GAP_LENGTH                  1
#define GINGLE_LENGTH               14
#define THREADS_COUNT               1
#define MAX_ALT_MATCHS              10
#define ALT_MATCHS_RATIO_DEFAULT    0.85
#define REF_FILE_NAME_DEFAULT       BASE_FILE_NAME ".fasta"
#define READS_FILE_NAME_DEFAULT     BASE_FILE_NAME ".reads.fastq"
#define BASE_FILE_NAME              "seqs/ecoli"

#define array_len(x)                ((sizeof(x)) / (sizeof((x)[0])))
#define CHUNK_SIZES_LEN             array_len(CHUNK_SIZES)

#define READ_INDEX                  true
#define WRITE_INDEX                 true


#endif //C_LONG_READ_ALIGNER_CONFIGS_H
