#ifndef C_LONG_READ_ALIGNER_CONFIGS_H
#define C_LONG_READ_ALIGNER_CONFIGS_H

const unsigned int CHUNK_SIZES[] = {1 << 9, 1 << 11, 1 << 13, 1 << 15};
const unsigned int SKETCH_SIZES[] = {1 << 7, 1 << 8, 1 << 9, 1 << 10};

extern int mismath_penalty, gap_open_penalty, gap_extend_penalty;


#define BIG_PRIME_NUMBER            14868587
#define SKETCH_SIZE                 SKETCH_SIZES[chunk_i]
#define LOG_MAX_BASENUMBER          5
#define GAP_LENGTH                  1
#define GINGLE_LENGTH               12
#define THREADS_COUNT               3
#define MAX_ALT_MATCHS              30
#define ALT_MATCHS_RATIO            0.70
#define REF_FILE_NAME               BASE_FILE_NAME ".fasta"
#define READS_FILE_NAME             BASE_FILE_NAME ".reads.fastq"
#define OUTPUT_FILE_NAME            READS_FILE_NAME ".sam"
#define BASE_FILE_NAME              "seqs/ecoli"
#define CHUNKS_FOLDER_NAME          "chunks"
#define SIMLORD_READS               false


#define array_len(x)                ((sizeof(x)) / (sizeof((x)[0])))
#define CHUNK_SIZES_LEN             array_len(CHUNK_SIZES)

#define READ_INDEX                  true
#define WRITE_INDEX                 true


#endif //C_LONG_READ_ALIGNER_CONFIGS_H
