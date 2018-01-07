#ifndef C_LONG_READ_ALIGNER_CONFIGS_H
#define C_LONG_READ_ALIGNER_CONFIGS_H

const unsigned int CHUNK_SIZES[] = {1 << 10, 1 << 11, 1 << 12, 1 << 13, 1 << 14, 1 << 15}; // TODO easier initiliaze

#define BIG_PRIME_NUMBER            14'868'629
#define SKETCH_RATIO                (1 << 6)
#define SKETCH_SIZE                 CHUNK_SIZES[chunk_i]/SKETCH_RATIO
#define LOG_MAX_BASENUMBER          5
#define GAP_LENGTH                  1
#define GINGLE_LENGTH               14
#define THREADS_COUNT               1
#define MAX_ALT_MATCHS              10
#define ALT_MATCHS_RATIO            0.75
#define BASE_FILE_NAME              "seqs/ecoli"
#define REF_FILE_NAME               BASE_FILE_NAME ".fasta"
#define READS_FILE_NAME             BASE_FILE_NAME ".reads.fastq"
#define RESULT_FILE_NAME            BASE_FILE_NAME ".results"
#define STAG_SIZE                   200
#define SSKETCH_SIZE                50
#define SSKETCH_WINDOW              10
#define array_len(x)                ((sizeof(x)) / (sizeof((x)[0])))
#define CHUNK_SIZES_LEN             array_len(CHUNK_SIZES)


#endif //C_LONG_READ_ALIGNER_CONFIGS_H
