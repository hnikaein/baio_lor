#include <vector>
#include <map>
#include "aryana/aryana_args.h"
#include "bio_utils/sequence.h"

#ifndef C_LONG_READ_ALIGNER_ARYANA_HELPER_H
#define C_LONG_READ_ALIGNER_ARYANA_HELPER_H

extern "C" void bwa_aln_core2(aryana_args *args);

void run_aryana(const char *ref_file_name, std::vector<Sequence> &reads, std::map<int, std::vector<int>> &results);

#endif //C_LONG_READ_ALIGNER_ARYANA_HELPER_H
