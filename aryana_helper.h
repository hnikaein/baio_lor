#include "aryana/aryana_args.h"

#ifndef C_LONG_READ_ALIGNER_ARYANA_HELPER_H
#define C_LONG_READ_ALIGNER_ARYANA_HELPER_H

extern "C" void bwa_aln_core2(aryana_args *args);

void run_aryana();

void create_aryana_index();

#endif //C_LONG_READ_ALIGNER_ARYANA_HELPER_H
