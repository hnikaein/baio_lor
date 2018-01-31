#ifndef C_LONG_READ_ALIGNER_MULTIPROC_H
#define C_LONG_READ_ALIGNER_MULTIPROC_H

#include <pthread.h>

void *consumer_thread(void *count);

int multiproc(int thread_count, int (*start_routine)(const int), int end, int from = 0);

#endif //C_LONG_READ_ALIGNER_MULTIPROC_H
