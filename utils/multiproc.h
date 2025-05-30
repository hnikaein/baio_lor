#include <vector>

#ifndef C_LONG_READ_ALIGNER_MULTIPROC_H
#define C_LONG_READ_ALIGNER_MULTIPROC_H

void *consumer_thread(void *count);

int multiproc(int thread_count, int (*start_routine)(const int), int end, int from = 0);

int multiproc(int thread_count, int (*start_routine)(const int), std::vector<int> ids);

#endif //C_LONG_READ_ALIGNER_MULTIPROC_H
