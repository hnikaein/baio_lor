#include "multiproc.h"
#include <pthread.h>

using namespace std;

vector<int> ids;
int last_unruned, result_count;
pthread_mutex_t mutex1, mutex2;

int (*start_routine)(const int);

void *consumer_thread(void *end) {
    while (true) {
        int my_id = -1;
        pthread_mutex_lock(&mutex1);
        if (last_unruned < *(int *) end)
            my_id = ids[last_unruned++];
        pthread_mutex_unlock(&mutex1);
        if (my_id == -1)
            return nullptr;
        int result = start_routine(my_id);
        pthread_mutex_lock(&mutex2);
        result_count += result;
        pthread_mutex_unlock(&mutex2);
    }
}

int multiproc(int thread_count, int (*start_routine1)(const int), int end, int from) {
    vector<int> ids1(static_cast<unsigned long>(end - from));
    for (int i = from; i < end; ++i)
        ids1[i - from] = i;
    return multiproc(thread_count, start_routine1, ids1);
}

int multiproc(int thread_count, int (*start_routine1)(const int), vector<int> ids1) {
    result_count = 0;
    last_unruned = 0;
    if (thread_count > 1) {
        ids = ids1;
        pthread_t thread[thread_count];
        start_routine = start_routine1;
        auto end = static_cast<int>(ids.size());
        for (int i = 0; i < thread_count; ++i)
            pthread_create(thread + i, nullptr, consumer_thread, &end);
        for (int i = 0; i < thread_count; ++i)
            pthread_join(thread[i], nullptr);
    } else
        for (int i : ids1)
            result_count += start_routine1(i);
    return result_count;
}
