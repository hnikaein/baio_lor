#include <utility>
#include <vector>
#include <map>
#include <unordered_set>

#ifndef C_LONG_READ_ALIGNER_HEAP_H
#define C_LONG_READ_ALIGNER_HEAP_H

class Heap {
public:
    Heap(int size);

    ~Heap();

    void inc_element(const int &element, int value = 1);

    std::pair<int, int> pop();

    bool has_element();

    int get_value(const int &element);

    int top_value();

private:
    int size;
    std::vector<std::pair<int, int>> h;
    int *index;
    std::unordered_set<int> indexx;
};

#endif //C_LONG_READ_ALIGNER_HEAP_H
