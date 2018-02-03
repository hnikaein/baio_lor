#include "heap.h"

using namespace std;

Heap::Heap(int size) : size(size) {
    h.emplace_back(1 << 30, -1);
    index = new int[size];
}

Heap::~Heap() {
    h.clear();
    indexx.clear();
    delete[] index;
}


void Heap::inc_element(const int &element) {
    if (!indexx.count(element)) {
        h.emplace_back(0, element);
        index[element] = static_cast<int>(h.size() - 1);
        indexx.insert(element);
    }
    h[index[element]].first += 1;
    while (h[index[element]].first > h[index[element] / 2].first) {
        pair<int, int> te1 = move(h[index[element]]);
        h[index[element]] = move(h[index[element] / 2]);
        h[index[element] / 2] = move(te1);
        index[h[index[element]].second] = index[element];
        index[element] /= 2;
    }
}

pair<int, int> Heap::pop() {
    if (h.size() == 1)
        throw;
    auto result = move(h[1]);
    h[1] = move(h[h.size() - 1]);
    h.pop_back();
    int element = h[1].second;
    index[element] = 1;
    while ((index[element] * 2 + 1 < h.size() && h[index[element]].first < h[index[element] * 2 + 1].first) ||
           (index[element] * 2 < h.size() && h[index[element]].first < h[index[element] * 2].first)) {
        auto bigger = index[element] * 2;
        if (index[element] * 2 + 1 < h.size() && h[index[element] * 2].first < h[index[element] * 2 + 1].first)
            bigger = index[element] * 2 + 1;
        auto te1 = move(h[index[element]]);
        h[index[element]] = move(h[bigger]);
        h[bigger] = move(te1);
        index[h[index[element]].second] = index[element];
        index[element] = bigger;
    }
    return result;
}

bool Heap::has_element() {
    return h.size() > 1;
}

