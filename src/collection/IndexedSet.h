//
// Created by 梦想家xixi on 2022/3/1.
//

#ifndef MUTECT2CPP_MASTER_INDEXEDSET_H
#define MUTECT2CPP_MASTER_INDEXEDSET_H

#include <vector>
#include <unordered_map>

template<typename E, typename _hash, typename _equal>
class IndexedSet{
private:
    std::vector<E> elements;
    std::unordered_map<E, int, _hash, _equal> indexByElement;

public:
    IndexedSet(int initialCapacity) {
        elements.reverse(initialCapacity);
        indexByElement.reserve(initialCapacity);
    }

    IndexedSet(const std::vector<E> & values) {
        int initialCapacity = values.size();
        elements.reverse(initialCapacity);
        indexByElement.reserve(initialCapacity);
        int nextIndex = 0;
        for(const E & value : values) {
            if(indexByElement.find(value) != indexByElement.end()) {
                continue;
            }
            indexByElement.insert({value, nextIndex++});
            elements.template emplace_back(value);
        }
    }

    std::vector<E> & asList() {
        return elements;
    }

    typename std::vector<E>::iterator iterator() {
        return elements.begin();
    }

    int size() const {
        return elements.size();
    }

    bool contains(const E & o) {
        return indexByElement.find(o) != indexByElement.end();
    }

    bool add(const E & o) {
        if(contains(o)) {
            return false;
        }
        int nextIndex = size();
        elements.template emplace_back(o);
        indexByElement.insert({o, nextIndex});
        return true;
    }

    bool remove(const E & o) {
        int index = indexOf(o);
        if(index == -1) {
            return false;
        }
        typename std::vector<E>::iterator viter = std::advance(elements.begin(), index);
        elements.erase(viter);
        indexByElement.erase(o);
        viter = std::advance(elements.begin(), index);
        int nextIndex = index;
        while(viter != elements.end()) {
            indexByElement.at(*viter) == nextIndex++;
            viter++;
        }
        return true;
    }

    int indexOf(const E & o) {
        typename std::unordered_map<E, int, _hash, _equal>::iterator iter = indexByElement.find(o);
        return iter != indexByElement.end() ? iter->second : -1;
    }

    void clear() {
        elements.clear();
        indexByElement.clear();
    }

    E get(const int index) {
        if(index >= size())
            throw std::invalid_argument("out of range");
        return elements[index];
    }
};

#endif //MUTECT2CPP_MASTER_INDEXEDSET_H
