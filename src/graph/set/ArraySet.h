//
// Created by 梦想家xixi on 2021/10/25.
//

#ifndef MUTECT2CPP_MASTER_ARRAYSET_H
#define MUTECT2CPP_MASTER_ARRAYSET_H

#include <vector>
#include <algorithm>
template<class T>
class ArraySet {
private:
    std::vector<T> arraySet;

public:
    ArraySet(ArraySet const & other) : arraySet(other.arraySet) {};
    ArraySet() = default;;
    std::pair<typename std::vector<T>::iterator, bool> insert(T t) {
        typename std::vector<T>::iterator iter = std::find(arraySet.begin(), arraySet.end(), t);
        if(iter != arraySet.end())
            return std::pair<typename std::vector<T>::iterator, bool>(iter, false);
        else {
            arraySet.emplace_back(t);
            return std::pair<typename std::vector<T>::iterator, bool>(arraySet.end()--, true);
        }
    }
    int size() const {return arraySet.size();};
    typename std::vector<T>::const_iterator find(T t) {
        return std::find(arraySet.begin(), arraySet.end(), t);
    }
    typename std::vector<T>::iterator begin() {
        return arraySet.begin();
    }
    typename std::vector<T>::const_iterator end() {
        return arraySet.end();
    }
    void erase(typename std::vector<T>::iterator iter) {
        arraySet.erase(iter);
    }
    void erase(T t) {
        typename std::vector<T>::iterator iter = std::find(arraySet.begin(), arraySet.end(), t);
        if(iter != arraySet.end())
            arraySet.erase(iter);
    }
    void clear() {arraySet.clear();}
    bool empty() {return arraySet.size() == 0;};

    T  operator[](int i) {return arraySet[i];}

    std::vector<T> & getArraySet() {return arraySet;}
};




#endif //MUTECT2CPP_MASTER_ARRAYSET_H
