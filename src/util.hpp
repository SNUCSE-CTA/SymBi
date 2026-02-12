#ifndef _UTIL_H_
#define _UTIL_H_

#include <algorithm>
#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#undef ASSERT
#define ASSERT(f) assert((f))
#define unlikely(x) __builtin_expect(!!(x), 0)
#define likely(x) __builtin_expect(!!(x), 1)
#define CACHE_LINE_SIZE 64

using namespace std;

#define MAX_BLOCK_SIZE 16777216

template <typename DataT>
class PODVector {
protected:
    size_t mSize;
    size_t mCapacity;
    DataT* mElements;
#ifdef MMAP
    bool dirty;
#endif

private:
    inline void reset_storage() {
        mSize = 0;
        mCapacity = 0;
        mElements = NULL;
    }

    inline DataT* allocate_block(size_t n) {
        return (n > 0) ? new DataT[n] : NULL;
    }

    inline void copy_prefix(DataT* dst, const DataT* src, size_t n) {
        if (n > 0) {
            memcpy((char*)dst, (char*)src, n * sizeof(DataT));
        }
    }

    inline void release_if_needed() {
        if (mCapacity > 0) {
            ASSERT(mElements != 0);
            delete[] mElements;
        }
    }

public:
    typedef DataT value_type;
    typedef DataT* pointer;
    typedef const DataT* const_pointer;
    typedef DataT* iterator;
    typedef const DataT* const_iterator;

    PODVector(const size_t iSize = 0) {
        mSize = iSize;
        mCapacity = iSize;
        mElements = allocate_block(mCapacity);
#ifdef MMAP
        dirty = false;
#endif
    }

    PODVector(const PODVector<DataT>& iRight) {
        mSize = iRight.mSize;
        mCapacity = iRight.mCapacity;
        mElements = allocate_block(mCapacity);
        copy_prefix(mElements, iRight.mElements, mSize);
#ifdef MMAP
        dirty = false;
#endif
    }

    inline PODVector<DataT>& operator=(const PODVector<DataT>& iRight) {
        release_if_needed();
        mCapacity = iRight.mCapacity;
        mSize = iRight.mSize;
        mElements = allocate_block(mCapacity);
        copy_prefix(mElements, iRight.mElements, mSize);
#ifdef MMAP
        dirty = false;
#endif
        return *this;
    }

    ~PODVector() {
#ifdef MMAP
        if (!dirty) {
#endif
            release_if_needed();
#ifdef MMAP
        }
#endif
    }

    void resize(const size_t iSize) {
        if (iSize > mCapacity) {
            if (mCapacity == 0) {
                mCapacity = iSize;
                mElements = allocate_block(mCapacity);
            } else {
                size_t next_capacity = std::max(iSize, mCapacity * 2);
                DataT* new_elements = allocate_block(next_capacity);
                copy_prefix(new_elements, mElements, mSize);
                delete[] mElements;
                mElements = new_elements;
                mCapacity = next_capacity;
            }
        }
        mSize = iSize;
        ASSERT(mSize <= mCapacity);
    }

    inline void resize_noerase(size_t new_size) {
        mSize = new_size;
    }

    inline void reserve_nocopy(size_t iCapacity) {
        if (iCapacity > mCapacity) {
            DataT* new_elements = allocate_block(iCapacity);
            if (mCapacity != 0) {
                delete[] mElements;
            }
            mElements = new_elements;
            mCapacity = iCapacity;
        }
    }

    void reserve(const size_t iCapacity) {
        if (iCapacity > mCapacity) {
            DataT* new_elements = allocate_block(iCapacity);
            copy_prefix(new_elements, mElements, mSize);
            if (mCapacity != 0) {
                delete[] mElements;
            }
            mElements = new_elements;
            mCapacity = iCapacity;
        }
    }

    void push_back(const DataT& iElem) {
        ASSERT(mSize <= mCapacity);
        if (mCapacity == 0) {
            ASSERT(mSize == 0);
            mCapacity = 2;
            mElements = allocate_block(mCapacity);
        } else if (mSize == mCapacity) {
            mCapacity *= 2;
            DataT* new_elements = allocate_block(mCapacity);
            memcpy((char*)new_elements, (char*)mElements, mSize * sizeof(DataT));
            delete[] mElements;
            mElements = new_elements;
        }
        mElements[mSize++] = iElem;
    }

    void insert(const DataT& item) {
        if (mSize == 0) {
            push_back(item);
            return;
        }

        if (mSize == mCapacity) {
            mCapacity *= 2;
            DataT* new_elements = allocate_block(mCapacity);
            memcpy((char*)new_elements, (char*)mElements, mSize * sizeof(DataT));
            delete[] mElements;
            mElements = new_elements;
        }

        int first = 0;
        int last = (int)mSize - 1;
        while (first <= last) {
            int mid = (first + last) / 2;
            if (item < mElements[mid]) {
                last = mid - 1;
            } else {
                first = mid + 1;
            }
        }

        size_t pos = (size_t)(last + 1);
        if (pos < mSize) {
            memmove(mElements + pos + 1, mElements + pos, sizeof(DataT) * (mSize - pos));
            mElements[pos] = item;
            ++mSize;
        } else {
            mElements[mSize++] = item;
        }
    }

    inline void pop_back() {
        ASSERT(mSize > 0);
        mSize -= 1;
    }

    inline DataT& operator[](const size_t iIndex) {
        ASSERT(iIndex < mSize);
        return mElements[iIndex];
    }

    inline const DataT& operator[](const size_t iIndex) const {
        ASSERT(iIndex < mSize);
        return mElements[iIndex];
    }

    bool operator==(const PODVector<DataT>& iRight) const {
        if (mSize != iRight.mSize) {
            return false;
        }
        for (size_t i = 0; i < mSize; ++i) {
            if (mElements[i] != iRight.mElements[i]) {
                return false;
            }
        }
        return true;
    }

    inline bool operator!=(const PODVector<DataT>& iRight) const {
        return !((*this) == iRight);
    }

    inline DataT& back() {
        ASSERT(mSize > 0);
        return mElements[mSize - 1];
    }

    inline const DataT& back() const {
        ASSERT(mSize > 0);
        return mElements[mSize - 1];
    }

    inline size_t size() const {
        return mSize;
    }

    inline size_t capacity() const {
        return mCapacity;
    }

    inline void clear() {
        mSize = 0;
    }

    inline void deallocate() {
        if (mCapacity > 0) {
            delete[] mElements;
        }
        reset_storage();
    }

    inline const_iterator begin() const {
        return mElements;
    }

    inline const_iterator end() const {
        return mElements + mSize;
    }

    inline iterator begin() {
        return mElements;
    }

    inline iterator end() {
        return mElements + mSize;
    }

    inline DataT* array() {
        return (size() > 0) ? mElements : NULL;
    }

    inline const DataT* array() const {
        return (const_pointer)((size() > 0) ? mElements : NULL);
    }

    iterator erase(iterator it) {
        ASSERT(mElements <= it && it < mElements + mSize && 0 < mSize);
        *it = mElements[mSize - 1];
        --mSize;
        return it;
    }

    inline void sort(const size_t i, const size_t j) {
        if (j - i > 1) {
            std::sort(mElements + i, mElements + j);
        }
    }

    void fill_n(const size_t iN, const DataT iValue) {
        ASSERT(iN <= size());
        for (size_t i = 0; i < iN; ++i) {
            mElements[i] = iValue;
        }
    }

    void fill_0(const size_t iN) {
        memset((char*)mElements, 0, iN * sizeof(DataT));
    }

    inline void encode_n(const size_t iN, DataT* oArray) const {
        ASSERT(iN <= size());
        memcpy((void*)oArray, (void*)mElements, iN * sizeof(DataT));
    }

    inline void decode_n(const DataT* iArray, const size_t iN) {
        ASSERT(iN <= size());
        memcpy((void*)mElements, (void*)iArray, iN * sizeof(DataT));
    }

#ifdef MMAP
    inline void decode_n_ptr(const DataT* iArray, const size_t iN) {
        delete mElements;
        mElements = (DataT*)iArray;
        mSize = iN;
        mCapacity = iN;
        dirty = true;
    }
#endif

    inline iterator find(const DataT& iValue, size_t bind, size_t eind) {
        iterator result = mElements + eind;
        while (bind < eind) {
            size_t mid = (bind + eind) / 2;
            if (mElements[mid] < iValue) {
                bind = mid + 1;
            } else if (iValue < mElements[mid]) {
                eind = mid;
            } else {
                return mElements + mid;
            }
        }
        return result;
    }

    inline const_iterator find(const DataT& iValue, size_t bind, size_t eind) const {
        const_iterator result = mElements + eind;
        while (bind < eind) {
            size_t mid = (bind + eind) / 2;
            if (mElements[mid] < iValue) {
                bind = mid + 1;
            } else if (iValue < mElements[mid]) {
                eind = mid;
            } else {
                return mElements + mid;
            }
        }
        return result;
    }

    inline bool contains(const DataT& iValue, size_t bind, size_t eind) const {
        while (bind < eind) {
            size_t mid = (bind + eind) / 2;
            if (mElements[mid] < iValue) {
                bind = mid + 1;
            } else if (iValue < mElements[mid]) {
                eind = mid;
            } else {
                return true;
            }
        }
        return false;
    }
};

char* inFileBuffer;
FILE* inFile;
long long leftFileSize;

long long getFileSize() {
    fseek(inFile, 0, SEEK_END);
    long long fileSize = ftell(inFile);
    fseek(inFile, 0, SEEK_SET);
    return fileSize;
}

void moveBufferForward(char** buf) {
    ++(*buf);
    if (leftFileSize == 0) {
        return;
    }

    if (**buf == 0) {
        long long readSize = MAX_BLOCK_SIZE;
        if (leftFileSize < readSize) {
            readSize = leftFileSize;
        }

        fread(inFileBuffer, sizeof(char), readSize, inFile);
        for (long long i = readSize; i <= MAX_BLOCK_SIZE; ++i) {
            inFileBuffer[i] = 0;
        }
        *buf = inFileBuffer;
        leftFileSize -= readSize;
    }
}

char* initFile(char* graphPath) {
    inFile = fopen(graphPath, "rb");
    leftFileSize = getFileSize();

    inFileBuffer = new char[MAX_BLOCK_SIZE + 1];
    inFileBuffer[0] = 0;

    char* inFileBufferPtr = inFileBuffer - 1;
    moveBufferForward(&inFileBufferPtr);
    return inFileBufferPtr;
}

void clearFile() {
    delete[] inFileBuffer;
    fclose(inFile);
    leftFileSize = 0;
}

int parseInt(char** buf) {
    while (**buf && **buf <= 32) {
        moveBufferForward(buf);
    }

    bool neg = false;
    if (**buf == '-') {
        neg = true;
        moveBufferForward(buf);
    } else if (**buf == '+') {
        moveBufferForward(buf);
    }

    int ret = 0;
    while (**buf >= '0' && **buf <= '9') {
        ret = (ret << 3) + (ret << 1) + (**buf - '0');
        moveBufferForward(buf);
    }

    if (neg) {
        ret *= -1;
    }
    return ret;
}

char parseChar(char** buf) {
    while (**buf && **buf <= 32) {
        moveBufferForward(buf);
    }

    char ret = **buf;
    moveBufferForward(buf);
    return ret;
}

inline long long factorization(int x) {
    switch (x) {
    case 1:
        return 1LL;
    case 2:
        return 2LL;
    case 3:
        return 6LL;
    case 4:
        return 24LL;
    case 5:
        return 120LL;
    case 6:
        return 720LL;
    case 7:
        return 5040LL;
    case 8:
        return 40320LL;
    case 9:
        return 362880LL;
    case 10:
        return 3628800LL;
    case 11:
        return 39916800LL;
    case 12:
        return 479001600LL;
    case 13:
        return 6227020800LL;
    case 14:
        return 87178291200LL;
    case 15:
        return 1307674368000LL;
    case 16:
        return 20922789888000LL;
    case 17:
        return 355687428096000LL;
    case 18:
        return 6402373705728000LL;
    case 19:
        return 121645100408832000LL;
    case 20:
        return 2432902008176640000LL;
    }

    long long result = 2432902008176640000LL;
    for (int i = 21; i <= x; i++) {
        result *= i;
    }
    return result;
}

#endif
