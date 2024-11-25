#ifndef SET_H
#define SET_H

#include "matrix.hpp"

class Set;

class Set{

    Matrix *vectors_;
    int n_vectors_;
public:
    Set();
    Set(Matrix *vectors, int n_vectors);
    Set(int n_vectors);
    ~Set();

    Set operator=(Set const &s);

    // void insert(Matrix const &v, int index);
    void append(Matrix const &v);
    void print();

friend Matrix;
};

#endif 