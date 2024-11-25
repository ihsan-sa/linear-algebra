#include "set.hpp"
#include <stdexcept>
#include <iostream>

Set::Set(){
    vectors_ = nullptr;
    n_vectors_ = 0;
}
Set::Set(Matrix *vectors, int n_vectors){
    vectors_ = new Matrix[n_vectors];
    for(int k{0}; k < n_vectors; k++){
        vectors_[k] = vectors[k];
    }
    n_vectors_ = n_vectors;
}
Set::Set(int n_vectors){
    vectors_ = new Matrix[n_vectors];
    n_vectors_ = n_vectors;

}
Set::~Set(){
    for(int k{0}; k < n_vectors_; k++){
        delete &vectors_[k];
    }
}

void Set::append(Matrix const &v){
    Set temp = *this;
    if(Matrix::is_vector(v) != vectors_[0].dimensions_.rows_){
        throw std::domain_error{
            "Err: Dimensions not appropriate"
        };
    }
    
    for(int k{0}; k < n_vectors_; k++){
        delete &vectors_[k];
    }

    vectors_ = new Matrix[n_vectors_ + 1];
    for(int k{0}; k < n_vectors_; k++){
        vectors_[k] = temp.vectors_[k];
    }
    vectors_[n_vectors_] = v;

}

Set Set::operator=(Set const &s){
    Set new_set();

}

void Set::print(){
    std::cout<<std::endl;
    for(int n{0}; n < vectors_[0].dimensions_.rows_; n++){
        for(int k{0}; k < n_vectors_; k++){
            std::cout<<vectors_[k].matrix_[n][0]<<"  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

