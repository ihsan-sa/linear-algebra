#ifndef MATRIX_H
#define MATRIX_H


struct Dimension{
    int rows_;
    int cols_;
};

class Matrix{
    float **matrix_;
    Dimension dimensions_;

public:

    //constructors/destructors
    Matrix(); //default constructor 
    Matrix(float **matrix, int m, int n); //constructor which takes a 2D array
    Matrix(float *matrix, int m, int n); //constructor which takes a 1D array
    ~Matrix(); //defuault destructor
    Matrix(int rows_init, int cols_init); //creates a matrix full of zeros

    static Matrix *new_matrix(); //take user input and create a matrix. Returns a pointer to the newly created matrix.


    //OPERATIONS
    void print(); //print the matrix
    
    //checking dimensions
    static bool is_same_dim(Matrix const &m1, Matrix const &m2);
    static bool is_square(Matrix const &m); 
    static bool is_mult_allowed(Matrix const &m1, Matrix const &m2); //returns true if multiplication is allowed (dimensions are appropriate)
    static int is_vector(Matrix const &m1); //returns the dimension of the vector, else returns -1



    //operations
    Matrix operator*(Matrix const &m); //matrix multiplication
    Matrix operator*(Matrix &m); //matrix multiplication -- this one uses the indexing operator, which is wht is cannot take a const reference
    Matrix operator*(float k); //scalar multiplication
    Matrix operator+(Matrix const &m); //adding matrices
    Matrix operator-(Matrix const &m); //subtracting matrices
    Matrix operator/(Matrix const &m); //solving a system
    Matrix operator[](int col); //returns a column
    float operator%(Matrix const &m); //dot product between two vectors
    Matrix operator|(Matrix const &m);

    Matrix vrow(int row); //returns a vector corresponding to the row
    Matrix row(int row); //returns row

    //more advanced operations
    Matrix transpose();
    Matrix ref();
    Matrix rref();
    Matrix ref(Matrix const &m); //returns the other matrix solved in parallel with first
    Matrix rref(Matrix const &m); //returns the other matrix solved in parallel with first
    Matrix inv();
    Matrix det();
    static Matrix eye(int n); //returns the identity matrix

};

#endif