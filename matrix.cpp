#include "matrix.hpp"
#include <iostream>
#include <stdexcept>


Matrix::Matrix(){
    matrix_ = nullptr; //there is no matrix

    //setting dimensions to 0
    dimensions_.rows_ = 0; 
    dimensions_.cols_ = 0;
}

Matrix::Matrix(float **matrix, int rows_init, int cols_init){

    //setting the appropriate dimensions
    dimensions_.rows_ = rows_init;
    dimensions_.cols_ = cols_init;

    //dynamically allocate the rows
    matrix_ = new float* [rows_init];
    

    for(int row{0}; row < rows_init; row++){

        //dynamically allocate each row
        matrix_[row] = new float[cols_init];

        //set each element in the rows to be equal
        for(int col{0}; col < cols_init; col++){
            matrix_[row][col] = matrix[row][col];
        }
    }
}
Matrix::Matrix(int rows_init, int cols_init){

    //setting the appropriate dimensions
    dimensions_.rows_ = rows_init;
    dimensions_.cols_ = cols_init;

    //dynamically allocate the rows
    matrix_ = new float* [rows_init];
    

    for(int row{0}; row < rows_init; row++){

        //dynamically allocate each row
        matrix_[row] = new float[cols_init];

        //set each element in the rows to be equal
        for(int col{0}; col < cols_init; col++){
            matrix_[row][col] = 0;
        }
    }
}
Matrix::Matrix(float *matrix, int rows_init, int cols_init){

    //setting the appropriate dimensions
    dimensions_.rows_ = rows_init;
    dimensions_.cols_ = cols_init;

    //dynamically allocate the rows
    matrix_ = new float* [rows_init];
    

    for(int row{0}; row < rows_init; row++){

        //dynamically allocate each row
        matrix_[row] = new float[cols_init];

        //set each element in the rows to be equal
        for(int col{0}; col < cols_init; col++){
            matrix_[row][col] = matrix[(row*cols_init)+col];
        }
    }
}
Matrix *Matrix::new_matrix(){

    //create a temporary array to store the user input
    int m;
    int n;
    std::cout<<"Create an mxn matrix...\n";
    std::cout<<"m: ";
    std::cin>>m;
    std::cout<<"n: ";
    std::cin>>n;
    std::cout<<"Creating an "<<m<<"x"<<n<<" matrix..."<<std::endl;

    float **arr = new float*[m];
    for(int i{0}; i<m; i++){
        arr[i] = new float[n];
        for(int j{0}; j<n;j++){
            std::cout<<i<<","<<j<<"th entry: ";
            std::cin>>arr[i][j];
        }
    }

    //create matrix
    Matrix* matrix  = new Matrix(arr, m, n);  //WHY DOES THIS NEED TO BE DYNAMICALLY ALLOCATED
    // Matrix matrix(arr, m, n);

    //delete array
    for(int i{0}; i<m; i++){
        for(int j{0}; j<n; j++){
            arr[i][j] = 0;
        }
        delete[] arr[i];
        arr[i] = nullptr;
    }
    delete arr;
    arr = nullptr;

    return matrix;

    // return &matrix; //could i do this or is it a local variable and will go out of scope?
}


Matrix::~Matrix(){
    //loop through each row, set each element and then delete the row.
    for(int row{0}; row < dimensions_.rows_; row++){

        //set all elements to zero
        for(int col{0}; col < dimensions_.cols_; col++){
            matrix_[row][col] = 0;
        }
        //delete each row
        delete[] matrix_[row];
        matrix_[row] = nullptr;
    }
    //delete the arr of rows
    delete[] matrix_;
    matrix_ = nullptr;
}

void Matrix::print(){
    for(int row{0}; row < dimensions_.rows_; row++){
        for(int col{0}; col < dimensions_.cols_; col++){
            std::cout<<matrix_[row][col]<<" ";
        }
        std::cout<<std::endl;
    }
}


Matrix Matrix::operator*(Matrix &m){ //we need to take a reference in order to use the indexing operator
    if(!is_mult_allowed(*this, m)){
        throw std::domain_error{
            "Dimensions do not agree for multiplication: m1: "+std::to_string(dimensions_.rows_)+"x"+std::to_string(dimensions_.cols_)+" m2: "+ std::to_string(m.dimensions_.rows_)+"x"+std::to_string(m.dimensions_.cols_)
        };
    }
    float *arr = new float[dimensions_.rows_ * m.dimensions_.cols_];

    for(int row_1{0}; row_1<dimensions_.rows_; row_1++){
        for(int col_2{0}; col_2<m.dimensions_.cols_; col_2++){
            arr[row_1*m.dimensions_.cols_ + col_2] = vrow(row_1) % (m[col_2]); //is equal to the dot product between the vector corresponding to the row of the first matrix and the column of the second.
        }
    }

    Matrix *matrix = new Matrix(arr, dimensions_.rows_, m.dimensions_.cols_);
    delete[] arr;
    arr = nullptr;
    return *matrix;
}

Matrix Matrix::operator[](int col){
    float *arr = new float[dimensions_.rows_];
    for(int i{0}; i<dimensions_.rows_; i++){
        arr[i] = matrix_[i][col];
    }
    Matrix *matrix = new Matrix(arr, dimensions_.rows_, 1);
    delete[] arr;
    arr = nullptr;
    return *matrix;
}



Matrix Matrix::operator*(float k){
    
    for(int i{0}; i<dimensions_.rows_; i++){
        for(int j{0}; j<dimensions_.cols_;j++){
            matrix_[i][j] *= k;
        }
    }
    return *this; //should i do this?? pun not intended
   
}

float Matrix::operator%(Matrix const &m){
    if(!is_vector(*this) || !is_vector(m)){
        throw std::domain_error{
            "Not a vector - dot product cannot be computed."
        };
    }
    if(!is_same_dim(*this, m)){
        throw std::domain_error{
            "Dimensions do not agree for dot product."
        };
    }
    
    float result{0};
    for(int i{0}; i<dimensions_.rows_; i++){
        result += matrix_[i][0] * m.matrix_[i][0];
    }
    return result;
}

Matrix Matrix::vrow(int row){
    float *arr = new float[dimensions_.cols_];
    for(int i{0}; i<dimensions_.cols_; i++){
        arr[i] = matrix_[row][i];
    }
    Matrix *matrix = new Matrix(arr, dimensions_.cols_, 1);
    delete[] arr;
    arr = nullptr;
    return *matrix;
}

Matrix Matrix::row(int row){
    float *arr = new float[dimensions_.cols_];

    for(int i{0}; i<dimensions_.cols_; i++){
        arr[i] = matrix_[row][i];
    }
    
    Matrix *matrix = new Matrix(arr, 1, dimensions_.cols_);
    delete[] arr;
    arr = nullptr;
    return *matrix;}


Matrix Matrix::operator+(Matrix const &m){
    if(!is_same_dim(*this, m)){
        throw std::domain_error{
            "Dimensions do not agree for addition: m1: "+std::to_string(dimensions_.rows_)+"x"+std::to_string(dimensions_.cols_)+" m2: "+ std::to_string(m.dimensions_.rows_)+"x"+std::to_string(m.dimensions_.cols_)
        };
    }
    for(int i{0}; i<dimensions_.rows_; i++){
        for(int j{0}; j<dimensions_.cols_;j++){
            matrix_[i][j] += m.matrix_[i][j];
        }
    }
    return *this; //should i do this?? pun not intended
}
Matrix Matrix::operator-(Matrix const &m){
    if(!is_same_dim(*this, m)){
        throw std::domain_error{
            "Dimensions do not agree for subtration: m1: "+std::to_string(dimensions_.rows_)+"x"+std::to_string(dimensions_.cols_)+" m2: "+ std::to_string(m.dimensions_.rows_)+"x"+std::to_string(m.dimensions_.cols_)
        };
    }
    for(int i{0}; i<dimensions_.rows_; i++){
        for(int j{0}; j<dimensions_.cols_;j++){
            matrix_[i][j] -= m.matrix_[i][j];
        }
    }
    return *this; //should i do this?? pun not intended
}

Matrix Matrix::operator|(Matrix const &m){
    return rref(m);
}
Matrix Matrix::eye(int n){
    Matrix *matrix = new Matrix(n, n);
    for(int i{0}; i < n; i++){
        matrix->matrix_[i][i] = 1;
    }
    return *matrix;
    
}

bool Matrix::is_same_dim(Matrix const &m1, Matrix const &m2){
    return (m1.dimensions_.rows_ == m2.dimensions_.rows_ && m1.dimensions_.cols_ == m2.dimensions_.cols_) ? true : false;
}
bool Matrix::is_square(Matrix const &m){
    return (m.dimensions_.rows_ == m.dimensions_.cols_) ? true : false;
}
bool Matrix::is_mult_allowed(Matrix const &m1, Matrix const &m2){
    return (m1.dimensions_.cols_ == m2.dimensions_.rows_) ? true : false;
}
int Matrix::is_vector(Matrix const &m1){
    return (m1.dimensions_.cols_ != 1) ? -1 : m1.dimensions_.rows_;
}

Matrix Matrix::operator/(Matrix const &m){
    
}
//more advanced operations
Matrix Matrix::ref(){

/*Algorithm for computing RREF

1. divide R1 by the first entry
2. subtract prev row * first entry of next row from next row.
3. divide
*/

//NOTE: Make assignmnent operator
Matrix matrix_ref(matrix_, dimensions_.rows_, dimensions_.cols_);


for(int row{0};  row < dimensions_.rows_-1; row++){
    
    float first_entry = matrix_ref.row(row).matrix_[0][row]; 
    // std::cout<<"Row: "<<row<<" first_entry: "<<first_entry<<std::endl;
    //divide row by first_entry

    // for(int i{row}; i < dimensions_.cols_; i++){
    //     matrix_ref.matrix_[row][i] *= (1/first_entry);
    // }

    for(int k{row + 1}; k < dimensions_.rows_; k++){
        float row_fe = matrix_ref.row(k).matrix_[0][row];
        // std::cout<<" next row: "<<k<<" row_fe: "<<row_fe<<std::endl;
        //subtract prev row * first entry next row from next row

        for(int i{row}; i < dimensions_.cols_; i++){
            // std::cout<<"   "<<matrix_ref.matrix_[k][i]<<" - "<< matrix_ref.matrix_[row][i] * row_fe<< " = "<<matrix_ref.matrix_[k][i] - matrix_ref.matrix_[row][i] * (row_fe/first_entry)<<std::endl;
            matrix_ref.matrix_[k][i] -= matrix_ref.matrix_[row][i] *  ((first_entry == 0) ? 0 : (row_fe/first_entry));
            
        }
    }

    // std::cout<<"\n";
    // matrix_ref.print();
    // std::cout<<std::endl;
}

return matrix_ref;

}

Matrix Matrix::rref(){
    Matrix matrix_rref(ref().matrix_, dimensions_.rows_, dimensions_.cols_);
    //start at bottom and work our way up
    // matrix_rref.print();
    // std::cout<<std::endl;
    int entry{0};
    for(int row{std::min(dimensions_.cols_, dimensions_.rows_ )-1};  row > 0; row--, entry++){
        
        float last_entry = matrix_rref.row(row).matrix_[0][std::min(dimensions_.cols_, dimensions_.rows_ ) - 1 - entry]; 
 
        for(int k{row - 1}; k >= 0; k--){
            float row_le = matrix_rref.row(k).matrix_[0][std::min(dimensions_.cols_, dimensions_.rows_ )- 1 - entry];
            // std::cout<<" next row: "<<k<<" row_fe: "<<row_le<<std::endl;
            // subtract prev row * first entry next row from next row

            for(int i{0}; i < dimensions_.cols_; i++){
                // std::cout<<"   "<<matrix_rref.matrix_[k][i]<<" - "<< matrix_rref.matrix_[row][i] * ((last_entry == 0) ? 0 : (row_le/last_entry))<< " = "<<matrix_rref.matrix_[k][i] - matrix_rref.matrix_[row][i] *  ((last_entry == 0) ? 0 : (row_le/last_entry))<<std::endl;
                matrix_rref.matrix_[k][i] -= matrix_rref.matrix_[row][i] * ((last_entry == 0) ? 0 : (row_le/last_entry));
                
            }
        }
        

        // std::cout<<"\n";
        // matrix_rref.print();
        // std::cout<<std::endl;
    }

    for(int i{0}; i < std::min(dimensions_.cols_, dimensions_.rows_); i++){
        float leading_entry = 1;
        for(int j{0}; j<matrix_rref.dimensions_.cols_; j++){
            if(matrix_rref.matrix_[i][j] != 0){
                leading_entry = matrix_rref.matrix_[i][j];
                break;
            }
        }
        for(int j{0}; j < matrix_rref.dimensions_.cols_;j++){
            // std::cout<<matrix_rref.matrix_[i][j] <<"/"<< leading_entry<<" = "<<matrix_rref.matrix_[i][j] / leading_entry<<std::endl;
            matrix_rref.matrix_[i][j] /= leading_entry;
            
            if(matrix_rref.matrix_[i][j] == -0){
                matrix_rref.matrix_[i][j] = 0;
            }

        }
        

        // std::cout<<std::endl;
    }
    return matrix_rref;

}

Matrix Matrix::ref(Matrix const &m){

    if(dimensions_.rows_ != m.dimensions_.rows_){
        throw std::domain_error{
            "Nbr of rows does not agree for solving."
        };
    }
    /*Algorithm for computing RREF

    1. divide R1 by the first entry
    2. subtract prev row * first entry of next row from next row.
    3. divide
    */

    //NOTE: Make assignmnent operator
    Matrix matrix_ref(matrix_, dimensions_.rows_, dimensions_.cols_);
    Matrix secondary_matrix(m.matrix_, m.dimensions_.rows_, m.dimensions_.cols_);


    for(int row{0};  row < dimensions_.rows_-1; row++){
        
        float first_entry = matrix_ref.row(row).matrix_[0][row]; 
        // std::cout<<"Row: "<<row<<" first_entry: "<<first_entry<<std::endl;
        //divide row by first_entry

        // for(int i{row}; i < dimensions_.cols_; i++){
        //     matrix_ref.matrix_[row][i] *= (1/first_entry);
        // }

        for(int k{row + 1}; k < dimensions_.rows_; k++){
            float row_fe = matrix_ref.row(k).matrix_[0][row];
            // std::cout<<" next row: "<<k<<" row_fe: "<<row_fe<<std::endl;
            //subtract prev row * first entry next row from next row

            for(int i{row}; i < dimensions_.cols_; i++){
                // std::cout<<"   "<<matrix_ref.matrix_[k][i]<<" - "<< matrix_ref.matrix_[row][i] * row_fe<< " = "<<matrix_ref.matrix_[k][i] - matrix_ref.matrix_[row][i] * (row_fe/first_entry)<<std::endl;
                matrix_ref.matrix_[k][i] -= matrix_ref.matrix_[row][i] *  ((first_entry == 0) ? 0 : (row_fe/first_entry));
        
                
            }
            for(int i{0}; i<m.dimensions_.cols_;i++ ){
                secondary_matrix.matrix_[k][i] -= secondary_matrix.matrix_[row][i] *  ((first_entry == 0) ? 0 : (row_fe/first_entry));
            }
        }

        // std::cout<<"\n";
        // matrix_ref.print();
        // std::cout<<std::endl;
    }

    return secondary_matrix;

}

Matrix Matrix::rref(Matrix const &m){

    if(dimensions_.rows_ != m.dimensions_.rows_){
        throw std::domain_error{
            "Nbr of rows does not agree for solving."
        };
    }
    
    Matrix matrix_rref(ref().matrix_, dimensions_.rows_, dimensions_.cols_);
    Matrix secondary_matrix(ref(m).matrix_, m.dimensions_.rows_, m.dimensions_.cols_);

    //start at bottom and work our way up
    // matrix_rref.print();
    // std::cout<<std::endl;
    int entry{0};
    for(int row{std::min(dimensions_.cols_, dimensions_.rows_ )-1};  row > 0; row--, entry++){
        
        float last_entry = matrix_rref.row(row).matrix_[0][std::min(dimensions_.cols_, dimensions_.rows_ ) - 1 - entry]; 
 
        for(int k{row - 1}; k >= 0; k--){
            float row_le = matrix_rref.row(k).matrix_[0][std::min(dimensions_.cols_, dimensions_.rows_ )- 1 - entry];
            // std::cout<<" next row: "<<k<<" row_fe: "<<row_le<<std::endl;
            // subtract prev row * first entry next row from next row

            for(int i{0}; i < dimensions_.cols_; i++){
                // std::cout<<"   "<<matrix_rref.matrix_[k][i]<<" - "<< matrix_rref.matrix_[row][i] * ((last_entry == 0) ? 0 : (row_le/last_entry))<< " = "<<matrix_rref.matrix_[k][i] - matrix_rref.matrix_[row][i] *  ((last_entry == 0) ? 0 : (row_le/last_entry))<<std::endl;
                matrix_rref.matrix_[k][i] -= matrix_rref.matrix_[row][i] * ((last_entry == 0) ? 0 : (row_le/last_entry));
                
            }
            for(int i{0}; i < m.dimensions_.cols_; i++){
                secondary_matrix.matrix_[k][i] -= secondary_matrix.matrix_[row][i] *  ((last_entry == 0) ? 0 : (row_le/last_entry));
            }
            
        }
        

        // std::cout<<"\n";
        // matrix_rref.print();
        // std::cout<<std::endl;
    }

    for(int i{0}; i < std::min(dimensions_.cols_, dimensions_.rows_); i++){
        float leading_entry = 1;
        for(int j{0}; j<matrix_rref.dimensions_.cols_; j++){
            if(matrix_rref.matrix_[i][j] != 0){
                leading_entry = matrix_rref.matrix_[i][j];
                break;
            }
        }
        for(int j{0}; j < matrix_rref.dimensions_.cols_;j++){
            // std::cout<<matrix_rref.matrix_[i][j] <<"/"<< leading_entry<<" = "<<matrix_rref.matrix_[i][j] / leading_entry<<std::endl;
            matrix_rref.matrix_[i][j] /= leading_entry;
            
            if(matrix_rref.matrix_[i][j] == -0){
                matrix_rref.matrix_[i][j] = 0;
            }

        }
        for(int j{0}; j < secondary_matrix.dimensions_.cols_;j++){
            // std::cout<<matrix_rref.matrix_[i][j] <<"/"<< leading_entry<<" = "<<matrix_rref.matrix_[i][j] / leading_entry<<std::endl;
            secondary_matrix.matrix_[i][j] /= leading_entry;
            
            if(secondary_matrix.matrix_[i][j] == -0){
                secondary_matrix.matrix_[i][j] = 0;
            }

        }
        // std::cout<<std::endl;
    }
    return secondary_matrix;

}


Matrix Matrix::inv(){
    return *this|Matrix::eye(dimensions_.cols_);
}
Matrix Matrix::det(){

}

Matrix Matrix::transpose(){
    float *arr = new float[dimensions_.rows_ * dimensions_.cols_];

    for(int i{0}; i < dimensions_.rows_ ; i++){
        for(int j{0}; j < dimensions_.cols_; j++){
            arr[(j*dimensions_.rows_) + i] = matrix_[i][j];
        }
    }
    Matrix *matrix = new Matrix(arr, dimensions_.cols_, dimensions_.rows_);
    delete[] arr;
    arr = nullptr;
    return *matrix;
}