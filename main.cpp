#include "matrix.hpp"
#include <iostream>

float **make_arr(){
    int m;
    int n;
    std::cout<<"Please enter the number of rows: ";
    std::cin>>m;
    std::cout<<"Please enter the number of columns: ";
    std::cin>>n;
    std::cout<<"Creating an "<<m<<"x"<<n<<" matrix\nEnter the following float values:"<<std::endl;

    float **arr = new float*[m];
    for(int i{0}; i<m; i++){
        arr[i] = new float[n];
        for(int j{0}; j<n;j++){
            std::cout<<i<<"-"<<j<<"th entry: ";
            std::cin>>arr[i][j];
        }
    }
    return arr; //can i do this??
}
void delete_arr(float **&arr, int m, int n){
    for(int i{0}; i<m; i++){
        for(int j{0}; j<n; j++){
            arr[i][j] = 0;
        }
        delete[] arr[i];
        arr[i] = nullptr;
    }
    delete arr;
    arr = nullptr;
    
}


int main(){


//some other stuff
{
    // int arr1[3][3] = {
    //     {1,2,3}, 
    //     {4,5,6}, 
    //     {7,8,9}
    // };
    // int *arr2 = reinterpret_cast<int *>(arr1);
    // for(int i{0}; i<9; i++){
    //     std::cout<<arr2[i]<<"\n";
    // }
   
    // Matrix B(reinterpret_cast<float *>(float [3][3]{
    //     {3,2,1},
    //     {6,5,4}, 
    //     {9,8,7}
    // }), 3, 3);
}
   
    float arr[12] = {
        1,2,12,2,
        2,5,2,4,
        7,2,2,4
    
    };
    float arr2[12] = {
        0,0,0,0,
        0,0,0,0,
        0,0,0,0    
    };
    float arr3[4] = {
        2, 3, 
        4, 5   
    };
    float eye[4] = {
        1, 0, 
        0, 1
    };
    Matrix I(eye, 2, 2);

    Matrix A(arr3, 2, 2);
    (A|Matrix::eye(2)).print();
    Matrix::eye(3).print();
    A.inv().print();

    return 0;
}



/*
//some more stuff
// {
 /*interesting experiment here:*/
    
    // for(int i{0}; i<3; i++){
    //     for(int j{0}; j<3;j++){
    //         // std::cout<<arr+i<<" "<<*(arr+i)<<" "<<**(arr+i)<<std::endl;
    //         *(*(arr + i)+i) = i;

    //         std::cout<<i<<","<<j<<" "<<&arr[i][j]<<" "<<*&arr[i][j]<<"\n";
            
        

    //     }


    // }
    // std::cout<<std::endl;
    // for(int i{0}; i<30; i++){
    //     std::cout<<i<<" "<<(*arr + i)<<" "<<*(*arr + i)<<std::endl;
    // }
    // std::cout<<std::endl;
