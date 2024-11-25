#include "matrix.hpp"
#include <iostream>

// float **make_arr(){
//     int m;
//     int n;
//     std::cout<<"Please enter the number of rows: ";
//     std::cin>>m;
//     std::cout<<"Please enter the number of columns: ";
//     std::cin>>n;
//     std::cout<<"Creating an "<<m<<"x"<<n<<" matrix\nEnter the following float values:"<<std::endl;

//     float **arr = new float*[m];
//     for(int i{0}; i<m; i++){
//         arr[i] = new float[n];
//         for(int j{0}; j<n;j++){
//             std::cout<<i<<"-"<<j<<"th entry: ";
//             std::cin>>arr[i][j];
//         }
//     }
//     return arr; //can i do this??
// }
// void delete_arr(float **&arr, int m, int n){
//     for(int i{0}; i<m; i++){
//         for(int j{0}; j<n; j++){
//             arr[i][j] = 0;
//         }
//         delete[] arr[i];
//         arr[i] = nullptr;
//     }
//     delete arr;
//     arr = nullptr;
    
// }


int main(){

    Matrix v1(new float[3]{
        1, 
        2,
        -1
    }, 3, 1);

    Matrix v2(new float[3]{
        3, 
        6,
        -1
    }, 3, 1);

    Matrix v3(new float[3]{
        10, 
        19,
        -6
    }, 3, 1);

    Matrix set[3]{
        v1, v2, v3
    };

    Matrix s1(set, 3);

    s1.rref().print();



   


    return 0;
}


