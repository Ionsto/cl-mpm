#include <iostream>
#include <array>
#include <vector>
#include <chrono>
#include <Eigen/Dense>
// #include "types.h"
// #include "hot.cuh"


void integrate_kernel(float *x,float* y,float* z,int size){
#pragma omp target teams parallel for
  for(int i = 0;i < size;++i)
    {
      z[i] = x[i] + y[i];
    }
}

// void add(mp * mps,int size){
// // #pragma omp target map(tofrom:mps[0,size])
// #pragma omp target teams parallel for
//   for(int i = 0;i < size;++i)
//     {
//       mps[i].stress = mps[i].de *  mps[i].strain;
//     }
// }

int main(int argc, char **args){
  typedef std::chrono::high_resolution_clock Clock;
  typedef Eigen::Matrix<float,Eigen::Dynamic,3> DofType;
  const int data_size = 10000;
  const int iters = 100;
  const int total_iters = iters;
  // std::vector<mp> x(data_size);

  DofType x = DofType::Zero(data_size,3);
  DofType y = DofType::Zero(data_size,3);
  DofType z = DofType::Zero(data_size,3);

  std::cout << "Mps: "<< data_size<<"\n";
  std::cout << "Iters: "<< iters <<"\n";
#pragma omp target enter data map(to:x.data()[0:data_size],y.data()[0:data_size])
  auto t1 = Clock::now();
  for(int i = 0;i < iters;++i){
    integrate_kernel(x.data(),y.data(),z.data(),data_size);
    // add(x.data(),data_size);
  }
  auto t2 = Clock::now();
#pragma omp target exit data map(from:z.data()[0:data_size])
  // int print_size = std::min(3,data_size);
  // for(int i = 0;i < print_size;++i){
  //   std::cout<<x[i].strain[0]<<"\n";
  // }
  std::cout << "Took: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " seconds \n";
  std::cout << "Throughput: " << total_iters / (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " GFLOPs \n";

}
