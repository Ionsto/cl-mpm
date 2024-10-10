#include <iostream>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
#include "src/constitutive.h"

int main(int args, char **argv){
  typedef std::chrono::high_resolution_clock Clock;
  std::cout<<"Starting test\n";
  Eigen::Matrix<double,6,1> strain;// {0.0,0.0,0.0,1.0,0.0,0.0};
  strain = (Eigen::Matrix<double,6,1>()<<
            10.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0).finished();
  double E = 10;
  double nu = 0.25;
  double phi = 0.1;
  double psi = 0.00;
  double c = 0;
  std::cout << "Drucker-prager" << "\n";
  auto t1 = Clock::now();
  int iters = 100000;
  // for (int i = 0; i< iters;++i)
  // {
  //     Eigen::Matrix<double,6,1> strainE = swizzle_voigt_coombs(DruckerPrager(swizzle_coombs_voigt(strain),E,nu,phi,psi,c));
  // }
  Eigen::Matrix<double,6,1> strainE = DruckerPrager(strain,E,nu,phi,psi,c);
  auto t2 = Clock::now();
  std::cout << "Took: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " seconds \n";
  std::cout << "Throughput: " << iters / (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " su/seconds \n";
  std::cout<<"finished test\n";
  std::cout << "Elastic strain" << "\n";
  // for(int i = 0;i < 6;++i) {
  //   std::cout << strain[i] << ", ";
  // }
  std::cout << "\nNew strain" << "\n";
  for(int i = 0;i < 6;++i) {
    std::cout << strainE[i] << ", ";
  }
  
  // std::cout << "Test Q\n";
  // Eigen::Matrix<double,3,3> strain = (Eigen::Matrix<double,3,3>() <<
  //                                     1.0,2.0,3.0,
  //                                     4.0,5.0,6.0,
  //                                     7.0,8.0,9.0).finished();
  // std::cout << strain <<"\n";
  // auto Q = AssembleQMatrix(strain);
  // std::cout << Q <<"\n";
}
