#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include "utils.h"
#include <iostream>

Eigen::Matrix<double,6,6> AssembleQMatrix(Eigen::Matrix<double,3,3> eigen_vectors){
  Eigen::Matrix<double,6,6> Q;
  for(int j = 0;j < 3;++j){
    for(int i = 0;i < 3;++i){
      Q(i,j)   = eigen_vectors(i,j) * eigen_vectors(i,j);
      Q(i+3,j) = eigen_vectors(i,j) * eigen_vectors((i+1)%3,j);
    }
  }
  // Q.block(0,0,3,3) = eigen_vectors.array().square().matrix();
  for(int j = 0;j < 3;++j){
    for(int i = 0;i < 3;++i){
      Q(i,j+3)   = 2*eigen_vectors(i,j) * eigen_vectors(i,(j+1)%3);

      Q(i+3,j+3) =
        (eigen_vectors(i,j) * eigen_vectors((i+1)%3,(j+1)%3))
        +
        eigen_vectors((i+1)%3,j) * eigen_vectors(i,(j+1)%3);
    }
  }
  return Q.transpose();
}

/*
  Because the performance of MAGICL is terrible for allocations it is actually highly
  benificial to run complicated constitutive models in c++ leaning on the optimised
  eigen library to do all linalg operations
*/
Eigen::Matrix<double,6,1> DruckerPrager(Eigen::Matrix<double,6,1> elastic_strain,
                    double E, double nu,
                    double phi, double psi, double c) {
  const double alfa = -std::tan(phi);
  const double bta = -std::tan(psi);
  const double xsic = std::sqrt(3)*(1.0 / std::tan(phi))*c;
  Eigen::Matrix<double,3,3> Ce = (Eigen::Matrix<double,3,3>()<<
                                  E-nu,-nu,-nu,
                                  -nu,E-nu,-nu,
                                  -nu,-nu,E-nu).finished();
  Eigen::Matrix<double,3,3> De3 =
    (E/((1+nu)/(1-(2*nu))))*
    (((1-(2*nu))*Eigen::Matrix<double,3,3>::Identity()) +
     Eigen::Matrix<double,3,3>::Constant(nu));
  // (Eigen::Matrix<double,3,3>()<<
  //                                E-nu,-nu,-nu,
  //                                -nu,E-nu,-nu,
  //                                -nu,-nu,E-nu).finished();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(voigt_to_matrix(elastic_strain));
  // elastic_strain = swizzle_voigt_coombs(elastic_strain);
  if (eigensolver.info() != Eigen::Success)
    {
      abort();
      /* return false; */
    }
  Eigen::Matrix<double,3,1> eigen_values = eigensolver.eigenvalues().reverse();
  std::cout <<"eigenvalues\n"<<eigen_values<<"\n";
  Eigen::Matrix<double,3,3> eigen_vectors = eigensolver.eigenvectors().rowwise().reverse();
  Eigen::Matrix<double,3,1> sig = ((De3 * eigen_values).array() - (xsic / std::sqrt(3))).matrix();
  const double tol = 1e-12;
  double xi = sig.sum()/std::sqrt(3);
  Eigen::Matrix<double,3,1> s = (sig.array() - (xi / std::sqrt(3))).matrix();
  double rho = std::sqrt(s.array().square().sum());
  double f = rho - alfa*xi;
  if (f>tol){
    Eigen::Matrix<double,3,1> epsE = Ce * sig;
    Eigen::Matrix<double,3,1> epsEtr = epsE;
    auto Q = AssembleQMatrix(eigen_vectors);
    double fap = rho * std::sqrt(1+nu) + (xi*std::sqrt(1-(2*nu))/(bta*std::sqrt(1+nu)/std::sqrt(1-2*nu)));
    if(fap<tol) {
      //Apex return
      sig=Eigen::Matrix<double,3,1>::Constant(xsic/sqrt(3));
    }
    else{
      //Apex return
      //Setup NR algorithm
      Eigen::Matrix<double,4,1> b;
      b<<0.0,0.0,0.0,f;
      int itnum = 0;
      double dgam = 0;
      Eigen::Matrix<double,3,1> df = ((s.array()/rho) - alfa/std::sqrt(3)).matrix();
      Eigen::Matrix<double,3,1> dg = ((s.array()/rho) - bta/std::sqrt(3)).matrix();
      Eigen::Matrix<double,3,3> ddg =
        (((1/(3*rho))
          *
          ((3*Eigen::Matrix<double,3,3>::Identity()) -
           Eigen::Matrix<double,3,3>::Constant(1.0)))
         - s*s.transpose()/std::pow(rho,3));
      std::cout<<"df\n"<<df<<"\n";
      std::cout<<"dg\n"<<dg<<"\n";
      std::cout<<"ddg\n"<<ddg<<"\n";
      const int maxit = 4;
      const double tolf = 1e-6;
      for(int iter = 0; (iter < maxit) && ((b.block(0,0,2,1).norm() > tol) || (std::abs(b(3)) > tolf));++iter){
        Eigen::Matrix<double,4,4> A;
        // A.block(0,0,3,3) << Eigen::Matrix<double,3,3>::Identity();
        // // A << Eigen::Matrix<double,3,3>::Identity() << dg;;// << df.transpose() * De3 << 0.0; 
        A.block(0,0,3,3) << Eigen::Matrix<double,3,3>::Identity() + (dgam * ddg * De3);
        A.block(0,3,3,1) << dg;
        A.block(3,0,1,3) << df.transpose() * De3;
        A(3,3) = 0.0;
        //A.transposeInPlace();

        auto dx = A.partialPivLu().solve(b) * -1.0;
        epsE += dx.block(0,0,3,1);
        dgam += dx(3);
        sig = De3*epsE;
        xi = sig.sum()/std::sqrt(3);
        s = (sig.array() - (xi / std::sqrt(3))).matrix();
        rho = std::sqrt(s.array().square().sum());
        df = ((s.array()/rho) - alfa/std::sqrt(3)).matrix();
        dg = ((s.array()/rho) - bta/std::sqrt(3)).matrix();
        ddg =
          (((1/(3*rho))
            *
            ((3*Eigen::Matrix<double,3,3>::Identity()) -
             Eigen::Matrix<double,3,3>::Constant(1.0)))
           - s*s.transpose()/std::pow(rho,3));

        b.block(0,0,3,1) =(epsE-epsEtr) + dg * dgam;
        b(3) = rho-alfa*xi;
      }
      sig.array() += xsic/std::sqrt(3);
    }
    epsE = Ce*sig;
    // Eigen::Matrix<double,6,1> eps_q;
    return swizzle_coombs_voigt(Q.partialPivLu().solve((Eigen::Matrix<double,6,1>()
                                         <<
                                         epsE[0],
                                         epsE[1],
                                         epsE[2],
                                         0.0,0.0,0.0).finished()));
  }
  return swizzle_coombs_voigt(elastic_strain);
}
