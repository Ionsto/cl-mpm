#include <iostream>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "src/utils.h"
#include "src/constitutive.h"

extern "C" {

  void test(double * strain_ptr){
    strain_ptr[0] = 1.0;
  }

  bool Kirchoff_Strain_Update(double * strain_ptr, double * df_ptr){
    Eigen::Map<Eigen::Matrix<double,6,1>> strain(strain_ptr);
    Eigen::Map<Eigen::Matrix<double,3,3>> df(df_ptr);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(voigt_to_matrix(strain));
    if (eigensolver.info() != Eigen::Success)
      {
        return false;
      }
    auto eigen_values = eigensolver.eigenvalues();
    auto eigen_vectors = eigensolver.eigenvectors();
    auto trial_lgs = df * (eigen_vectors
                           * (eigen_values.array() * 2.0).exp().matrix().asDiagonal()
                           * eigen_vectors.transpose()) * df.transpose();
    //0.5 * (trial_lgs + trial_lgs.transpose())
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> trialeigensolver(trial_lgs);
    if (trialeigensolver.info() != Eigen::Success)
      {
        return false;
      }
    auto l = trialeigensolver.eigenvalues();
    auto v = trialeigensolver.eigenvectors();
    if ((l.array() <= 0.0).any())
      {
        //Bad news we've got negative or zero eigenvalues - this is not good
        return false;
      }
    //Nasty hack we may get very small negative eigenvalues - so we take max of 0
    strain = (matrix_to_voigt(v * l.array().log().matrix().asDiagonal() * v.transpose()).array() * 0.5).matrix();
    //No clue if this works
    return true;
  }

  // bool MC_plastic_model(double * strain_ptr, double )


  bool CppDruckerPrager(double * strain_ptr,double E, double nu, double phi, double psi, double c)
  {
    Eigen::Map<Eigen::Matrix<double,6,1>> strain(strain_ptr);
    Eigen::Matrix<double,6,1> strainE = DruckerPrager(strain,E,nu,phi,psi,c);
    strain = strainE;
    return true;
  }
  bool CppMohrCoulomb(double * strain_ptr,double & f,double & psinc,double E, double nu, double phi, double psi, double c)
  {
    Eigen::Map<Eigen::Matrix<double,6,1>> strain(strain_ptr);
    MohrCoulombReturn result = MohrCoulomb(strain,E,nu,phi,psi,c);
    Eigen::Matrix<double,6,1> strainE = std::get<0>(result);
    //std::cout<<"strain\n"<<strainE<<"\n";
    f = std::get<1>(result);
    psinc = std::get<2>(result);
    strain = strainE;
    return std::get<3>(result);
  }
  bool CppViscoelastic(double * strain_ptr,double E, double nu,double viscosity, double dt) 
  {
    Eigen::Map<Eigen::Matrix<double,6,1>> strain(strain_ptr);
    Eigen::Matrix<double,6,1> strainE = Viscoelastic(strain,E,nu,viscosity,dt);
    strain = strainE;
    return true;
  }
  bool MatrixSqrt(double * input_ptr,double * output_ptr){
    Eigen::Map<Eigen::Matrix<double,3,3>> in(input_ptr);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(in);
    Eigen::Map<Eigen::Matrix<double,3,3>> out(output_ptr);
    out = eigensolver.operatorSqrt();
    return true;
  }
  bool VoigtSqrt(double * input_ptr,double * output_ptr){
    Eigen::Map<Eigen::Matrix<double,6,1>> in(input_ptr);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(voigt_to_matrix_stress(in));
    Eigen::Map<Eigen::Matrix<double,6,1>> out(output_ptr);
    out = matrix_to_voigt_stress(eigensolver.operatorSqrt());
    return true;
  }
}
int main(int arc,char **args){
  typedef std::chrono::high_resolution_clock Clock;
  std::cout<<"Starting test\n";
  auto t1 = Clock::now();
  double strains[6] = {0};
  double df[9] = {1, 0, 0,
                  0, 1, 0,
                  0, 0, 1};

  // double test = 10.0;
  // Eigen::Matrix<double,6,1> strain;
  // strain << 1, 0, 0, 0, 0, 0;
  // Eigen::Matrix<double,3,3> df;
  // df << 1, 0, 0,
  //   0, 1, 0,
  //   0, 0, 1;
  int iters = 10000000;
  for (int i = 0; i< iters;++i)
  {
    Kirchoff_Strain_Update(strains,df);
  //   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(voigt_to_matrix(strain));
  //   if (eigensolver.info() != Eigen::Success)
  //     {
  //       //console_->error("No eigenvectors in strain matrix?\n");
  //       abort();
  //     }
  //   auto eigen_values = eigensolver.eigenvalues();
  //   auto eigen_vectors = eigensolver.eigenvectors();
  //   auto trial_lgs = df * (eigen_vectors
  //                          * (eigen_values.array() * 2.0).exp().matrix().asDiagonal()
  //                          * eigen_vectors.transpose()) * df.transpose();

  //   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> trialeigensolver(trial_lgs);
  //   if (trialeigensolver.info() != Eigen::Success)
  //     {
  //       //console_->error("No eigenvectors in trial strain matrix?\n");
  //       abort();
  //     }
  //   auto l = trialeigensolver.eigenvalues();
  //   auto v = trialeigensolver.eigenvectors();
  //   strain = (matrix_to_voigt(v * l.array().log().matrix().asDiagonal() * v.transpose()).array() * 0.5).matrix();
  //   //strain_(3) *= 2.0;
  //   //strain_(4) *= 2.0;
  //   //strain_(5) *= 2.0;
  }
  auto t2 = Clock::now();
  std::cout << "Took: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " seconds \n";
  std::cout << "Throughput: " << iters / (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " su/seconds \n";
  std::cout<<"finished test\n";
  return 0;
}
