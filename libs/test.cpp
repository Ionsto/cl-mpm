#include <iostream>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

Eigen::Matrix<double,6,1> matrix_to_voigt(Eigen::Matrix<double,3,3> mat) {
  return (Eigen::Matrix<double,6,1>() <<
          mat(0,0), mat(1,1),mat(2,2),
          2.0*mat(1,2), 2.0*mat(0,2),2.0*mat(0,1)
          ).finished();
}

Eigen::Matrix<double,3,3> voigt_to_matrix(Eigen::Matrix<double,6,1> voigt) {
  return (Eigen::Matrix3d() <<
          voigt(0), 0.5*voigt(5), 0.5*voigt(4),
          0.5*voigt(5), voigt(1), 0.5*voigt(3),
          0.5*voigt(4), 0.5*voigt(3), voigt(2)).finished();
}
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
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> trialeigensolver(trial_lgs);
    if (trialeigensolver.info() != Eigen::Success)
      {
        return false;
      }
    auto l = trialeigensolver.eigenvalues();
    auto v = trialeigensolver.eigenvectors();
    strain = (matrix_to_voigt(v * l.array().log().matrix().asDiagonal() * v.transpose()).array() * 0.5).matrix();
    //No clue if this works
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
