#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using VoigtMatrix = Eigen::Matrix<double,6,1>;

Eigen::Matrix<double,6,1> swizzle_voigt_coombs(Eigen::Matrix<double,6,1> voigt)
{
  return (Eigen::Matrix<double,6,1>() <<
          voigt[0],
          voigt[1],
          voigt[2],
          voigt[5],
          voigt[3],
          voigt[4]
          ).finished();
}

Eigen::Matrix<double,6,1> swizzle_coombs_voigt(Eigen::Matrix<double,6,1> voigt)
{
  return (Eigen::Matrix<double,6,1>() <<
          voigt[0],
          voigt[1],
          voigt[2],
          voigt[4],
          voigt[5],
          voigt[3]
          ).finished();
}

Eigen::Matrix<double,6,1> strain_to_stress(Eigen::Matrix<double,6,1> strain) {
  return (Eigen::Matrix<double,6,1>() <<
          1.0,
          1.0,
          1.0,
          0.5,
          0.5,
          0.5
          ).finished().cwiseProduct(strain);
}

Eigen::Matrix<double,6,1> stress_to_strain(Eigen::Matrix<double,6,1> stress) {
  return (Eigen::Matrix<double,6,1>() <<
          1.0,
          1.0,
          1.0,
          2.0,
          2.0,
          2.0
          ).finished().cwiseProduct(stress);
}

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

Eigen::Matrix<double,6,6> elastic_matrix(double E, double nu) {
  const double g = 0.5 * (1.0 - (2 * nu));
  const double k = 1.0 - nu;
  return (Eigen::Matrix<double,6,6>() <<
          k, nu, nu, 0, 0, 0,
          nu, k, nu, 0, 0, 0,
          nu, nu, k, 0, 0, 0,
          0,0,0,g,0,0,
          0,0,0,0,g,0,
          0,0,0,0,0,g
          ).finished() * (E / ((1+nu)*(1 - (2 * nu))));
}
