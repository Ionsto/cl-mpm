#pragma once
#include "kernel.h"

namespace gpu{
  void setup_sim(Sim &sim){
    Mesh & m = *sim.mesh;
    int data_size = (m.node_count * 3);
    double * mass = m.mass.data();
    double * disp = m.displacement.data();
    double * vel = m.velocity.data();
    double * force = m.force.data();
#pragma omp target enter data map(to:mass[0:data_size], disp[0:data_size], vel[0:data_size], force[0:data_size], mass[0:data_size])
  }
  void sync_sim(Sim &sim){
    Mesh & m = *sim.mesh;
    int data_size = (m.node_count * 3);
    double * mass = m.mass.data();
    double * disp = m.displacement.data();
    double * vel = m.velocity.data();
    double * force = m.force.data();
#pragma omp target update from (mass[0:data_size], disp[0:data_size], vel[0:data_size], force[0:data_size], mass[0:data_size])
  }
  void close_sim(Sim &sim){
    Mesh & m = *sim.mesh;
    int data_size = (m.node_count * 3);
    double * mass = m.mass.data();
    double * disp = m.displacement.data();
    double * vel = m.velocity.data();
    double * force = m.force.data();
#pragma omp target exit data map(from:mass[0:data_size], disp[0:data_size], vel[0:data_size], force[0:data_size], mass[0:data_size])
  }
  void update_mps(Sim &sim){
    update_mps_kernel(sim.mps.data(),sim.mesh->displacement.data(),sim.dt,sim.mps.size());
  }
  void update_mps_kernel(mp* mps,double *disp,double dt, int mp_count){
#pragma omp target teams parallel for
    for(int i = 0; i < mp_count;++i){
      mps[i].strain[0] = 0;
    }
  }

  void integrate(Sim &sim){
    auto &m = *sim.mesh;
    integration_kernel(
                       m.mass.data(),
                       m.displacement.data(),
                       m.velocity.data(),
                       m.force.data(),
                       sim.damping_factor,
                       sim.dt,
                      (m.node_count * 3));
  }

  void integration_kernel(double * mass, double * disp, double * vel,double * force,double damping_factor, double dt, int data_size){
  #pragma omp target enter data map(to:mass[0:data_size], disp[0:data_size], vel[0:data_size], force[0:data_size], mass[0:data_size])
  #pragma omp target teams parallel for
    for(int i = 0; i < data_size;++i){
      force[i] += mass[i]*vel[i] * -1.0 * damping_factor;
      vel[i] += force[i]/mass[i] * dt;
      disp[i] += disp[i]/mass[i] * dt;
    }
  #pragma omp target exit data map(from:disp[0:data_size],vel[0:data_size])
  }

  void reset_force(Sim &sim){
    auto &m = *sim.mesh;
    reset_force_kernel(m.force.data(), (m.node_count * 3));
  }
  void reset_force_kernel(double * force, int data_size){
#pragma omp target teams parallel for
    for(int i = 0; i < data_size;++i){
      force[i] = 0.0;
    }
  }


  void apply_bcs(Sim &sim){
    auto &m = *sim.mesh;
    apply_bcs_kernel(m.bcs.data(),
              m.displacement.data(),
              m.velocity.data(),
              m.force.data(),
              (m.node_count * 3));
  }

  void apply_bcs_kernel(double * bcs,double * disp, double *vel,double * force,int data_size){
#pragma omp target enter data map(to:disp[0:data_size], vel[0:data_size], force[0:data_size])
#pragma omp target teams parallel for
    for(int i = 0; i < data_size;++i){
      vel[i] *= bcs[i];
      disp[i] *= bcs[i];
      force[i] *= bcs[i];
    }
#pragma omp target exit data map(from:disp[0:data_size],vel[0:data_size],force[0:data_size])
  }
}

/* Eigen::Matrix<double,6,1> log_strain_update_gpu(const Eigen::Matrix<double,6,1> & strain, const Eigen::Matrix3d& df){ */
/*   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(voigt_to_matrix(strain)); */
/*   if (eigensolver.info() != Eigen::Success) */
/*     { */
/*       // std::cout<<"Eigensolve failed\n"; */
/*       // abort(); */
/*       return Eigen::Matrix<double,6,1>::Zero(); */
/*     } */
/*   auto eigen_values = eigensolver.eigenvalues(); */
/*   auto eigen_vectors = eigensolver.eigenvectors(); */
/*   auto trial_lgs = df * (eigen_vectors */
/*                          * (eigen_values.array() * 2.0).exp().matrix().asDiagonal() */
/*                          * eigen_vectors.transpose()) * df.transpose(); */
/*   //0.5 * (trial_lgs + trial_lgs.transpose()) */
/*   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> trialeigensolver(trial_lgs); */
/*   if (trialeigensolver.info() != Eigen::Success) */
/*     { */
/*       // std::cout<<"Eigensolve failed\n"; */
/*       // abort(); */
/*       return Eigen::Matrix<double,6,1>::Zero(); */
/*     } */
/*   auto l = trialeigensolver.eigenvalues(); */
/*   auto v = trialeigensolver.eigenvectors(); */
/*   if ((l.array() <= 0.0).any()) */
/*     { */
/*       // std::cout<<"Eigensolve failed\n"; */
/*       // abort(); */
/*       // return false; */
/*       return Eigen::Matrix<double,6,1>::Zero(); */
/*     } */
/*   return (matrix_to_voigt(v * l.array().log().matrix().asDiagonal() * v.transpose()).array() * 0.5).matrix(); */
/* } */

/* void gpu_log_strain_update(mp * mps, int size){ */
/* #pragma omp target teams parallel for */
/*   for(int i = 0; i < size;++i){ */
/*     mps[i].strain = log_strain_update_gpu(mps[i].strain_n,mps[i].df); */
/*   } */
/* } */
