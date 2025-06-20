#pragma once
#include <Eigen/Dense>
#include "../src/utils.h"
/* #include "../src/constitutive.h" */
#include "types.h"
namespace gpu{
  void setup_sim(Sim &sim);
  void sync_sim(Sim &sim);
  void close_sim(Sim &sim);
  void update_mps(Sim &sim);
  void update_mps_kernel(mp* mps,double *disp,double dt,int mp_count);
  void integrate(Sim &sim);
  void integration_kernel(double * mass, double * disp, double *vel,double * force,double damping_factor, double dt, int data_size);
  void reset_force(Sim&sim);
  void reset_force_kernel(double * force, int data_size);
  void p2g_gather_force(mp * mps,double * force,super_node_cache * snc, nodecache *nc, int data_size);
  void apply_bcs(Sim &sim);
  void apply_bcs_kernel(double *bcs,double * disp, double *vel,double * force,int data_size);

}
