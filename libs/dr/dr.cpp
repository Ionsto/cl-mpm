#pragma once
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "../src/utils.h"
#include "../src/constitutive.h"
#include <chrono>
#include "types.h"
#include "kernel.h"
#include <functional>
#include <mutex>
#include <memory>


Dsvp assemble_dsvp(Vector grads){
  Dsvp dsvp_alloc = Dsvp::Zero();
  dsvp_alloc(0,0) = grads[0];
  dsvp_alloc(1,1) = grads[1];
  dsvp_alloc(2,2) = grads[2];

  dsvp_alloc(3,1) = grads[2];
  dsvp_alloc(3,2) = grads[1];

  dsvp_alloc(4,0) = grads[2];
  dsvp_alloc(4,2) = grads[0];

  dsvp_alloc(5,0) = grads[1];
  dsvp_alloc(5,1) = grads[0];
  return dsvp_alloc;
}

template <typename F>
void iterate_over_neighbours(Sim & sim, mp & mp, F&& func){
  for (auto &n : mp.nc){
    func(n);
  }
}

void combi_assemble_dsvp(Eigen::Ref<Eigen::Matrix3d> df, const Vector & disp, const Vector & grads){
  df(0,0) += grads[0]*disp[0];
  df(1,1) += grads[1]*disp[1];
  df(2,2) += grads[2]*disp[2];

  df(0,1) += grads[1]*disp[0];
  df(0,2) += grads[0]*disp[1];
  df(1,0) += grads[2]*disp[0];

  df(2,0) += grads[0]*disp[2];
  df(1,2) += grads[2]*disp[1];
  df(2,1) += grads[1]*disp[2];
}

// void combi_det_stress(Vector force, Voigt stress, Vector grads){
//   force += grads * stress;
// }

void calculate_strain(Sim & sim, mp & mp){
  Mesh & mesh = *sim.mesh;
  mp.df = Eigen::Matrix3d::Identity();
  iterate_over_neighbours(sim,mp,
                          [&](auto & n){
                            combi_assemble_dsvp(mp.df,mesh.displacement.row(n.node).transpose(),n.grads);
                          });
  mp.f = mp.df * mp.f_n;
  mp.volume = mp.volume_n * mp.df.determinant();
  mp.strain = constitutive::log_strain_update(mp.strain_n,mp.df);
}

void calculate_force(Mesh & mesh, mp & mp){
  for(auto &n : mp.nc){
    std::lock_guard<std::mutex> guard((*mesh.node_locks)[n.node]);
    Dsvp dsvp = assemble_dsvp(mp.df.inverse()*n.grads);
    mesh.force.row(n.node) += (mp.mass * n.svp * mesh.gravity.transpose()) + (-1.0 * mp.volume*(dsvp.transpose() * mp.stress).transpose());
  }
}



int position_to_index(Mesh& m, Eigen::Array3i index){
  return index[0] + (index[1] * m.element_count[0]) + (index[2] * m.element_count[0] * m.element_count[1]);
}

int in_bounds(Mesh& m, Eigen::Array3i index){
  return index[0]>=0 && index[1]>=0  && index[2]>=0 &&
    index[0]<m.element_count[0]&&
    index[1]<m.element_count[1]&&
    (index[2]<m.element_count[2]);
}

double shape_linear(double distance, double h){
  return 1.0 - std::abs(distance/h);
}

double shape_linear_dsvp(double distance, double h){
  // if(distance > 0){
  //   return -1.0;
  // } else if(distance < 0){
  //   return 1.0;
  // }
  // return 0.0;
  return std::copysign(1.0,distance)/h;
}

void make_nc(Mesh & m, mp & mp){
  Eigen::Array3i bottom_index = (mp.position.array()/m.h).floor().cast<int>();
  const double h = m.h;
  for(int dx = 0;dx < 2;++dx){
    for(int dy = 0;dy < 2;++dy){
      for(int dz = 0;dz < 2;++dz){
        Eigen::Array3i dindex;
        dindex<<dx,dy,dz;
        Eigen::Array3i node_index = bottom_index+dindex;
        if(in_bounds(m,node_index)){
          Eigen::Array3d node_position = node_index.cast<double>()*m.h;
          Eigen::Array3d distance = mp.position.array()-node_position;
          Eigen::Array3d weights = distance.matrix().unaryExpr([&](double x){return shape_linear(x,h);});
          Eigen::Array3d grads = distance.matrix().unaryExpr([&](double x){return shape_linear_dsvp(x,h);});
          nodecache nc;
          nc.mp = mp.index;
          nc.svp = weights[0]*weights[1]   ;//*weights[2];
          nc.grads[0] = grads[0]*weights[1];//*weights[2];
          nc.grads[1] = grads[1]*weights[0];//*weights[2];
          nc.grads[2] = 0.0;//grads[2]*weights[0]*weights[1];
          nc.node = position_to_index(m,node_index);
          // nc.dsvp = assemble_dsvp(nc.grads);
          // std::cout <<"mp pos\n";
          // std::cout << mp.position <<"\n";
          // std::cout <<"node pos\n";
          // std::cout << node_position <<"\n";
          // std::cout <<"dsvp\n";
          // std::cout << nc.dsvp<<"\n";
          mp.nc.push_back(nc);

          {
            std::lock_guard<std::mutex> guard((*m.node_locks)[nc.node]);
            m.nodes[nc.node].nc.push_back(nc);
          }
        }
      }
    }
  }
}


void make_gpu_nc(Sim & sim, mp & mp){
  Mesh &m = *sim.mesh;
  Eigen::Array3i bottom_index = (mp.position.array()/m.h).floor().cast<int>();
  const double h = m.h;
  for(int dx = 0;dx < 2;++dx){
    for(int dy = 0;dy < 2;++dy){
      for(int dz = 0;dz < 2;++dz){
        Eigen::Array3i dindex;
        dindex<<dx,dy,dz;
        Eigen::Array3i node_index = bottom_index+dindex;
        if(in_bounds(m,node_index)){
          Eigen::Array3d node_position = node_index.cast<double>()*m.h;
          Eigen::Array3d distance = mp.position.array()-node_position;
          Eigen::Array3d weights = distance.matrix().unaryExpr([&](double x){return shape_linear(x,h);});
          Eigen::Array3d grads = distance.matrix().unaryExpr([&](double x){return shape_linear_dsvp(x,h);});
          nodecache nc;
          nc.mp = mp.index;
          nc.svp = weights[0]*weights[1]   ;//*weights[2];
          nc.grads[0] = grads[0]*weights[1];//*weights[2];
          nc.grads[1] = grads[1]*weights[0];//*weights[2];
          nc.grads[2] = 0.0;//grads[2]*weights[0]*weights[1];
          nc.node = position_to_index(m,node_index);
          // nc.dsvp = assemble_dsvp(nc.grads);

          sim.global_nodecache.push_back(nc);

          // mp.nc.push_back(nc);
          // {
          //   std::lock_guard<std::mutex> guard((*m.node_locks)[nc.node]);
          //   m.nodes[nc.node].nc.push_back(nc);
          // }
        }
      }
    }
  }
}




void setup_mps(Sim & sim){
#pragma omp parallel for
  for(int i = 0;i < sim.mps.size();++i){
    auto & mp = sim.mps[i];
    mp.index = i;
    // make_nc(*sim.mesh,mp);
    mp.set_elastic(1e4,0.3);
    double density = 1;
    mp.mass = density * mp.volume;
  }
  for(int i = 0;i < sim.mps.size();++i){
    make_gpu_nc(sim,sim.mps[i]);
  }

}
void p2g(Sim & sim){
#pragma omp parallel for
    for(int i = 0;i < sim.mps.size();++i){
      auto & mp = sim.mps[i];
      iterate_over_neighbours(sim,mp,
                              [&](auto & nc){
                                std::lock_guard<std::mutex> guard((*sim.mesh->node_locks)[nc.node]);
                                sim.mesh->mass.row(nc.node) += Vector::Constant(mp.mass * nc.svp);
                              });
    }
}

void g2p(Sim & sim){
#pragma omp parallel for
  for(int i = 0;i < sim.mps.size();++i){
    auto & mp = sim.mps[i];
    mp.disp_inc *= 0.0;
    iterate_over_neighbours(sim,mp,
                            [&](auto & nc){
                              //mp.disp_inc += nc.svp * sim.mesh->displacement.row(node.index).transpose();
                              mp.disp_inc += nc.svp * sim.mesh->displacement.row(nc.node).transpose();
                            });
  }
}

void update_mps(Sim & sim, double dt){
#pragma omp parallel for
  for(int i = 0;i < sim.mps.size();++i){
    auto & mp = sim.mps[i];
    calculate_strain(sim,mp);
    mp.stress = mp.de * mp.strain;
    mp.stress = mp.stress/mp.f.determinant();

  }
}


void p2g_node_force(Sim & sim, node & node){
  for(auto &n : node.nc){
    auto & mp = sim.mps[n.mp];
    Dsvp dsvp = assemble_dsvp(mp.df.inverse()*n.grads);
    sim.mesh->force.row(n.node) += (-1.0*mp.volume*(dsvp.transpose() * mp.stress).transpose()) + (mp.mass * n.svp * sim.mesh->gravity.transpose());
  }
}

void p2g_scatter_force(Sim &sim){
#pragma omp parallel for
  for(int i = 0;i < sim.mps.size();++i){
    auto & mp = sim.mps[i];
    calculate_force(*sim.mesh,mp);
  }
}

void p2g_gather_force(Sim &sim){
#pragma omp parallel for
  for(int i = 0;i < sim.mesh->nodes.size();++i){
    auto & node = sim.mesh->nodes[i];
    p2g_node_force(sim,node);
  }
}

void integrate(Sim&sim){
  auto &m = *sim.mesh;
  m.force += (m.mass.array()*sim.mesh->velocity.array() * (-1.0 * sim.damping_factor)).matrix();
  m.velocity += sim.dt*(m.force.array()/m.mass.array()).matrix();
  m.displacement += sim.dt*m.velocity;
}
void reset_force(Sim& sim){
  sim.mesh->force *= 0.0;
}
void compute_residuals(){
}

void apply_bcs(Sim & sim){
  Mesh &mesh = *sim.mesh;
  mesh.displacement = (mesh.displacement.array() * mesh.bcs.array()).matrix();
  mesh.force = (mesh.force.array() * mesh.bcs.array()).matrix();
  mesh.velocity = (mesh.velocity.array() * mesh.bcs.array()).matrix();
}

void save_vtk(std::string filename,std::vector<mp> &mps){
  std::ofstream myfile;
  myfile.open(filename, std::ios::out);
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "C++ generated vtk file, SJVS\n";
  myfile << "ASCII\n";
  myfile << "DATASET UNSTRUCTURED_GRID\n";
  myfile << "POINTS " << mps.size() << " double\n";
  for(auto & mp : mps){
    Eigen::Vector3d disp_pos = mp.position + mp.disp_inc;
    myfile << disp_pos[0] << " " << disp_pos[1] << " " << disp_pos[2] << "\n";
  }
  myfile << "POINT_DATA " << mps.size() << "\n";
  myfile << "SCALARS s_yy FLOAT 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(auto & mp : mps){
    myfile << mp.stress[1] << "\n";
  }
  myfile << "SCALARS e_yy FLOAT 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(auto & mp : mps){
    myfile << mp.strain[1] << "\n";
  }
  myfile.close();
}

// void save_param(std::ofstream & stream, std::string name, std::function lambda){

// }

void save_vtk_nodes(std::string filename,Mesh &mesh){
  std::vector<node> & nodes = mesh.nodes;
  std::ofstream myfile;
  myfile.open(filename, std::ios::out);
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "C++ generated vtk file, SJVS\n";
  myfile << "ASCII\n";
  myfile << "DATASET UNSTRUCTURED_GRID\n";
  myfile << "POINTS " << nodes.size() << " double\n";
  for(auto & node : nodes){
    Eigen::Vector3d disp_pos = node.position + mesh.displacement.row(node.index).transpose();
    myfile << disp_pos[0] << " " << disp_pos[1] << " " << disp_pos[2] << "\n";
  }
  myfile << "\nPOINT_DATA " << nodes.size() << "\n";
  myfile << "SCALARS mass FLOAT 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(auto & n : mesh.nodes){
    myfile << mesh.mass(n.index,0) << "\n";
  }
  myfile << "VECTORS force FLOAT\n";
  // myfile << "LOOKUP_TABLE default\n";
  for(auto & node : nodes){
    Eigen::Vector3d f = mesh.force.row(node.index).transpose();
    myfile << f[0] << " " << f[1] << " " << f[2] << "\n";
  }
  myfile << "VECTORS velocity FLOAT\n";
  for(auto & node : nodes){
    Eigen::Vector3d f = mesh.velocity.row(node.index).transpose();
    myfile << f[0] << " " << f[1] << " " << f[2] << "\n";
  }


  myfile.close();
}
void save(Sim &sim, int index){
  std::ostringstream ss;
  ss << "./output/test_"<<index<<".vtk";
  save_vtk(ss.str(),sim.mps);
  std::ostringstream ssn;
  ssn << "./output/test_n_"<<index<<".vtk";
  save_vtk_nodes(ssn.str(),*sim.mesh);
}


int main(int argc,char ** args) {
  typedef std::chrono::high_resolution_clock Clock;
  const double h = 1;
  Eigen::Vector3d block_size = (Eigen::Vector3d()<<100.0,100.0,1.0).finished();
  Sim sim(h,std::round(block_size[0]/h)+1,std::round(block_size[1]/h)+1,+1);
  double mps_per_cell = 2;
  double mp_res = (h/(mps_per_cell));
  double offset = mp_res*0.5;//(h/(1+mps_per_cell));
  int x_count = std::floor(block_size[0]/mp_res);
  int y_count = std::floor(block_size[1]/mp_res);
  // int x_count = 2;
  // int y_count = 2;
  int mp_count = x_count * y_count;
  sim.mps = std::vector<mp>(mp_count);
  int i = 0;
  for(int x = 0; x < x_count;++x){
    for(int y = 0; y < y_count;++y){
      sim.mps[i].position[0] = (x * mp_res) + offset;
      sim.mps[i].position[1] = (y * mp_res) + offset;
      sim.mps[i].volume = mp_res * mp_res;
      sim.mps[i].volume_n = sim.mps[i].volume;
      i++;
    }
  }
  std::cout<<"Setup\n";
  setup_mps(sim);
  save(sim,0);
  // save_vtk("./test_0.vtk",mps);
  // save_vtk_nodes("./test_n_0.vtk",m);

  std::cout << "Node count: " << sim.mesh->node_count << "\n";
  std::cout << "NDOFs: " << (3 * sim.mesh->node_count) << "\n";
  std::cout << "MPs: " << mp_count << "\n";


  int iters = 100;
  int substeps = 100;
  std::cout << "Iters: " << iters << "\n";
  std::cout << "Substeps: " << substeps << "\n";
  sim.dt = 0.25 * sim.estimate_elastic_dt();
  sim.damping_factor = sim.estimate_critical_damping();
  std::cout<<"Estimated timestep: " << sim.dt<<"\n";
  std::cout<<"Estimated damping: " << sim.damping_factor<<"\n";
  p2g(sim);
  auto t1 = Clock::now();
  gpu::setup_sim(sim);
  for(int i = 0;i < iters;++i){
    std::cout<<"Step: " << i <<"\n";
    for(int substep = 0;substep < substeps;++substep){
      gpu::reset_force(sim);
      // update_mps(sim,sim.dt);
      // p2g_gather_force(sim);
      // gpu::integrate(sim);
      // gpu::apply_bcs(sim);
    }
    // g2p(sim);

    gpu::sync_sim(sim);
    save(sim,i+1);
  }
  gpu::close_sim(sim);
  auto t2 = Clock::now();
  std::cout << "Took: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " seconds \n";
  // std::cout << "Throughput: " << iters / (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " su/seconds \n";
  std::cout << "MP time: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count()/(iters *substeps* mp_count)  << " su/seconds \n";
  return 0;
}
