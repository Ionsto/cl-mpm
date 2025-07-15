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
  // for(int i = 0;i < mp.snc.count;++i){
  //   auto &n = sim.global_nodecache[sim.global_nodecache_mp_indirection[mp.snc.index + i]];
  //   func(n);
  // }
  for (auto &n : mp.nc){
    func(n);
  }
}
template <typename F>
void iterate_over_node_neighbours(Sim & sim, node & node, F&& func){
  // for(int i = 0;i < node.snc.count;++i){
  //   auto &n = sim.global_nodecache[sim.global_nodecache_node_indirection[node.snc.index + i]];
  //   func(n);
  // }
  for (auto &n : node.nc){
    func(n);
  }
}

void combi_assemble_dsvp(Eigen::Ref<Eigen::Matrix3d> df, const Vector & disp, const Vector & grads){
  df(0,0) += grads[0]*disp[0];
  df(1,1) += grads[1]*disp[1];
  df(2,2) += grads[2]*disp[2];

  //Correct?
  df(0,1) += grads[0]*disp[1];
  df(0,2) += grads[0]*disp[2];
  df(1,0) += grads[1]*disp[0];
  df(2,0) += grads[2]*disp[0];
  df(1,2) += grads[1]*disp[2];
  df(2,1) += grads[2]*disp[1];

  // df(1,0) += grads[1]*disp[0];
  // df(2,0) += grads[0]*disp[1];
  // df(0,1) += grads[0]*disp[1];
  // df(0,2) += grads[0]*disp[2];
  // df(2,1) += grads[2]*disp[1];
  // df(1,2) += grads[1]*disp[2];
}

// void combi_det_stress(Vector force, Voigt stress, Vector grads){
//   force += grads * stress;
// }

void calculate_strain(Sim & sim, mp & mp){
  Mesh & mesh = *sim.mesh;
  mp.df = Eigen::Matrix3d::Identity();
  iterate_over_neighbours(sim,mp,
                          [&](auto & n){
                            //Stretch increment on df
                            combi_assemble_dsvp(mp.df,mesh.displacement.row(n.node).transpose(),n.grads);
                          });
  mp.f = mp.df * mp.f_n;
  mp.volume = mp.volume_n * mp.df.determinant();
  // mp.strain = constitutive::log_strain_update(mp.strain_n,mp.df);
  mp.strain = mp.strain_n + utils::matrix_to_voigt(0.5*(mp.df + mp.df.transpose() - (2*Eigen::Matrix3d::Identity())));
  // mp.strain(3) *= 0.5;
  // mp.strain(4) *= 0.5;
  // mp.strain(5) *= 0.5;
}

void calculate_force(Mesh & mesh, mp & mp){
  for(auto &n : mp.nc){
    std::lock_guard<std::mutex> guard((*mesh.node_locks)[n.node]);
    Dsvp dsvp = assemble_dsvp(mp.df.inverse()*n.grads);
    // Dsvp dsvp = assemble_dsvp(n.grads);
    Vector internal_force = (-1.0 * mp.volume*(dsvp.transpose() * mp.stress).transpose());
    Vector external_force = (mp.mass * n.svp * mesh.gravity.transpose());
    // mesh.force.row(n.node) += (mp.mass * n.svp * mesh.gravity.transpose())
    //   + (-1.0 * mp.volume*(dsvp.transpose() * mp.stress).transpose());
    mesh.force.row(n.node) += external_force + internal_force;
    mesh.force_int.row(n.node) += internal_force;
    mesh.force_ext.row(n.node) += external_force;
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
  if(distance > 0){
    return -1/h;
  } else if (distance < 0){
    return 1/h;
  } else {
    return 0;
  }
}

template <typename F>
void iterate_over_shape_linear(Mesh & mesh, mp & mp, F && func){
  const double h = mesh.h;
  Eigen::Array3i bottom_index = (mp.position.array()/h).floor().cast<int>();
  for(int dx = 0;dx < 2;++dx){
    for(int dy = 0;dy < 2;++dy){
      for(int dz = 0;dz < 2;++dz){
        Eigen::Array3i dindex;
        dindex<<dx,dy,dz;
        Eigen::Array3i node_index = bottom_index+dindex;
        if(in_bounds(mesh,node_index)){
          Eigen::Array3d node_position = node_index.cast<double>()*h;
          Eigen::Array3d distance = mp.position.array()-node_position;
          Eigen::Array3d weights = distance.matrix().unaryExpr([&](double x){return shape_linear(x,h);});
          Eigen::Array3d grads = distance.matrix().unaryExpr([&](double x){return shape_linear_dsvp(x,h);});
          double svp = weights[0]*weights[1];//*weights[2];
          grads[0] = grads[0]*weights[1];//*weights[2];
          grads[1] = grads[1]*weights[0];//*weights[2];
          grads[2] = 0.0;//grads[2]*weights[0]*weights[1];
          func(mesh.nodes[position_to_index(mesh,node_index)],svp,grads.matrix());
        }
      }
    }
  }
}

double shape_gimp(double distance,double l, double h){
  const double ax = std::abs(distance);
  if(ax <= l){
    return 1 - (((ax*ax)+(l*l))/(2 * h * l));
  } else if (ax <= (h-l)){
    return 1-(ax/h);
  } else if (ax <= h+l){
    return std::pow((h+l)-ax,2)/(4*h*l);
  } else{
    return 0;
  }
}

double shape_gimp_dsvp(double distance,double l, double h){
  const double ax = std::abs(distance);
  if(ax <= l){
    return distance/(h * l);
  } else if (ax <= (h-l)){
    return std::copysign(1.0,distance) * 1/h;
  } else if (ax <= h+l){
    return std::copysign(1.0,distance) * (((-ax)+h+l)/(2*h*l));
  } else{
    return 0;
  }
}

template <typename F>
void iterate_over_shape_gimp(Mesh & mesh, mp & mp, F && func){
  const double h = mesh.h;
  Eigen::Array3i bottom_index = (mp.position.array()/h).round().cast<int>();
  for(int dx = -1;dx <= 1;++dx){
    for(int dy = -1;dy <= 1;++dy){
      for(int dz = -1;dz <= 1;++dz){
        Eigen::Array3i dindex;
        dindex<<dx,dy,dz;
        Eigen::Array3i node_index = bottom_index+dindex;
        if(in_bounds(mesh,node_index)){
          Eigen::Array3d node_position = node_index.cast<double>()*h;
          Eigen::Array3d distance = mp.position.array()-node_position;
          Eigen::Array3d l = Eigen::Array3d::Constant(h/2);
          Eigen::Array3d weights = Eigen::Array3d::Zero();
          Eigen::Array3d grads = Eigen::Array3d::Zero();
          const int nd = 2;
          for(int i = 0;i < nd;++i){
            weights[i] = shape_gimp(distance[i],l[i],h);
            grads[i] = shape_gimp_dsvp(distance[i],l[i],h);
          }
          double svp = weights[0]*weights[1];//*weights[2];
          grads[0] = grads[0]*weights[1];//*weights[2];
          grads[1] = grads[1]*weights[0];//*weights[2];
          grads[2] = 0.0;//grads[2]*weights[0]*weights[1];
          func(mesh.nodes[position_to_index(mesh,node_index)],svp,grads.matrix());
        }
      }
    }
  }
}


// double shape_linear(double distance, double h){
//   return 1.0 - std::abs(distance/h);
// }

// double shape_linear_dsvp(double distance, double h){
//   return -std::copysign(1.0,distance)/h;
// }

void make_nc(Sim & sim, mp & mp){
  Mesh & m = *sim.mesh;
  mp.nc.clear();
  iterate_over_shape_gimp(*sim.mesh,mp,
                          [&](node &n, double svp, Vector grads){
                            nodecache nc;
                            nc.mp = mp.index;
                            nc.svp = svp;
                            nc.grads = grads;
                            nc.node = n.index;
                            nc.nc_index = sim.global_nodecache.size();
                            mp.nc.push_back(nc);
                            {
                              std::lock_guard<std::mutex> guard((*m.node_locks)[nc.node]);
                              m.nodes[nc.node].nc.push_back(nc);
                            }
                            sim.global_nodecache.push_back(nc);
  });

  // Eigen::Array3i bottom_index = (mp.position.array()/m.h).floor().cast<int>();
  // const double h = m.h;
  // for(int dx = 0;dx < 2;++dx){
  //   for(int dy = 0;dy < 2;++dy){
  //     for(int dz = 0;dz < 2;++dz){
  //       Eigen::Array3i dindex;
  //       dindex<<dx,dy,dz;
  //       Eigen::Array3i node_index = bottom_index+dindex;
  //       if(in_bounds(m,node_index)){
  //         Eigen::Array3d node_position = node_index.cast<double>()*m.h;
  //         Eigen::Array3d distance = mp.position.array()-node_position;
  //         Eigen::Array3d weights = distance.matrix().unaryExpr([&](double x){return shape_linear(x,h);});
  //         Eigen::Array3d grads = distance.matrix().unaryExpr([&](double x){return shape_linear_dsvp(x,h);});
  //         nodecache nc;
  //         nc.mp = mp.index;
  //         nc.svp = weights[0]*weights[1];//*weights[2];
  //         nc.grads[0] = grads[0]*weights[1];//*weights[2];
  //         nc.grads[1] = grads[1]*weights[0];//*weights[2];
  //         nc.grads[2] = 0.0;//grads[2]*weights[0]*weights[1];
  //         nc.node = position_to_index(m,node_index);
  //         nc.nc_index = sim.global_nodecache.size();

  //         mp.nc.push_back(nc);
  //         {
  //           std::lock_guard<std::mutex> guard((*m.node_locks)[nc.node]);
  //           m.nodes[nc.node].nc.push_back(nc);
  //         }

  //         sim.global_nodecache.push_back(nc);
  //       }
  //     }
  //   }
  // }
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


void make_nc_mp_ind(Sim &sim, mp & mp){
  int start_index = sim.global_nodecache_mp_indirection.size();
  mp.snc.index = start_index;
  mp.snc.count = mp.nc.size();
  for(auto &nc : mp.nc){
    sim.global_nodecache_mp_indirection.push_back(nc.nc_index);
  }
}

void make_nc_node_ind(Sim &sim, node & node){
  int start_index = sim.global_nodecache_node_indirection.size();
  node.snc.index = start_index;
  node.snc.count = node.nc.size();
  for(auto &nc : node.nc){
    sim.global_nodecache_node_indirection.push_back(nc.nc_index);
  }
}

void setup_svp(Sim &sim){
  for(auto & n : sim.mesh->nodes){
    n.nc.clear();
  }
  sim.global_nodecache.clear();
  for(int i = 0;i < sim.mps.size();++i){
    make_nc(sim,sim.mps[i]);
  }
  for(auto &mp : sim.mps){
    make_nc_mp_ind(sim,mp);
  }
  for(auto &node : sim.mesh->nodes){
    make_nc_node_ind(sim,node);
  }
}

double density = 100;
void setup_mps(Sim & sim){
#pragma omp parallel for
  for(int i = 0;i < sim.mps.size();++i){
    auto & mp = sim.mps[i];
    mp.index = i;
    // make_nc(*sim.mesh,mp);
    mp.set_elastic(1e4,0.2);
    mp.mass = density * mp.volume;
  }
  setup_svp(sim);
}

void update_particles(Sim &sim){
#pragma omp parallel for
  for(int i = 0;i < sim.mps.size();++i){
    auto & mp = sim.mps[i];
    mp.f_n = mp.f;
    mp.strain_n = mp.strain;
    mp.position += mp.disp_inc;
    mp.disp_inc *= 0;
    mp.volume_n = mp.volume;
    mp.nc.clear();
    // calculate_strain(sim,mp);
    // mp.stress = mp.de * mp.strain;
    // mp.stress = mp.stress/mp.f.determinant();
  }
}

void p2g(Sim & sim){
  sim.mesh->mass *= 0.0;
  Eigen::Matrix<double,Eigen::Dynamic,3> volume = sim.mesh->mass;
#pragma omp parallel for
  for(int i = 0;i < sim.mps.size();++i){
    auto & mp = sim.mps[i];
    iterate_over_neighbours(sim,mp,
                            [&](auto & nc){
                              std::lock_guard<std::mutex> guard((*sim.mesh->node_locks)[nc.node]);
                              // sim.mesh->nodes[nc.node] = true;
                              sim.mesh->mass.row(nc.node) += Vector::Constant(mp.mass * nc.svp).transpose();
                              // sim.mesh->mass.row(nc.node) += Vector::Constant(mp.p_mod * mp.volume * nc.svp).transpose();
                              // volume.row(nc.node) += Vector::Constant(mp.volume * nc.svp).transpose();
                                // sim.mesh->mass.row(nc.node) += Vector::Constant(mp.mass * nc.svp).transpose();
                              });
    }
  // auto & fds = sim.mesh->fds;
  // sim.mesh->mass(fds,Eigen::all) = (sim.mesh->mass(fds,Eigen::all).array() / volume(fds,Eigen::all).array());
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

void update_stress(Sim & sim){
#pragma omp parallel for
  for(int i = 0;i < sim.mps.size();++i){
    auto & mp = sim.mps[i];
    calculate_strain(sim,mp);
    mp.stress = mp.de * mp.strain;
    mp.stress = mp.stress/mp.f.determinant();
  }
}


void p2g_node_force(Sim & sim, node & node){
  iterate_over_node_neighbours(sim,node,[&](auto & n){
    auto & mp = sim.mps[n.mp];
    Dsvp dsvp = assemble_dsvp(mp.df.inverse()*n.grads);
    sim.mesh->force.row(n.node) += (-1.0*mp.volume*(dsvp.transpose() * mp.stress).transpose()) + (mp.mass * n.svp * sim.mesh->gravity.transpose());
  });
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
  auto & fds = sim.mesh->fds;
  auto &m = *sim.mesh;
  m.force += (m.mass.array()*sim.mesh->velocity.array() * (-1.0 * sim.damping_factor)).matrix();
  m.velocity(fds,Eigen::all) += sim.dt*(m.force(fds,Eigen::all).array()/m.mass(fds,Eigen::all).array()).matrix();
  // m.velocity += sim.dt*(m.force.array()/m.mass.array()).matrix();
  m.displacement += sim.dt*m.velocity;
}
void reset_force(Sim& sim){
  sim.mesh->force *= 0.0;
  sim.mesh->force_int *= 0.0;
  sim.mesh->force_ext *= 0.0;
}
void compute_residuals(){
}

void apply_bcs(Sim & sim){
  Mesh &mesh = *sim.mesh;
  mesh.displacement = (mesh.displacement.array() * mesh.bcs.array()).matrix();
  mesh.force = (mesh.force.array() * mesh.bcs.array()).matrix();
  mesh.velocity = (mesh.velocity.array() * mesh.bcs.array()).matrix();
  mesh.force_int = (mesh.force_int.array() * mesh.bcs.array()).matrix();
  mesh.force_ext = (mesh.force_ext.array() * mesh.bcs.array()).matrix();
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
    Eigen::Vector3d disp_pos = mp.position;// + mp.disp_inc;
    myfile << disp_pos[0] << " " << disp_pos[1] << " " << disp_pos[2] << "\n";
  }
  myfile << "POINT_DATA " << mps.size() << "\n";
  myfile << "SCALARS s_xx FLOAT 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(auto & mp : mps){
    myfile << mp.stress[0] << "\n";
  }
  myfile << "SCALARS s_yy FLOAT 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(auto & mp : mps){
    myfile << mp.stress[1] << "\n";
  }
  myfile << "SCALARS s_xy FLOAT 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(auto & mp : mps){
    myfile << mp.stress[5] << "\n";
  }
  // myfile << "SCALARS e_yy FLOAT 1\n";
  // myfile << "LOOKUP_TABLE default\n";
  // for(auto & mp : mps){
  //   myfile << mp.strain[1] << "\n";
  // }
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
  myfile << "VECTORS disp FLOAT\n";
  for(auto & node : nodes){
    Eigen::Vector3d f = mesh.displacement.row(node.index).transpose();
    myfile << f[0] << " " << f[1] << " " << f[2] << "\n";
  }


  myfile.close();
}


void setup_sim(Sim &sim){
  sim.mesh->compact_mesh();
}

void save(Sim &sim, int index){
  std::ostringstream ss;
  ss << "./output/test_"<<index<<".vtk";
  save_vtk(ss.str(),sim.mps);
  std::ostringstream ssn;
  ssn << "./output/test_n_"<<index<<".vtk";
  save_vtk_nodes(ssn.str(),*sim.mesh);
}

int total_iters = 0;
void solve(Sim & sim,int step){
  const double oobf_crit = 1e-3;
  sim.mesh->reset();
  setup_svp(sim);
  p2g(sim);
  setup_sim(sim);
  sim.dt = 0.25 * sim.estimate_elastic_dt();
  sim.damping_factor = 0.0;
  const int iters = 100;
  const int substeps = 100;
  double oobf = oobf_crit;
  double ke_0,ke_1 = 0;
  for(int i = 0;(i < iters) && (oobf >= oobf_crit);++i){
    std::cout<<"Solve step: " << i <<"\n";
    for(int substep = 0;substep < substeps;++substep){
      total_iters++;
      reset_force(sim);
      update_stress(sim);
      // p2g_scatter_force(sim);
      p2g_gather_force(sim);
      integrate(sim);
      apply_bcs(sim);
      // time += sim.dt;
      ke_1 = ke_0;
      ke_0 = sim.calculate_ke();
      if(ke_0 < ke_1){
        sim.mesh->velocity *= 0;
      }
    }
    oobf = sim.calculate_oobf();
    std::cout<<"OOBF: " << oobf << " - E: " << sim.calculate_ke() <<"\n";
    // save(sim,i);
  }
  g2p(sim);
  update_particles(sim);
  if(oobf >= oobf_crit){
    save(sim,step);
    std::cout << "Failed to converge!\n";
    std::abort();
  }
}


int main(int argc,char ** args) {
  typedef std::chrono::high_resolution_clock Clock;
  const double h = 10;
  const double size = 100.0;
  Eigen::Vector3d block_size = (Eigen::Vector3d()<<size,size,1.0).finished();
  Sim sim(h,std::round(block_size[0]/h)+1,std::round(block_size[1]/h)+1,+1);
  double mps_per_cell = 3;
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
  sim.setup_mass_filter(density,1e-4);
  save(sim,0);
  // save_vtk("./test_0.vtk",mps);
  // save_vtk_nodes("./test_n_0.vtk",m);

  std::cout << "Node count: " << sim.mesh->node_count << "\n";
  std::cout << "NDOFs: " << (3 * sim.mesh->node_count) << "\n";
  std::cout << "MPs: " << mp_count << "\n";


  int iters = 1000;
  int substeps = 1;
  std::cout << "Iters: " << iters << "\n";
  std::cout << "Substeps: " << substeps << "\n";
  // std::cout<<sim.mesh->mass;
  auto t1 = Clock::now();
  // gpu::setup_sim(sim);
  p2g(sim);
  setup_sim(sim);
  sim.dt = 0.25 * sim.estimate_elastic_dt();
  // sim.dt = 0.25;
  sim.damping_factor = 0.0;
  // sim.damping_factor = 0.1*sim.estimate_critical_damping();
  std::cout<<"Estimated timestep: " << sim.dt<<"\n";
  std::cout<<"Estimated damping: " << sim.damping_factor<<"\n";
  // sim.dt = 1e-3;
  sim.mesh->gravity *= 0.01;
  // sim.damping_factor = 0;
  // sim.mps[0].stress[1] = 1;
  double time = 0;
  double ke_0,ke_1 = 0;
  for(int i = 0;i < iters;++i){
    std::cout<<"Step: " << i <<"\n";
    for(int substep = 0;substep < substeps;++substep){
      reset_force(sim);
      update_stress(sim);
      p2g_scatter_force(sim);
      // p2g_gather_force(sim);
      integrate(sim);
      apply_bcs(sim);
      time += sim.dt;
      ke_1 = ke_0;
      ke_0 = sim.calculate_ke();
      if(ke_0 < ke_1){
        sim.mesh->velocity *= 0;
      }
    }
    // gpu::sync_sim(sim);
    // g2p(sim);
    save(sim,i+1);
    std::cout<<"OOBF: " << sim.calculate_oobf() << " - E: " << sim.calculate_ke() <<"\n";

    // std::cout <<"\n" << sim.mps[0].f<<"\n";
    // std::cout <<"\n" << sim.mps[0].df<<"\n";
  }
  // gpu::close_sim(sim);
  auto t2 = Clock::now();
  std::cout << "Took: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " seconds \n";
  std::cout << "MP time: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count()/(iters *substeps* mp_count)  << " su/seconds \n";
  return 0;
}



// int main(int argc,char ** args) {
//   typedef std::chrono::high_resolution_clock Clock;
//   const double h = 1;
//   const double dsize = 1000.0;
//   const double size = 100.0;
//   Eigen::Vector3d domain_size = (Eigen::Vector3d()<<dsize,dsize,1.0).finished();
//   Eigen::Vector3d block_size = (Eigen::Vector3d()<<size,size,1.0).finished();
//   Sim sim(h,std::round(domain_size[0]/h)+1,std::round(domain_size[1]/h)+1,+1);
//   sim.setup_mass_filter(density,1e-4);
//   double mps_per_cell = 2;
//   double mp_res = (h/(mps_per_cell));
//   double offset = mp_res*0.5;//(h/(1+mps_per_cell));
//   int x_count = std::floor(block_size[0]/mp_res);
//   int y_count = std::floor(block_size[1]/mp_res);
//   // int x_count = 2;
//   // int y_count = 2;
//   int mp_count = x_count * y_count;
//   sim.mps = std::vector<mp>(mp_count);
//   int i = 0;
//   for(int x = 0; x < x_count;++x){
//     for(int y = 0; y < y_count;++y){
//       sim.mps[i].position[0] = (x * mp_res) + offset;
//       sim.mps[i].position[1] = (y * mp_res) + offset;
//       sim.mps[i].volume = mp_res * mp_res;
//       sim.mps[i].volume_n = sim.mps[i].volume;
//       i++;
//     }
//   }
//   std::cout<<"Setup\n";
//   setup_mps(sim);
//   save(sim,0);
//   // save_vtk("./test_0.vtk",mps);
//   // save_vtk_nodes("./test_n_0.vtk",m);

//   std::cout << "Node count: " << sim.mesh->node_count << "\n";
//   std::cout << "NDOFs: " << (3 * sim.mesh->node_count) << "\n";
//   std::cout << "MPs: " << mp_count << "\n";

//   int lstps = 100;
//   auto t1 = Clock::now();
//   Vector gravity;
//   gravity << 0.0,-9.8,0.0;
//   // gravity *= 1;
//   for(int i = 1;i <= lstps;++i){
//     std::cout<<"Step: " << i <<"\n";
//     sim.mesh->gravity = gravity * ((double)i / (double)lstps);
//     solve(sim,i);
//     save(sim,i);
//   }
//   auto t2 = Clock::now();
//   std::cout << "Took: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " seconds \n";
//   std::cout << "MP time: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count()/(total_iters* mp_count)  << " su/seconds \n";
//   return 0;
// }

// int main(int argc,char ** args) {
//   double h = 1;
//   double l = h/2;
//   for(double x = -2;x<=2;x+=0.1){
//     std::cout <<x<<","<< shape_gimp(x,l,h)<<"," << shape_gimp_dsvp(x,l,h)<<"\n";
//   }
// }
