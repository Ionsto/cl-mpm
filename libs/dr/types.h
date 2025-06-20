#pragma once
#include <Eigen/Dense>
#include "../src/constitutive.h"
#include <mutex>
#include <memory>
#include <vector>

using Voigt = Eigen::Matrix<double,6,1>;
using Vector = Eigen::Matrix<double,3,1>;
using Dsvp = Eigen::Matrix<double,6,3>;
using Vec_mass = Eigen::Matrix<double,Eigen::Dynamic,1>;
using Vec_vel = Eigen::Matrix<double,Eigen::Dynamic,3>;

struct nodecache {
public:
  int node = 0;
  int mp = 0;
  double svp = 0;
  Vector grads = Vector::Zero();
};

struct gpu_mp{
  int index = 0;
  Vector position = Vector::Zero();
  Vector disp_inc = Vector::Zero();
  Voigt strain_n = Voigt::Zero();
  Voigt strain = Voigt::Zero();
  Voigt stress = Voigt::Zero();
  Eigen::Matrix3d f_n = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d f = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d df = Eigen::Matrix3d::Identity();
  double volume_n = 1;
  double volume = 1;
};

struct mp {
  int index = 0;
  Vector position = Vector::Zero();
  Vector disp_inc = Vector::Zero();
  Voigt strain_n = Voigt::Zero();
  Voigt strain = Voigt::Zero();
  Voigt stress = Voigt::Zero();
  double mass = 1;
  double e = 1;
  double nu = 1;
  double p_mod = 1;
  Eigen::Matrix<double,6,6> de = Eigen::Matrix<double,6,6>::Ones();
  Eigen::Matrix3d f_n = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d f = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d df = Eigen::Matrix3d::Identity();
  double volume_n = 1;
  double volume = 1;
  std::vector<nodecache> nc;
  mp() = default;
  void set_elastic(double e, double nu){
    de = constitutive::AssembleDE(e,nu);
    this->e = e;
    this->nu = e;
    p_mod = e/((1+nu)*(1-nu));
  }
};


struct node{
  int index = 0;
  std::vector<nodecache> nc;
  Eigen::Matrix<double,3,1> position;
};

struct super_node_cache{
  int index = 0;
  int count = 0;
};


struct Mesh{
public:
  double h;
  int node_count;
  Eigen::Array3i element_count;
  Vector gravity;
  Eigen::Matrix<double,Eigen::Dynamic,3> mass;
  Eigen::Matrix<double,Eigen::Dynamic,3> velocity;
  Eigen::Matrix<double,Eigen::Dynamic,3> force;
  Eigen::Matrix<double,Eigen::Dynamic,3> displacement;
  Eigen::Matrix<double,Eigen::Dynamic,3> bcs;

  std::vector<int> compact_bcs_index;
  std::vector<Eigen::Matrix<double,3,1>> compact_bcs_values;

  std::vector<int> index_indirection;

  std::vector<node> nodes;
  std::unique_ptr<std::vector<std::mutex>> node_locks;

  int location_to_index(Eigen::Vector3i index){
    return index[0] + (index[1] * element_count[0]) + (index[2] * element_count[0] * element_count[1]);
  }
  Mesh(double h_0, int element_count_x, int element_count_y, int element_count_z){
    h = h_0;
    element_count << element_count_x,element_count_y,element_count_z;
    node_count = element_count_x * element_count_y * element_count_z;
    node_locks = std::make_unique<std::vector<std::mutex>>(node_count);
    std::cout<<"Make locks: " << node_locks->size();
    gravity << 0.0,-9.8,0.0;
    mass = Vec_vel::Zero(node_count,3);
    force = Vec_vel::Zero(node_count,3);
    bcs = Vec_vel::Ones(node_count,3);
    velocity = Vec_vel::Zero(node_count,3);
    displacement = Vec_vel::Zero(node_count,3);
    nodes = std::vector<node>(node_count);

    for(int x = 0;x < element_count_x;++x){
      for(int y = 0;y < element_count_y;++y){
        for(int z = 0;z < element_count_z;++z){
          int i = x + (y*element_count_x) + (z*element_count_x*element_count_y);
          nodes[i].index = i;
          nodes[i].position[0] = x*h;
          nodes[i].position[1] = y*h;
          nodes[i].position[2] = z*h;
        }
      }
    }

    for(int x = 0;x< element_count_x;++x){
      bcs(location_to_index((Eigen::Vector3i()<<x,0,0).finished()), 1) = 0.0;
      bcs(location_to_index((Eigen::Vector3i()<<0,x,0).finished()), 0) = 0.0;
      // bcs(location_to_index((Eigen::Vector3i()<<element_count_x-1,x,0).finished()), 0) = 0.0;
      // add_bc((Eigen::Vector3i()<<0,x,0).finished(),(Eigen::Vector3i()<<0,1,1).finished());
      // add_bc((Eigen::Vector3i()<<0,x,0).finished(),(Eigen::Vector3i()<<0,1,1).finished());
    }
  }
  void add_bc(Eigen::Vector3i index, Eigen::Vector3d values){
    compact_bcs_index.push_back(location_to_index(index));
    compact_bcs_values.push_back(values);
  }
  void compact_mesh(){
  }
};


struct Sim{
  double dt = 1;
  double damping_factor = 0;
  std::unique_ptr<Mesh> mesh;
  std::vector<mp> mps;

  std::vector<super_node_cache> global_nodecache_mp_indirection;
  std::vector<super_node_cache> global_nodecache_node_indirection;
  std::vector<nodecache> global_nodecache;

  Sim(double h_0, int element_count_x, int element_count_y, int element_count_z){
    mesh = std::make_unique<Mesh>(h_0,element_count_x,element_count_y,element_count_z);
  }

  double estimate_elastic_dt_mp(mp & mp){
    return mesh->h * std::sqrt(mp.mass/(mp.volume*mp.p_mod));
  }
  double estimate_elastic_dt(){
    if(mps.size()==0){
      return dt;
    }
    else
    {
      return estimate_elastic_dt_mp(*std::min_element(mps.begin(),mps.end(),
                                                   [&](auto l,auto r){
                                                     return this->estimate_elastic_dt_mp(l) < this->estimate_elastic_dt_mp(r);
                                                   }
                                                   ));
    }
  }

  double estimate_critical_damping_mp(mp & mp){
    return 3.14 * 2 * std::sqrt(mp.e / (mesh->h*mesh->h * mp.mass/(mp.volume)));
  }

  double estimate_critical_damping(){
    if(mps.size()==0){
      return dt;
    }
    else
      {
        double damping = estimate_critical_damping_mp(mps[0]);
        for(auto &mp :mps){
          double d = estimate_critical_damping_mp(mp);
          if (d < damping){
            damping = d;
          }
        }
        return damping;
      }
  }

  // inline
  // template <typename F>
  // void iterate_over_neighbours(mp & mp, F&& func){
  //   for (auto &n : mp.nc){
  //     func(mesh->nodes[n.node],n);
  //   }
  // }
};
