#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "src/utils.h"
#include "src/constitutive.h"
#include <chrono>
#include <mutex>
#include <memory>

using Voigt = Eigen::Matrix<double,6,1>;
using Vector = Eigen::Matrix<double,3,1>;
using Dsvp = Eigen::Matrix<double,6,3>;

struct nodecache {
public:
  int node = 0;
  int mp = 0;
  double svp = 0;
  Vector grads = Vector::Zero();
  Dsvp dsvp;
};

struct mp {
  int index = 0;
  Vector position = Vector::Zero();
  Voigt strain_n = Voigt::Zero();
  Voigt strain = Voigt::Zero();
  Voigt stress = Voigt::Zero();
  double mass = 1;
  Eigen::Matrix<double,6,6> de = Eigen::Matrix<double,6,6>::Ones();
  Eigen::Matrix3d f_n = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d f = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d df = Eigen::Matrix3d::Identity();
  double volume_n = 1;
  double volume = 1;
  std::vector<nodecache> nc;
};

using Vec_mass = Eigen::Matrix<double,Eigen::Dynamic,1>;
using Vec_vel = Eigen::Matrix<double,Eigen::Dynamic,3>;

struct node{
  int index = 0;
  std::vector<nodecache> nc;
  Eigen::Matrix<double,3,1> position;
};



struct mesh{
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
  std::vector<node> nodes;
  std::unique_ptr<std::vector<std::mutex>> node_locks;

  int location_to_index(Eigen::Vector3i index){
    return index[0] + (index[1] * element_count[0]) + (index[2] * element_count[0] * element_count[1]);
  }
  mesh(double h_0, int element_count_x, int element_count_y, int element_count_z){
    h = h_0;
    element_count << element_count_x,element_count_y,element_count_z;
    node_count = element_count_x * element_count_y * element_count_z;
    // node_locks = std::vector<std::unique_ptr<std::mutex>>(node_count);
    node_locks = std::make_unique<std::vector<std::mutex>>(node_count);
    std::cout<<"Make locks: " << node_locks->size();
    gravity << 0.0,-9.8,0.0;
    mass = Vec_vel::Ones(node_count,3);
    force = Vec_vel::Zero(node_count,3);
    bcs = Vec_vel::Ones(node_count,3);
    velocity = Vec_vel::Zero(node_count,3);
    displacement = Vec_vel::Zero(node_count,3);
    nodes = std::vector<node>(node_count);
    for(int x = 0;x< element_count_x;++x){
      for(int y = 0;y< element_count_y;++y){
        for(int z = 0;z< element_count_z;++z){
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
    }
  }
};

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

void combi_assemble_dsvp(Eigen::Matrix3d df, Vector disp, Vector grads){
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


void calculate_strain(mesh & mesh, mp & mp){
  // mp.strain_n = mp.strain;
  mp.df = Eigen::Matrix3d::Identity();
  for(auto &n : mp.nc){
    combi_assemble_dsvp(mp.df,mesh.displacement.row(n.node).transpose(),n.grads);
  }
  mp.f_n = mp.df * mp.f;
  mp.volume = mp.volume_n * mp.df.determinant();
  mp.strain = log_strain_update(mp.strain_n,mp.df);
}

void calculate_force(mesh & mesh, mp & mp){
  for(auto &n : mp.nc){
    std::lock_guard<std::mutex> guard((*mesh.node_locks)[n.node]);
    mesh.force.row(n.node) += (mp.mass * n.svp * mesh.gravity.transpose())
     +(mp.volume*(n.dsvp.transpose() * mp.stress).transpose());
  }
}



int position_to_index(mesh& m, Eigen::Array3i index){
  return index[0] + (index[1] * m.element_count[0]) + (index[2] * m.element_count[0] * m.element_count[1]);
}

int in_bounds(mesh& m, Eigen::Array3i index){
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
    return -1.0;
  } else if(distance < 0){
    return 1.0;
  }
  return 0.0;
}

void make_nc(mesh & m, mp & mp){
  Eigen::Array3i bottom_index = (mp.position.array()/m.h).floor().cast<int>();
  const double h = m.h;
  for(int dx = 0;dx < 2;++dx){
    for(int dy = 0;dy < 2;++dy){
      for(int dz = 0;dz < 2;++dz){
        Eigen::Array3i dindex;
        dindex<<dx,dy,dz;
        Eigen::Array3i node_index = bottom_index+dindex;
        if(in_bounds(m,node_index)){
          Eigen::Array3d node_position = bottom_index.cast<double>()*m.h;
          Eigen::Array3d distance = mp.position.array()-node_position;
          Eigen::Array3d weights = distance;
          weights.matrix().unaryExpr([&](double x){return shape_linear(x,h);});
          Eigen::Array3d grads = distance;
          grads.matrix().unaryExpr([&](double x){return shape_linear_dsvp(x,h);});
          nodecache nc;
          nc.mp = mp.index;
          nc.svp = weights[0]*weights[1]   ;//*weights[2];
          nc.grads[0] = grads[0]*weights[1];//*weights[2];
          nc.grads[1] = grads[1]*weights[0];//*weights[2];
          nc.grads[2] = 0.0;//grads[2]*weights[0]*weights[1];
          nc.node = position_to_index(m,node_index);
          nc.dsvp = assemble_dsvp(nc.grads);
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


void setup_mps(std::vector<mp> & mps, mesh & mesh){
#pragma omp parallel for
  for(int i = 0;i < mps.size();++i){
    auto & mp = mps[i];
    mp.index = i;
    make_nc(mesh,mp);
  }
}

void update_mps(std::vector<mp> & mps, mesh & mesh, double dt){
#pragma omp parallel for
  for(int i = 0;i < mps.size();++i){
    auto & mp = mps[i];
    calculate_strain(mesh,mp);
    mp.stress = mp.de * mp.strain;
  }
}


void p2g_node_force(mesh & mesh,std::vector<mp> &mps, node & node){
  for(auto &n : node.nc){
    auto & mp = mps[n.mp];
    mesh.force.row(n.node) += (mp.volume*(n.dsvp.transpose() * mp.stress).transpose()) + (mp.mass * n.svp * mesh.gravity.transpose());
  }
}

void p2g_scatter_force(mesh & mesh, std::vector<mp> & mps){
#pragma omp parallel for
  for(int i = 0;i < mps.size();++i){
    auto & mp = mps[i];
    calculate_force(mesh,mp);
  }
}

void p2g_gather_force(mesh & mesh, std::vector<mp> & mps){
#pragma omp parallel for
  for(int i = 0;i < mesh.nodes.size();++i){
    auto & node = mesh.nodes[i];
    p2g_node_force(mesh,mps,node);
  }
}

void apply_bcs(mesh & mesh){
  mesh.displacement = (mesh.displacement.array() * mesh.bcs.array()).matrix();
  // mesh.velocity *= mesh.bcs;
  // mesh.force *= mesh.bcs;
// #pragma omp parallel for
//   for(int i = 0;i < mesh.bcs.size();++i){
//     for(int d = 0;d < 3;++d){
//       if(mesh.bcs(i,d) == 1){
//         mesh.displacement(i,d) = 0.0;
//         mesh.velocity(i,d) = 0.0;
//         mesh.force(i,d) = 0.0;
//       }
//     }
//   }
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
    myfile << mp.position[0] << " " << mp.position[1] << " " << mp.position[2] << "\n";
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

void save_vtk_nodes(std::string filename,mesh &mesh){
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
  // myfile << "SCALARS FORCE FLOAT 1\n";
  // myfile << "LOOKUP_TABLE default\n";
  // for(auto & node : nodes){
  //   Eigen::Vector3d f = mesh.force.row(node.index).transpose();
  //   myfile << f[0] << " " << f[1] << " " << f[2] << "\n";
  // }
  myfile << "VECTORS FORCE FLOAT\n";
  // myfile << "LOOKUP_TABLE default\n";
  for(auto & node : nodes){
    Eigen::Vector3d f = mesh.force.row(node.index).transpose();
    myfile << f[0] << " " << f[1] << " " << f[2] << "\n";
  }

  myfile.close();
}


int main(int argc,char ** args) {
  typedef std::chrono::high_resolution_clock Clock;
  const double h = 1;
  mesh m(h,100,100,1);
  Eigen::Vector3d block_size = (Eigen::Vector3d()<<100.0,100.0,0.0).finished();
  double mps_per_cell = 2;
  double mp_res = (h/(mps_per_cell));
  double offset = (h/(1+mps_per_cell));
  int x_count = std::round(block_size[0]/mp_res);
  int y_count = std::round(block_size[1]/mp_res);
  // int x_count = 2;
  // int y_count = 2;
  int mp_count = x_count * y_count;
  std::vector<mp> mps(mp_count);
  int i = 0;
  for(int x = 0; x < x_count;++x){
    for(int y = 0; y < y_count;++y){
      mps[i].position[0] = (x * mp_res) + offset;
      mps[i].position[1] = (y * mp_res) + offset;
      i++;
    }
  }
  std::cout<<"Setup\n";
  setup_mps(mps,m);

  // for(auto &mp :mps){
  //   std::cout << mp.nc[0].dsvp<<"\n";
  // }
  save_vtk("./test_0.vtk",mps);
  save_vtk_nodes("./test_n_0.vtk",m);

  std::cout << "Node count: " << m.node_count << "\n";
  std::cout << "NDOFs: " << (3 * m.node_count) << "\n";
  std::cout << "MPs: " << mp_count << "\n";


  int iters = 100;
  std::cout << "Iters: " << iters << "\n";
  double dt = 1e-2;
  // init_mps(mps,m,dt);
  auto t1 = Clock::now();
  for(int i = 0;i < iters;++i){
    m.force.array() *= 0.0;
    update_mps(mps,m,dt);
    p2g_scatter_force(m,mps);
    // p2g_gather_force(m,mps);
    m.velocity += dt*(m.force.array()/m.mass.array()).matrix();
    m.displacement += dt*m.velocity;
    apply_bcs(m);
    std::ostringstream ss;
    ss << "./test_"<<i<<".vtk";
    save_vtk(ss.str(),mps);
    std::ostringstream ssn;
    ssn << "./test_n_"<<i<<".vtk";
    save_vtk_nodes(ssn.str(),m);
  }
  // save_vtk("./test_1.vtk",mps);
  auto t2 = Clock::now();
  std::cout << "Took: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " seconds \n";
  // std::cout << "Throughput: " << iters / (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count() << " su/seconds \n";
  std::cout << "MP time: " << (std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1)).count()/(iters * mp_count)  << " su/seconds \n";
  return 0;
}
