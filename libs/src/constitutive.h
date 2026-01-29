#pragma once
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include "utils.h"
#include <iostream>

namespace constitutive{
  using namespace utils;
  inline
    Eigen::Matrix<double,6,6> AssembleDE(double E, double nu){
      const double nm = (1.0 - nu);
      const double gf = 0.5*(1.0 - (2 * nu));
      Eigen::Matrix<double,6,6> De;
      De<<
        nm,nu,nu,0.0,0.0,0.0,
        nu,nm,nu,0.0,0.0,0.0,
        nu,nu,nm,0.0,0.0,0.0,
        0.0,0.0,0.0,gf,0.0,0.0,
        0.0,0.0,0.0,0.0,gf,0.0,
        0.0,0.0,0.0,0.0,0.0,gf;
      De *= E/((1+nu)*(1-(2*nu)));
      return De;
    }

  inline
    int rotate(int i){
      return (i+1)%3;
    }
  inline
    Eigen::Matrix<double,6,6> AssembleQMatrix(Eigen::Matrix<double,3,3> eigen_vectors){
      Eigen::Matrix<double,6,6> Q;
      for(int j = 0;j < 3;++j){
        for(int i = 0;i < 3;++i){
          Q(i,j)   = eigen_vectors(i,j) * eigen_vectors(i,j);
          Q(i+3,j) = eigen_vectors(i,j) * eigen_vectors(rotate(i),j);
        }
      }
      // Q.block(0,0,3,3) = eigen_vectors.array().square().matrix();
      for(int j = 0;j < 3;++j){
        for(int i = 0;i < 3;++i){
          Q(i,j+3) = 2*eigen_vectors(i,j) * eigen_vectors(i,rotate(j));
          Q(i+3,j+3) =
            (eigen_vectors(i,j) * eigen_vectors(rotate(i),rotate(j)))
            +
            eigen_vectors(rotate(i),j) * eigen_vectors(i,rotate(j));
          /* Q(i+3,j+3) = */
          /*   (eigen_vectors(i,j) * eigen_vectors((i+1)%3,(j+1)%3)) */
          /*   + */
          /*   eigen_vectors((i+1)%3,j) * eigen_vectors(i,(j+1)%3); */
        }
      }
      return Q.transpose();
    }

    /*
      Because the performance of MAGICL is terrible for allocations it is actually highly
      benificial to run complicated constitutive models in c++ leaning on the optimised
      eigen library to do all linalg operations
    */

  inline
    Eigen::Matrix<double,6,1> DruckerPrager(Eigen::Matrix<double,6,1> elastic_strain,
                        double E, double nu, double phi, double psi, double c) {

      const double alfa = -std::tan(phi);
      const double bta = -std::tan(psi);
      const double xsic = std::sqrt(3)*(1.0 / std::tan(phi))*c;
      // Eigen::Matrix<double,3,3> Ce = (Eigen::Matrix<double,3,3>()<<
      //                                 E,-nu,-nu,
      //                                 -nu,E,-nu,
      //                                 -nu,-nu,E).finished()/E;
      Eigen::Matrix<double,3,3> De3 =
        (E/((1+nu) * (1-(2*nu))))*
        (((1-(2*nu))*Eigen::Matrix<double,3,3>::Identity()) +
        Eigen::Matrix<double,3,3>::Constant(nu));
      Eigen::Matrix<double,3,3> Ce = (Eigen::Matrix<double,3,3>()<<
                                      1,-nu,-nu,
                                      -nu,1,-nu,
                                      -nu,-nu,1).finished()/E;
      /* Eigen::Matrix<double,3,3> Ce = De3.inverse(); */

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
      //std::cout <<"eigenvalues\n"<<eigen_values<<"\n";
      Eigen::Matrix<double,3,3> eigen_vectors = eigensolver.eigenvectors().rowwise().reverse();
      //std::cout <<"eigenvectors\n"<<eigen_vectors<<"\n";
      Eigen::Matrix<double,3,1> sig = ((De3 * eigen_values).array() - (xsic / std::sqrt(3))).matrix();
      //std::cout<<"De3\n"<<De3<<"\n";
      const double tol = 1e-12;
      double xi = sig.sum()/std::sqrt(3);
      Eigen::Matrix<double,3,1> s = (sig.array() - (xi / std::sqrt(3))).matrix();
      double rho = std::sqrt(s.array().square().sum());
      double f = rho - alfa*xi;
      // std::cout<<"f\n"<<f<<"\n";
      if (f>tol){
        Eigen::Matrix<double,3,1> epsE = Ce * sig;
        Eigen::Matrix<double,3,1> epsEtr = epsE;
        auto Q = AssembleQMatrix(eigen_vectors);
        double fap = 0.0;
        if(bta != 0.0){
          fap = rho * std::sqrt(1+nu) + (xi*std::sqrt(1-(2*nu))/(bta*std::sqrt(1+nu)/std::sqrt(1-2*nu)));
        }
        else{
          fap = rho * std::sqrt(1+nu);
        }
        if(fap<tol)
          {
            //Apex return
            // std::cout<<"xsic:"<<xsic<<"\n";
            sig=Eigen::Matrix<double,3,1>::Constant(xsic/sqrt(3));
          }
        else{
          //Surface return
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
          // std::cout<<"df\n"<<df<<"\n";
          // std::cout<<"dg\n"<<dg<<"\n";
          // std::cout<<"ddg\n"<<ddg<<"\n";
          const int maxit = 4;
          const double tolf = 1e-6;
          for(int iter = 0; (iter < maxit) && ((b.block(0,0,2,1).norm() > tol) || (std::abs(b(3)) > tolf));++iter){
            Eigen::Matrix<double,4,4> A;
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
        //std::cout<<"Sig\n"<<sig<<"\n";
        epsE = Ce*sig;
        Eigen::Matrix<double,3,1> pinc = epsE - epsEtr;
        const double psinc = (1/6) * (std::pow(pinc[0] - pinc[1],2) +
                                            std::pow(pinc[1] - pinc[2],2) +
                                            std::pow(pinc[2] - pinc[0],2));
        //std::cout<<"EpsE\n"<<epsE<<"\n";
        // std::cout<<"Ce\n"<<Ce<<"\n";
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


  inline
    double MC_princ_yield_func(Eigen::Matrix<double,3,1> sig, double phi, double c)
    {
      const double k = (1.0 + std::sin(phi)) / (1.0 - std::sin(phi));
      const double sigc = 2 * c * std::sqrt(k);
      return (k * sig[0]) - (sig[2] + sigc);
    }


  using MohrCoulombReturn = std::tuple<Eigen::Matrix<double,6,1>,float,float,bool,float>;

  inline
    MohrCoulombReturn MohrCoulomb(Eigen::Matrix<double,6,1> elastic_strain, double E, double nu, double phi, double psi, double c) {


      Eigen::Matrix<double,3,3> Ce = (Eigen::Matrix<double,3,3>()<<
                                      1,-nu,-nu,
                                      -nu,1,-nu,
                                      -nu,-nu,1).finished()/E;
      Eigen::Matrix<double,3,3> De3 =
        (E/((1+nu) * (1-(2*nu))))*
        (((1-(2*nu))*Eigen::Matrix<double,3,3>::Identity()) +
         Eigen::Matrix<double,3,3>::Constant(nu));
      /* Eigen::Matrix<double,3,3> Ce = De3.inverse(); */

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(voigt_to_matrix(elastic_strain));
      if (eigensolver.info() != Eigen::Success)
        {
          abort();
        }
      Eigen::Matrix<double,3,1> eigen_values = eigensolver.eigenvalues().reverse();
      Eigen::Matrix<double,3,3> eigen_vectors = eigensolver.eigenvectors().rowwise().reverse();//.array().reverse().matrix();
      Eigen::Matrix<double,3,1> sig = (De3 * eigen_values);
      const double tol = 1e-12;
      double f = MC_princ_yield_func(sig,phi,c);
      if (f>tol){
        // Eigen::Matrix<double,3,1> sigTr = sig;
        const double k = (1+std::sin(phi))/(1-std::sin(phi));
        const double sigc = 2*c*std::sqrt(k);
        const double m = (1+std::sin(psi))/(1-std::sin(psi));

        Eigen::Matrix<double,3,1> siga = Eigen::Matrix<double,3,1>::Constant(sigc / (k - 1)); Eigen::Matrix<double,3,1> epsE = Ce * sig;
        Eigen::Matrix<double,3,1> epsEtr = epsE;
        Eigen::Matrix<double,3,1> r1 = (Eigen::Matrix<double,3,1>() << 1.0, 1.0, k).finished();
        Eigen::Matrix<double,3,1> r2 = (Eigen::Matrix<double,3,1>() << 1.0, k, k).finished();
        Eigen::Matrix<double,3,1> rg1 = (Eigen::Matrix<double,3,1>() << 1.0, 1.0, m).finished();
        Eigen::Matrix<double,3,1> rg2 = (Eigen::Matrix<double,3,1>() << 1.0, m, m).finished();

        Eigen::Matrix<double,3,1> df = (Eigen::Matrix<double,3,1>() << k, 0.0, -1.0).finished();
        Eigen::Matrix<double,3,1> dg = (Eigen::Matrix<double,3,1>() << m, 0.0, -1.0).finished();

        Eigen::Matrix<double,3,1> rp = (De3 * dg) * (1.0 / (dg.transpose() * De3 * df)(0,0));
        const double t1 = (rg1.transpose() * Ce * (sig - siga))[0] / (rg1.transpose() * Ce * r1)[0];
        const double t2 = (rg2.transpose() * Ce * (sig - siga))[0] / (rg2.transpose() * Ce * r2)[0];
        const double f12 = (rp.cross(r1).transpose() * (sig - siga))[0] / (rg2.transpose() * Ce * r2)[0];
        const double f13 = (rp.cross(r2).transpose() * (sig - siga))[0] / (rg2.transpose() * Ce * r2)[0];
        Eigen::Matrix<double,6,6> dep = Eigen::Matrix<double,6,6>::Zero();
        auto Q = AssembleQMatrix(eigen_vectors);
        if((t1 > tol) && (t2 > tol)){
          sig = siga;
        }
        else if ((f12 < tol) && (f13 < tol)){
          sig = siga + (r1 * t1);
          dep.block<3,3>(0,0) = r1*rg1.transpose() / (r1.transpose() * Ce * rg1)[0];
          dep.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Identity()*(E/(2*(1+nu)));
        }
        else if ((f12 > tol) && (f13 > tol)){
          sig = siga + (r2 * t2);
          dep.block<3,3>(0,0) = r2*rg2.transpose() / (r2.transpose() * Ce * rg2)[0];
          dep.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Identity()*(E/(2*(1+nu)));
        }
        else{
          sig = sig - (rp * f);
          dep.block<3,3>(0,0) = De3 - ((De3 * (dg*df.transpose())* De3) / (df.transpose() * De3 * dg)(0,0));
          dep.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Identity()*(E/(2*(1+nu)));
        }
        Eigen::Matrix<double,3,3> T = Eigen::Matrix<double,3,3>::Zero();
        epsE = Ce*sig;
        Eigen::Matrix<double,3,1> sigTr = (De3 * eigen_values);
        if(std::abs(sigTr[0]-sigTr[1])>1e-3){
          T(0,0) = (sig[0]-sig[1])/(sigTr[0]-sigTr[1]);
        }
        if(std::abs(sigTr[1]-sigTr[2])>1e-3){
          T(1,1) = (sig[1]-sig[2])/(sigTr[1]-sigTr[2]);
        }
        if(std::abs(sigTr[0]-sigTr[2])>1e-3){
          T(2,2) = (sig[0]-sig[2])/(sigTr[0]-sigTr[2]);
        }
        dep.block<3,3>(3,3) = T * dep.block<3,3>(3,3);
        dep = Q.transpose() * dep * Q;
        // Eigen::Matrix<double,3,1> pinc = (epsE - epsEtr).reverse();
        Eigen::Matrix<double,3,1> pinc = (epsE - epsEtr);
        const double psinc = std::sqrt(0.5 *
                                      (std::pow(pinc[0] - pinc[1],2) +
                                        std::pow(pinc[1] - pinc[2],2) +
                                        std::pow(pinc[2] - pinc[0],2)));
        Eigen::Matrix<double,6,1> outstrain = swizzle_coombs_voigt(Q.partialPivLu().solve((Eigen::Matrix<double,6,1>()
                                                                              <<
                                                                              epsE[0],
                                                                              epsE[1],
                                                                              epsE[2],
                                                                              0.0,0.0,0.0).finished()));
        double pmod_0 = (((1-nu)*E)/((1 + nu) * (1 - (2 * nu))));
        double pmod = pmod_0*1e-9;
        // for(int i = 0;i < 3;++i){
        //   Eigen::Matrix<double,6,1> n = Eigen::Matrix<double,6,1>::Zero();
        //   n(i) = 1.0;
        //   pmod = std::max(pmod,(n.transpose() * dep * n)(0,0));
        // }
        Eigen::Matrix<double,6,1> n = Eigen::Matrix<double,6,1>::Zero();
        n(0) = 1.0;
        pmod = std::max(pmod,(n.transpose() * dep * n)(0,0));
        return MohrCoulombReturn(outstrain,f,psinc,true,pmod);
      }
      double pmod = (((1-nu)*E)/((1 + nu) * (1 - (2 * nu))));
      return MohrCoulombReturn(elastic_strain,f,0.0,false,pmod);
    }



  inline
    Eigen::Matrix<double,6,1> Viscoelastic(Eigen::Matrix<double,6,1> elastic_strain,
                                          double E, double nu, double viscosity,double dt) {
      Eigen::Matrix<double,3,3> dev = Eigen::Matrix<double,3,3>::Identity() - Eigen::Matrix<double,3,3>::Constant(1.0/3.0);
      Eigen::Matrix<double,3,3> De3 =
        (E/((1+nu) * (1-(2*nu))))*
        (((1-(2*nu))*Eigen::Matrix<double,3,3>::Identity()) +
        Eigen::Matrix<double,3,3>::Constant(nu));

      double G = (E/(2*(1+nu)));
      double K = (E/(3*(1-(2*nu))));

      Eigen::Matrix<double,3,3> d_neq = De3;//2*G*dev;
      Eigen::Matrix<double,3,3> C = d_neq.inverse();
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(voigt_to_matrix(elastic_strain));
      if (eigensolver.info() != Eigen::Success)
        {
          abort();
          /* return false; */
        }
      Eigen::Matrix<double,3,1> eigen_values = eigensolver.eigenvalues();//.reverse();
      Eigen::Matrix<double,3,3> eigen_vectors = eigensolver.eigenvectors();//.rowwise().reverse();

      Eigen::Matrix<double,3,1> EpsTr = eigen_values;
      Eigen::Matrix<double,3,1> en = EpsTr;
      Eigen::Matrix<double,3,1> beta;

      Eigen::Matrix<double,3,3> a = C * (C + (dev * (dt/viscosity))).inverse();
      // //std::cout<<"De3\n"<<De3<<"\n";
      // std::cout<<a<<"\n";
      const double ftol = 1e-5;
      double f = ftol;
      const int maxsteps = 1000;
      for (int i = 0;(i < maxsteps) && (f >= ftol); ++i){
        beta = d_neq * en;
        // std::cout<<beta<<"\n";
        Eigen::Matrix<double,3,1> r = EpsTr - (en + ((dev * beta) * (dt / (1.0 * viscosity))));
        f = r.norm();
        if(f >= ftol){
          en += (a * r);
        }
      }
      if(f > ftol){
        std::cout<<"Viscoelastic didn't converge "<<f<<"\n";
        abort();
      }
      elastic_strain = matrix_to_voigt(eigen_vectors * en.asDiagonal() * eigen_vectors.transpose());
      return elastic_strain;
    }


  inline
    Eigen::Matrix<double,6,1> log_strain_update(const Eigen::Matrix<double,6,1> & strain, const Eigen::Matrix3d& df){
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(voigt_to_matrix(strain));
      if (eigensolver.info() != Eigen::Success)
        {
          std::cout<<"Eigensolve failed\n";
          abort();
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
          std::cout<<"Eigensolve failed\n";
          abort();
          // return false;
        }
      auto l = trialeigensolver.eigenvalues();
      auto v = trialeigensolver.eigenvectors();
      if ((l.array() <= 0.0).any())
        {
          std::cout<<"Eigensolve failed\n";
          abort();
          // return false;
        }
      return (matrix_to_voigt(v * l.array().log().matrix().asDiagonal() * v.transpose()).array() * 0.5).matrix();
    }
}
