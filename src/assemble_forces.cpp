#include <assemble_forces.h>
#include "dV_spring_particle_particle_dq.h"
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
        

    int r = q.rows();
    f.resize(r);
    //f.setZero();

    for(int i = 0; i < E.rows(); i++)
    {
        Eigen::Vector3i index_q0; // store the index to the vector q to reach [q0_x, q0_y, q0_z]
        Eigen::Vector3i index_q1; 
        index_q0 << E(i, 0)*3, E(i, 0)*3+1, E(i, 0)*3+2;
        index_q1 << E(i, 1)*3, E(i, 1)*3+1, E(i, 1)*3+2;

        Eigen::Vector3d q0, q1;
        q0 <<  q(index_q0(0)), q(index_q0(1)), q(index_q0(2));
        q1 <<  q(index_q1(0)), q(index_q1(1)), q(index_q1(2));

        Eigen::VectorXd f_01;
        f_01.resize(6);
        dV_spring_particle_particle_dq(f_01, q0, q1, l0(i), k);

        //std::cout << "f_01" << f_01 << std::endl;
        f(index_q0(0)) -= f_01(0);
        f(index_q0(1)) -= f_01(1);
        f(index_q0(2)) -= f_01(2);

        f(index_q1(0)) -= f_01(3);
        f(index_q1(1)) -= f_01(4);
        f(index_q1(2)) -= f_01(5);
    }
        
}