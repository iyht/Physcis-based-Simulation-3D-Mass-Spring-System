#include <V_spring_particle_particle.h>
#include <iostream>
//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    Eigen::MatrixXd B(3, 6);
    Eigen::VectorXd q(6);
    //B << -Eigen::MatrixXd::Identity(q0.rows(), q0.rows()), Eigen::MatrixXd::Identity(q0.rows(), q0.rows());
    B << -1., 0., 0., 1., 0., 0.,
          0., -1., 0., 0., 1., 0.,
          0., 0., -1., 0., 0., 1.;
    //q << q0, q1;
    q << q0(0),q0(1), q0(2), q1(0), q1(1), q1(2);
    V = pow(sqrt(q.transpose()*B.transpose()*B*q) - l0, 2)*stiffness/2.0;
    //V = pow((q0 - q1).norm() - l0, 2)*stiffness/2.0;
    //std::cout << V << std::endl;
}