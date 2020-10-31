#include <T_particle.h>
#include <iostream>

void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {


    Eigen::MatrixXd M(qdot.rows(), qdot.rows());
    M = Eigen::MatrixXd::Identity(qdot.rows(), qdot.rows());
    M = mass*M;
    T = (qdot.transpose()*M*qdot/2.0)(0);
}