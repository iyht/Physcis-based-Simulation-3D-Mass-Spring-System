#include <dV_spring_particle_particle_dq.h>
#include <iostream>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    Eigen::MatrixXd B(q0.rows(), q0.rows()*2);
    Eigen::VectorXd q(6);
    //B << -Eigen::MatrixXd::Identity(q0.rows(), q0.rows()), Eigen::MatrixXd::Identity(q0.rows(), q0.rows());
    B << -1., 0., 0., 1., 0., 0.,
          0., -1., 0., 0., 1., 0.,
          0., 0., -1., 0., 0., 1.;
    //q << q0, q1;
    q << q0(0),q0(1), q0(2), q1(0), q1(1), q1(2);

    double qBBq = sqrt(q.transpose()*B.transpose()*B*q);
    //f << (l0 - pow(q.transpose()*B.transpose()*B*q, -0.5))*(B.transpose()*B*q)*(stiffness*0.5);

    //f << ((stiffness*(qBBq - l0))/qBBq) * (q.transpose()*B.transpose()*B);
    //f(0) = -stiffness*(q(0)*2.0+q(3)*2.0)*(l0-sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5))))*1.0/sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5)))*(-1.0/2.0);
    //f(1) = -stiffness*(q(1)*2.0+q(4)*2.0)*(l0-sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5))))*1.0/sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5)))*(-1.0/2.0);
    //f(2) = -stiffness*(q(2)*2.0+q(5)*2.0)*(l0-sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5))))*1.0/sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5)))*(-1.0/2.0);
    //f(3) = -stiffness*(q(0)*2.0+q(3)*2.0)*(l0-sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5))))*1.0/sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5)))*(-1.0/2.0);
    //f(4) = -stiffness*(q(1)*2.0+q(4)*2.0)*(l0-sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5))))*1.0/sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5)))*(-1.0/2.0);
    //f(5) = -stiffness*(q(2)*2.0+q(5)*2.0)*(l0-sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5))))*1.0/sqrt(q(0)*(q(0)+q(3))+q(1)*(q(1)+q(4))+q(3)*(q(0)+q(3))+q(2)*(q(2)+q(5))+q(4)*(q(1)+q(4))+q(5)*(q(2)+q(5)))*(-1.0/2.0);

    //f(0) = stiffness*(q(0)*2.0-q(3)*2.0)*(l0-sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))*1.0/sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5)))*(-1.0/2.0);
    //f(1) = stiffness*(q(1)*2.0-q(4)*2.0)*(l0-sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))*1.0/sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5)))*(-1.0/2.0);
    //f(2) = stiffness*(q(2)*2.0-q(5)*2.0)*(l0-sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))*1.0/sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5)))*(-1.0/2.0);
    //f(3) = (stiffness*(q(0)*2.0-q(3)*2.0)*(l0-sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))*1.0/sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))/2.0;
    //f(4) = (stiffness*(q(1)*2.0-q(4)*2.0)*(l0-sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))*1.0/sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))/2.0;
    //f(5) = (stiffness*(q(2)*2.0-q(5)*2.0)*(l0-sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))*1.0/sqrt(q(0)*(q(0)-q(3))+q(1)*(q(1)-q(4))-q(3)*(q(0)-q(3))+q(2)*(q(2)-q(5))-q(4)*(q(1)-q(4))-q(5)*(q(2)-q(5))))/2.0;
    //f.setZero();

    //f(0) = -stiffness*(l0-sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))*(q(0)-q(3))*1.0/sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5));
    //f(1) = -stiffness*(l0-sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))*(q(1)-q(4))*1.0/sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5));
    //f(2) = -stiffness*(l0-sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))*(q(2)-q(5))*1.0/sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5));
    //f(3) = (stiffness*(q(0)*2.0-q(3)*2.0)*(l0-sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))*1.0/sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))/2.0;
    //f(4) = (stiffness*(q(1)*2.0-q(4)*2.0)*(l0-sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))*1.0/sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))/2.0;
    //f(5) = (stiffness*(q(2)*2.0-q(5)*2.0)*(l0-sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))*1.0/sqrt(q(0)*q(3)*-2.0-q(1)*q(4)*2.0-q(2)*q(5)*2.0+q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4)+q(5)*q(5)))/2.0;
    //std::cout << "f" << f << std::endl;
    //f = f;

    f(0) = -stiffness*(q0(0)-q1(0))*(l0-sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0)))*1.0/sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0));
    f(1) = -stiffness*(q0(1)-q1(1))*(l0-sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0)))*1.0/sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0));
    f(2) = -stiffness*(q0(2)-q1(2))*(l0-sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0)))*1.0/sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0));
    f(3) = stiffness*(q0(0)-q1(0))*(l0-sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0)))*1.0/sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0));
    f(4) = stiffness*(q0(1)-q1(1))*(l0-sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0)))*1.0/sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0));
    f(5) = stiffness*(q0(2)-q1(2))*(l0-sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0)))*1.0/sqrt(pow(q0(0)-q1(0),2.0)+pow(q0(1)-q1(1),2.0)+pow(q0(2)-q1(2),2.0));

    //f << 0, 0, 0, 0, 0, 0;
}