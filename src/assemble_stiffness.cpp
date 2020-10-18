#include <assemble_stiffness.h>
#include "d2V_spring_particle_particle_dq2.h"
#include <iostream>

bool is_symmetric(Eigen::MatrixXd K, Eigen::MatrixXd B) {
    //Eigen::MatrixXd K_tr = K.transpose();
    for(int i = 0; i < K.rows(); i ++)
    {
        for(int j = 0; j < K.cols(); j ++)
        {
            if(fabs((K.coeff(i, j) - B.coeff(i, j))) < 0.001)
            {
                std::cout << "false:" <<K.coeff(i,j) <<","<< B.coeff(i,j)<< std::endl;
                return false;
            }

            else
            {
                std::cout << "true" << std::endl;
            }
            

        }
    }
}

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    int r = q.rows(); // r = 3n
    K.setZero();
    K.resize(r, r);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripleList;
    tripleList.reserve(r*r);

    for(int i = 0; i < E.rows(); i++)
    {
        //Eigen::Vector3d V_q0 = V.row(E.row(i)(0));
        //Eigen::Vector3d V_q1 = V.row(E.row(i)(1));
        int Vindex_q0 = E(i, 0); 
        int Vindex_q1 = E(i, 1); 

        Eigen::Vector3i index_q0; // store the index to the vector q to reach [q0_x, q0_y, q0_z]
        Eigen::Vector3i index_q1; 
        //index_q0 << E(i, 0)*3, E(i, 0)*3+1, E(i, 0)*3+2;
        //index_q1 << E(i, 1)*3, E(i, 1)*3+1, E(i, 1)*3+2;

        index_q0 << E.row(i)(0)*3, E.row(i)(0)*3+1, E.row(i)(0)*3+2;
        index_q1 << E.row(i)(1)*3, E.row(i)(1)*3+1, E.row(i)(1)*3+2;
        Eigen::Vector3d q0, q1;
        q0 <<  q(index_q0(0)), q(index_q0(1)), q(index_q0(2));
        q1 <<  q(index_q1(0)), q(index_q1(1)), q(index_q1(2));

        Eigen::Matrix66d H_i;
        d2V_spring_particle_particle_dq2(H_i, q0, q1, l0(i), k);
        //d2V_spring_particle_particle_dq2(H_i, V_q0, V_q1, l0(i), k);

        Eigen::Matrix3d H_iAA = H_i.block(0, 0, 3, 3);
        Eigen::Matrix3d H_iAB = H_i.block(0, 2, 3, 3);
        Eigen::Matrix3d H_iABT = H_i.block(2, 0, 3, 3);
        Eigen::Matrix3d H_iBB = H_i.block(2, 2, 3, 3);

        //bool t = is_symmetric(H_iAB, H_iABT);
        //tripleList.push_back(T(E(i, 0)+0, E(i, 0)+0, H_i(0, 0)));
        //tripleList.push_back(T(E(i, 0)+0, E(i, 0)+1, H_i(0, 1)));
        //tripleList.push_back(T(E(i, 0)+0, E(i, 0)+2, H_i(0, 2)));

        //tripleList.push_back(T(E(i, 0)+1, E(i, 0)+1, H_i(1, 1)));
        //tripleList.push_back(T(E(i, 0)+2, E(i, 0)+2, H_i(2, 2)));
        //tripleList.push_back(T(E(i, 0)+0, E(i, 0)+2, H_i(0, 2)));

        //tripleList.push_back(T(E(i, 0)+1, E(i, 0)+0, H_i(1, 0)));
        //tripleList.push_back(T(E(i, 0)+2, E(i, 0)+0, H_i(2, 0)));
        //tripleList.push_back(T(E(i, 0)+1, E(i, 0)+1, H_i(1, 1)));
        for(int ii = 0; ii < 3; ii++)
        {
            for (int jj = 0; jj < 3; jj++)
            {
                //std::cout << "q0:" << 3*E(i, 0) + ii << std::endl;
                int index = 3*E(i, 0);
                tripleList.push_back(T(3*Vindex_q0+ii, 3*Vindex_q0+jj, H_iAA(ii, jj)));
                tripleList.push_back(T(3*Vindex_q0+ii, 3*Vindex_q1+jj, H_iAB(ii, jj)));
                tripleList.push_back(T(3*Vindex_q1+ii, 3*Vindex_q0+jj, H_iABT(ii, jj)));
                tripleList.push_back(T(3*Vindex_q1+ii, 3*Vindex_q1+jj, H_iBB(ii, jj)));

                //tripleList.push_back(T(3*E(i, 0)+ii, 3*E(i, 0)+jj, -H_i(ii, jj)));
                //tripleList.push_back(T(3*E(i, 0)+ii, 3*E(i, 1)+jj, -H_i(ii, jj+2)));
                //tripleList.push_back(T(3*E(i, 1)+ii, 3*E(i, 0)+jj, -H_i(ii+2, jj)));
                //tripleList.push_back(T(3*E(i, 1)+ii, 3*E(i, 1)+jj, -H_i(ii+2, jj+2)));
            }
        }

        //tripleList.push_back(T(E(i, 0), E(i, 0), H_i.block(0, 0, 3, 3)));
        //tripleList.push_back(T(E(i, 0), E(i, 1), H_i.block(2, 0, 3, 3)));
        //tripleList.push_back(T(E(i, 1), E(i, 0), H_i.block(0, 2, 3, 3)));
        //tripleList.push_back(T(E(i, 1), E(i, 1), H_i.block(2, 2, 3, 3)));
    }
    K.setFromTriplets(tripleList.begin(), tripleList.end());



}