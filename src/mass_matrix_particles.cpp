#include <mass_matrix_particles.h>


void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {

    int r = q.rows();
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripleList;
    tripleList.reserve(r);
    for(int i = 0; i < r; i++)
    {
        tripleList.push_back(T(i, i, mass));
    }

    M.resize(r, r);
    M.setFromTriplets(tripleList.begin(), tripleList.end());


    
}
