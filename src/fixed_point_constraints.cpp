#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {


    int l = indices.size();

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripleList;
    tripleList.reserve(q_size*q_size);

    int j=0;

    for(int i=0; i<q_size; i+=3) {
        
        if(j < indices.size() && indices[j] == i/3) 
        {
            j += 1;
            //triplets.push_back(Eigen::Triplet<double>(i-3*j, i, 1.));
        } 
        else 
        {
            for(int k=0; k<3; ++k) {
                tripleList.push_back(T(i-3*j +k ,i+k, 1.));
                //tripleList.push_back(T(i +k ,i+k, 1.));
            }
        }
        
    }

    P.resize(q_size - (3*l), q_size);
    P.setFromTriplets(tripleList.begin(), tripleList.end());
    


}