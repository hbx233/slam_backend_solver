#include "solver/vertex.hh"
#include "solver/edge.hh"
#include <iostream>

using namespace std;


namespace SLAMSolver {

Edge::Edge(const int errors_dimension, const int num_vertices) : 
errors_dimension_(errors_dimension), num_vertices_(num_vertices)
{
  static unsigned long uuid_edge= 0;
  id_ = uuid_edge++;
  //Resize the Vertex pointer vector to number of vertices 
  vertices_interface_ptr_.resize(num_vertices_);
  //Resize Error Vector
  errors_.resize(errors_dimension_,1);
  //Resize Jacobian Matrices
  jacobians_.resize(num_vertices_);
  //Set information matrix
  Eigen::MatrixXd information(errors_dimension_, errors_dimension_);
  information.setIdentity();
  information_ = information;
}


shared_ptr< Vertex > Edge::get_vertex_interface(const int i) {
    //Check the vertex index
    if (i >= num_vertices_) {
        std::cout << "[ERROR] Access Vertex Number out of Bound" << std::endl;
        exit(0);
    }
    return vertices_interface_ptr_[i];
}


double Edge::chi2() {
    // TODO::  we should not Multiply information here, because we have computed Jacobian = sqrt_info * Jacobian
    return errors_.transpose() * information_ * errors_;
//    return residual_.transpose() * residual_;   // 当计算 residual 的时候已经乘以了 sqrt_info, 这里不要再乘
}


}
