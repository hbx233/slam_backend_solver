#ifndef MYSLAM_BACKEND_EDGE_H
#define MYSLAM_BACKEND_EDGE_H

#include <memory>
#include <string>
#include "solver/eigen_types.h"

namespace SLAMSolver{

class Vertex;

class Edge {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Edge(const int errors_dimension, const int num_vertices);

  virtual ~Edge(){};

  unsigned long id() const { return id_; }

  int errors_dimension() const {return errors_dimension_;}
    
  int num_vertices() const { return num_vertices_; }
    
  std::shared_ptr<Vertex> get_vertex_interface(const int i);

  virtual void compute_errors() = 0;

  virtual void compute_jacobians() = 0;

  //errors_.transpose() * information_ * errors_
  double chi2();

  void set_information(const MatXX &information) {
    information_ = information;
  }

  MatXX information() const {
    return information_;
  }

  int ordering_id() const { return ordering_id_; }

  void set_ordering_id(unsigned long id) { ordering_id_ = id; };
    
  VecX errors_; // errors
  std::vector<MatXX> jacobians_;
protected:
    const int errors_dimension_; //Dimension of errors
    const int num_vertices_; //Number of vertices
    unsigned long id_;  // unique edge id
    int ordering_id_;   //ordering edge id in problem
    std::vector<std::shared_ptr<Vertex>> vertices_interface_ptr_;
    MatXX information_;
};

}

#endif
