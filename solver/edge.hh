#ifndef MYSLAM_BACKEND_EDGE_H
#define MYSLAM_BACKEND_EDGE_H

#include <memory>
#include <string>
#include "solver/eigen_types.h"
#include <iostream>
namespace SLAMSolver{

class Vertex;

/*!
 * @brief Base class for all edges
 */
class Edge {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /*!
   * @brief Constructor of Edge base class
   * @param errors_dimension Dimension of errors vector
   * @param num_vertices Number of vertices connected by this edge
   */
  Edge(const int errors_dimension, const int num_vertices);

  /// @brief Virtual destructor of Edge base class
  virtual ~Edge(){};

  /// @return ID of the Edge
  unsigned long id() const { return id_; }

  /// @return Dimension of errors vector
  int errors_dimension() const {return errors_dimension_;}

  /// @return Number of vertices connected by this edge
  int num_vertices() const { return num_vertices_; }

  /*!
   * @brief Get Vertex class pointer(not BaseVertex template)
   * @param i The index of Vertex
   * @return Shared pointer to Vertex
   */
  std::shared_ptr<Vertex> get_vertex_interface(const int i);

  /*!
   * @brief Pure virtual function for computing the errors, every Edge class need to implement this
   */
  virtual void compute_errors() = 0;

  /*!
   * @brief Pure virtual function for computing the jacobians, every Edge class need to implement this
   */
  virtual void compute_jacobians() = 0;

  /// @return cost
  //errors_.transpose() * information_ * errors_
  double chi2();

  /// @brief Set information matrix in the Edge
  void set_information(const MatXX &information) {
    information_ = information;
    std::cout<<information_<<std::endl;
  }

  /// @return Information matrix
  Eigen::MatrixXd information() const {
    return information_;
  }

  /// Errors vector and jacobians matrix
  Eigen::VectorXd errors_; // errors
  std::vector<Eigen::MatrixXd> jacobians_; //every vertex that the edge connected to has one jacobian matrix
protected:
    const int errors_dimension_; //Dimension of errors
    const int num_vertices_; //Number of vertices
    unsigned long id_;  // unique edge id
    std::vector<std::shared_ptr<Vertex>> vertices_interface_ptr_; //interface pointer of type Vertex
    Eigen::MatrixXd information_; //Information matrix
};

}

#endif
