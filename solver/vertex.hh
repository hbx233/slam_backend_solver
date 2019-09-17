#pragma once 

#include "solver/eigen_types.h"

namespace SLAMSolver{
/*!
 * Base interface class for all vertex
 */
class Vertex {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  /*!
   * @brief Constructor of Vertex
   * @param minimal_dimension Dimension of minimal parametrization
   */
  Vertex(const int minimal_dimension);

  /// @brief Virtual destructor
  virtual ~Vertex();

  /// @return Dimension of minimal parametrization
  int minimal_dimension() const{return minimal_dimension_;}

  /// @return ID of the vertex
  unsigned long id() const { return id_; }

  /*!
   * @brief Pure virtual function of updating the parameter(over parametrization)
   * @param delta update value need to add to parameter
   */
  virtual void plus(const Eigen::VectorXd &delta) = 0;

  /*!
   * @brief Set the Vertex to be fixed or changable
   * @param fixed True if want to set this vertex fixed
   */
  void set_fixed(bool fixed = true) {
      fixed_ = fixed;
  }

  /// @return True if the Vertex is fixed
  bool is_fixed() const { return fixed_; }

protected:
  const int minimal_dimension_; //Dimension of minimal parametrization
  unsigned long id_; //vertex id

  bool fixed_ = false;   //True if the Vertex is fixed
};
}