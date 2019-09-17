#pragma once

#include "solver/base_vertex.hh"

namespace SLAMSolver{
/*!
 * @brief Graph optimization vertex representing 3D point position
 */
class Point3DVertex : public BaseVertex<3,Eigen::Vector3d>{
public:
  /// @brief Default Constructor
  Point3DVertex() = default;
  /// @brief Default Destructor
  ~Point3DVertex() = default;
  /// @brief Trivial add
  void plus(const Eigen::VectorXd &delta) override{
    parameters_(0) += delta(0);
    parameters_(1) += delta(1);
    parameters_(2) += delta(2);
  }
};
}
