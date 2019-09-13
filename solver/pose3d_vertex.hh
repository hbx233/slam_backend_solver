#pragma once
#include "solver/base_vertex.hh"
#include "sophus/se3.hpp"

namespace SLAMSolver {
/*!
 * @brief Graph optimization vertex representing the 3D pose
 */
class Pose3DVertex : public BaseVertex<6,Sophus::SE3>{
public:
  /// @brief Default constructor
  Pose3DVertex() = default;
  /// @brief Default destructor
  ~Pose3DVertex() = default;
  /// @brief Override the update function for 3D pose
  /// @note Use left multiplication of perturbation and separately update rotation and translation
  void plus(const VecX &delta) override;
};
}
