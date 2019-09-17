#pragma once
#include "solver/base_vertex.hh"
#include "sophus/se3.hpp"

namespace SLAMSolver {
/*!
 * @brief Graph optimization vertex representing the 3D pose
 */
class Pose3DVertex : public BaseVertex<6,Sophus::SE3>{
public:
  /*!
   * @brief Default constructor
   * @param combined True if combine the translation and rotation in optimization,
   */
  Pose3DVertex(const bool combined);
  /// @brief Default destructor
  ~Pose3DVertex() = default;
  /// @brief Override the update function for 3D pose
  /// @note Use left multiplication of perturbation and separately update rotation and translation
  void plus(const Eigen::VectorXd &delta) override;
private:
  /// True if want to use full se3
  /// False if want to optimize on rotation and translation separately
  bool combined_ = true;
};
}
