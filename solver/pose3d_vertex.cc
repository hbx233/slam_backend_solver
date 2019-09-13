#include "pose3d_vertex.hh"
namespace SLAMSolver {
void Pose3DVertex::plus(const VecX &delta) {
  //delta(0) delta(1) delta(2) is so3, lie algebra for rotation
  //delta(3) delta(4) delta(5) is translation
  //construct delta SO3
  Sophus::SO3 delta_SO3(delta(0), delta(1), delta(2));
  //delta translation
  Eigen::Vector3d delta_translation(delta(3), delta(4), delta(5));
#if 0
  //delta transformation
  Sophus::SE3 delta_SE3(original_SO3_original_plus_delta, delta_translation);
  //apply delta transformation
  parameters_ = delta_SE3 * parameters_;
#endif
  //Left multiplication model, all edge connected to this Pose3DVertex also need to
  //use left multiplication model
  Sophus::SO3 new_SO3 = delta_SO3 * parameters_.so3();
  Eigen::Vector3d new_translation = delta_translation + parameters_.translation();
  parameters_ = Sophus::SE3(new_SO3, new_translation);
}
}