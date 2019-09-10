#include "pose3d_vertex.hh"
namespace SLAMSolver {
void Pose3DVertex::plus(const VecX &delta) {
  //delta(0) delta(1) delta(2) is so3, lie algebra for rotation
  //delta(3) delta(4) delta(5) is translation
  //construct delta SO3
  Sophus::SO3 original_SO3_original_plus_delta(delta(0), delta(1), delta(2));
  //delta translation
  Eigen::Vector3d delta_translation(delta(3), delta(4), delta(5));
  //delta transformation
  Sophus::SE3 delta_SE3(original_SO3_original_plus_delta, delta_translation);
  //apply delta transformation
  parameters_ = delta_SE3 * parameters_;
}
}