#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <map>

// double matricies
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatXX;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;
// Quaternions
typedef Eigen::Quaterniond Qd;
typedef Eigen::Quaternionf Qf;

