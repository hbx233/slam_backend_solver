cmake_minimum_required(VERSION 3.2)
project(slam_backend_solver)

set(DEFAULT_BUILD_TYPE "Release")

set(CMAKE_CXX_FLAGS "-std=c++14 -g -Wall")


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wno-reorder" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -DNDEBUG -Wno-reorder -O2" CACHE STRING "" FORCE)


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# third party libs
# eigen
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})


# sophus
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/Sophus)

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(solver)
add_subdirectory(app)
