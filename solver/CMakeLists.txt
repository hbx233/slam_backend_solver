add_library(SLAMSolver SHARED
        edge.hh
        edge.cc
        base_edge.hh
        vertex.hh
        base_vertex.hh
        pose3d_vertex.hh
        pose3d_vertex.cc
        problem.hh
        problem.cc
        relative_pose3d_edge.cc
        relative_pose3d_edge.hh
        base_solver.hh
        plain_solver.cc
        plain_solver.hh
        base_solver.cc
        minimizer_levenberg_marquardt.cc
        minimizer_levenberg_marquardt.hh
        sparse_schur_solver.cc
        sparse_schur_solver.hh
        sparse_cholesky_solver.cc
        sparse_cholesky_solver.hh
        vertex.cc
        reprojection_edge.hh
        reprojection_edge.cc
        )