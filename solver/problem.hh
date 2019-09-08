#pragma once

#include <unordered_map>
#include <map>
#include <memory>

#include "solver/eigen_types.h"
#include "solver/edge.hh"
#include "solver/vertex.hh"

namespace SLAMSolver{
class Problem {
public:
    using IDType = unsigned long;
    using HashVertexIdToEdge = std::unordered_multimap<unsigned long, std::shared_ptr<Edge>>;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    Problem();

    ~Problem();

    bool add_vertex(std::shared_ptr<Vertex> vertex);

    
    bool remove_vertex(std::shared_ptr<Vertex> vertex);

    bool add_edge(std::shared_ptr<Edge> edge);

    bool remove_edge(std::shared_ptr<Edge> edge);
    
    void get_outlier_edges(std::vector<std::shared_ptr<Edge>> &outlier_edges);

    bool solve(int max_iterations);

private:
    void compute_start_index();

    void build_gauss_newton();

    void solve_normal_equation();

    void update_vertices();

    void rollback_vertices();

    double compute_cost();

    double compute_gain_factor(const double cost_before_update, const double lambda);

    void compute_prior();

    double compute_initial_lambda(const Eigen::MatrixXd& Hessian, const double tau);

    void add_lambda_to_hessian(const double lambda);

    void remove_lambda_from_hessian(const double lambda);

    bool is_good_step_in_LM();

    double currentLambda_;
    double currentChi_;
    double stopThresholdLM_;    // 
    double ni_;                 //
    
    int problem_variable_size_ = 0;
    
    MatXX Hessian_;
    VecX b_;
    VecX delta_x_;

    // Hash map from Vertex's ID to Vertex pointer 
    std::map<IDType, std::shared_ptr<Vertex>> vertex_id_map_vertex_ptr_;
    
    // Hash map from Vertex's ID to Vertex's start index in Problem's 
    std::map<IDType, int> vertex_id_map_start_index_;

    /// all edges
    std::map<IDType, std::shared_ptr<Edge>> edge_id_map_edge_ptr_;
};

}
