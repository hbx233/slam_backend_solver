#pragma once

#include <unordered_map>
#include <map>
#include <memory>

#include "solver/eigen_types.h"
#include "solver/edge.hh"
#include "solver/vertex.hh"

namespace SLAMSolver{
/*!
 * @brief Problem class containing all vertices connected by edges
 */
class Problem {
public:
    using IDType = unsigned long;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    /// @brief Constructor
    Problem();

    /// @brief Destructor
    ~Problem();

    /*!
     * @brief Add vertex to Problem
     * @param vertex Pointer to vertex
     * @return True if the vertex is added
     */
    bool add_vertex(std::shared_ptr<Vertex> vertex);

    /*!
     * @brief Remove vertex from problem
     * @param vertex Pointer to vertex that want to remove
     * @return True if removed
     */
    bool remove_vertex(std::shared_ptr<Vertex> vertex);

    /*!
     * @brief Add edge to problem
     * @param edge Pointer to edge that want to add
     * @return True if added
     */
    bool add_edge(std::shared_ptr<Edge> edge);

    /*!
     * @brief Remove edge from problem
     * @param edge Pointer to edge that want to remove
     * @return True if removed
     */
    bool remove_edge(std::shared_ptr<Edge> edge);
    
public:
    /*!
     * @brief Cost of the whole problem by summing cost from all edges
     * @return Accumulated total cost of problem 
     */
    double compute_cost();

    // Hash map from Vertex's ID to Vertex pointer 
    std::map<IDType, std::shared_ptr<Vertex>> vertex_id_map_vertex_ptr_;

    // Hash map from Edge's ID to Edge pointer
    std::map<IDType, std::shared_ptr<Edge>> edge_id_map_edge_ptr_;

    // Hash map from Vertex's ID to Edges that it connected to
    std::unordered_multimap<IDType, std::shared_ptr<Edge>>  vertex_id_map_edges_ptr_;
};

}
