#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "solver/problem.hh"
#include "utils/tic_toc.h"

namespace SLAMSolver {

Problem::Problem()
{
}

Problem::~Problem() {}

bool Problem::add_vertex(std::shared_ptr<Vertex> vertex) {
  if (vertex_id_map_vertex_ptr_.find(vertex->id()) != vertex_id_map_vertex_ptr_.end()) {
    return false;
  } else {
    vertex_id_map_vertex_ptr_[vertex->id()] = vertex;
    return true;
  }
}


bool Problem::add_edge(std::shared_ptr<Edge> edge) {
  //Add Edge Interface pointer to hash map with its ID as key 
  if (edge_id_map_edge_ptr_.find(edge->id()) != edge_id_map_edge_ptr_.end()) {
    //Already been added
    return false;
  } else {
    //Has not been added to the problem
    edge_id_map_edge_ptr_[edge->id()] = edge;
    for (int i = 0; i < edge->num_vertices(); i++){
      auto vertex_interface_ptr = edge->get_vertex_interface(i);
      add_vertex(vertex_interface_ptr);
      //Associate the vertex id with edge
      vertex_id_map_edges_ptr_.insert(std::pair<IDType,std::shared_ptr<Edge>>(vertex_interface_ptr->id(),edge));
    }
    return true;
  }
}

double Problem::compute_cost() {
  //Compute total cost
  double total_cost = 0;
  for (auto &edge: edge_id_map_edge_ptr_) {
    //Compute errors and error function's jacobians
    edge.second->compute_errors();
    total_cost += edge.second->chi2();
  }
  return total_cost;
}

}






