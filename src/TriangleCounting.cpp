/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* ******************************************************************************/
/* Narayanan Sundaram (Intel Corp.)
 * ******************************************************************************/

#include "GraphMatRuntime.h"
#include <vector>
#include <algorithm>
#include <assert.h>
#include <boost/serialization/vector.hpp>
#include <iostream>

#include "/home/dzhang/zsim/misc/hooks/zsim_hooks.h"

class TC : public GraphMat::Serializable {
  public:
    int id;
    int triangles;
    std::vector<int> neighbors;

  public:
    TC() {
      triangles = 0;
    }
    int operator!=(const TC& t) const {
      return (t.triangles != this->triangles); //dummy
    }
    void print() const {
      printf("%d : %d : ", id, triangles);
      printf("\n");
    }
    friend std::ostream& operator<<(std::ostream& out, const TC& t) {
      out << t.triangles;
      return out;
    }

    friend boost::serialization::access;
    template<class Archive> 
    void serialize(Archive& ar, const unsigned int version) {
      ar & id;
      ar & triangles;
      ar & neighbors;
    }
};

template<typename T>
class serializable_vector : public GraphMat::Serializable {
  public:
    std::vector<T> v;
  public:
    friend boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
      ar & v;
    }
};

class GetNeighbors : public GraphMat::GraphProgram<int, serializable_vector<int>, TC> {

  public:

  GetNeighbors() {
    this->activity = GraphMat::ALL_VERTICES;
    this->order = GraphMat::IN_EDGES;
    this->process_message_requires_vertexprop = false;
  }

  void reduce_function(serializable_vector<int>& a, const serializable_vector<int>& b) const {
    a.v.insert(a.v.end(), b.v.begin(), b.v.end()); 
  }

  void process_message(const int& message, const int edge_val, const TC& vertexprop, serializable_vector<int>& res) const {
    res.v.clear(); 
    res.v.push_back(message);
  }

  bool send_message(const TC& vertexprop, int& message) const {
    message = vertexprop.id;
    return true;
  }

  void apply(const serializable_vector<int>& message_out, TC& vertexprop) {
    vertexprop.neighbors = message_out.v;
    std::sort(vertexprop.neighbors.begin(), vertexprop.neighbors.end());
  }

};


class CountTriangles: public GraphMat::GraphProgram<TC, int, TC> {

  public:

  CountTriangles() {
    this->activity = GraphMat::ALL_VERTICES;
    this->order = GraphMat::OUT_EDGES;
  }

  void reduce_function(int& v, const int& w) const {
    v += w;
  }

  void process_message(const TC& message, const int edge_val, const TC& vertexprop, int& res) const {
    res = 0;
    int it1 = 0, it2 = 0;
    int it1_end = message.neighbors.size();
    int it2_end = vertexprop.neighbors.size();

    /*
    std::cout << "Vertexprop id: " << vertexprop.id-1 << ", message id: " << message.id-1 << "\n";
    std::cout << "Vertexprop neighbors: ";
    for(int i = 0; i < vertexprop.neighbors.size(); i++)
        std::cout << vertexprop.neighbors[i]-1 << " ";
    std::cout << "\n";
    std::cout << "Message neighbors: ";
    for(int i = 0; i < message.neighbors.size(); i++)
        std::cout << message.neighbors[i]-1 << " ";
    std::cout << "\n";
    */

    // Only check one of the two pairs
    if(message.id < vertexprop.id)
        return;

    // To avoid counting the same triangle 3 (or 6) times, enforce vertexprop.id  < message.id < shared_neighbor
    while (it1 != it1_end && it2 != it2_end){
      if (message.neighbors[it1] == vertexprop.neighbors[it2]) {
        if(message.id < vertexprop.neighbors[it2])
          res++;
        ++it1; ++it2;
      } else if (message.neighbors[it1] < vertexprop.neighbors[it2]) {
        ++it1;
      } else {
        ++it2;
      }
    } 
    //std::cout << "Res: " <<res<<"\n\n";
    return;
  }

  bool send_message(const TC& vertexprop, TC& message) const {
    message = vertexprop;
    return true;
  }

  void apply(const int& message_out, TC& vertexprop) {
    vertexprop.triangles += message_out;
  }


};

void return_triangles(TC* v, unsigned long int* out, void* params) {
  *out = v->triangles;
}

void run_triangle_counting(char* filename) {
  GraphMat::Graph<TC> G;
  GraphMat::edgelist_t<int> E;
  // Loads with <nvertices nvertices nedges> header
  GraphMat::load_edgelist(filename, &E, false, true, true);
  G.ReadEdgelist(E);
  E.clear();

  int numberOfVertices = G.getNumberOfVertices();
  GetNeighbors gn;
  CountTriangles ct;

  auto gn_tmp = GraphMat::graph_program_init(gn, G);
  auto ct_tmp = GraphMat::graph_program_init(ct, G);
  
  struct timeval start, end;

  for (int i = 1; i <= numberOfVertices; i++) {
    if (G.vertexNodeOwner(i)) {
      TC vp = G.getVertexproperty(i);
      vp.id = i;
      G.setVertexproperty(i, vp);
    }
  }
  gettimeofday(&start, 0);

  GraphMat::run_graph_program(&gn, G, 1, &gn_tmp);

  zsim_roi_begin();
  GraphMat::run_graph_program(&ct, G, 1, &ct_tmp);
  zsim_roi_end();

  gettimeofday(&end, 0);
  printf("Time = %.3f ms \n", (end.tv_sec-start.tv_sec)*1e3+(end.tv_usec-start.tv_usec)*1e-3);

  GraphMat::graph_program_clear(gn_tmp);
  GraphMat::graph_program_clear(ct_tmp);  

  unsigned long int ntriangles = 0;
  G.applyReduceAllVertices(&ntriangles, return_triangles, GraphMat::AddFn);
  if(GraphMat::get_global_myrank() == 0) printf("Total triangles = %lu \n", ntriangles);
  
  for (int i = 1; i <= std::min(10, numberOfVertices); i++) {
    if (G.vertexNodeOwner(i)) {
      G.getVertexproperty(i).print();
    }
  }

}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    printf("Correct format: %s A.mtx\n", argv[0]);
    return 0;
  }
  MPI_Init(&argc, &argv);
 
  run_triangle_counting(argv[1]); 
  MPI_Finalize();
}

