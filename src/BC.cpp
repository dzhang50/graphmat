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
#include <climits>
#include <ostream>

#include "/home/dzhang/zsim/misc/hooks/zsim_hooks.h"

typedef int color_type;

color_type UNKNOWN = std::numeric_limits<color_type>::max();

class BCD2 {
  public: 
    color_type color;
  public:
    BCD2() {
      color = UNKNOWN;
    }
    bool operator != (const BCD2& p) {
      return (this->color != p.color);
    }
    void print() const {
        if(color == UNKNOWN) {
            std::cout << "UNKNOWN" << std::endl;
        }
        else {
            std::cout << color << std::endl;
        }
    }

  friend std::ostream &operator<<(std::ostream &outstream, const BCD2 & val)
    {
      outstream << val.color; 
      return outstream;
    }
};

class BC2 : public GraphMat::GraphProgram<unsigned long long int, unsigned long long int, BCD2> {

  public:

  BC2() {
    this->order = GraphMat::OUT_EDGES;
    this->process_message_requires_vertexprop = false;
  }

  void reduce_function(unsigned long long int& a, const unsigned long long int& b) const {
      if(a != b) {
          std::cout << "Not bipartite graph!!" << std::endl;
          exit(0);
      }
      a = b;
  }

  void process_message(const unsigned long long int& message, const int edge_val, const BCD2& vertexprop, unsigned long long int &res) const {
    res = message;
  }

  bool send_message(const BCD2& vertexprop, unsigned long long int& message) const {
      if(vertexprop.color == 0)
          message = 1;
      else if(vertexprop.color == 1)
          message = 0;
    return true;
  }

  void apply(const unsigned long long int& message_out, BCD2& vertexprop)  {
    if (vertexprop.color == UNKNOWN) {
      vertexprop.color = message_out;
    }
    else if(vertexprop.color != message_out) {
        std::cout << "Not bipartite graph!" << std::endl;
        exit(0);
    }
  }
};

void reachable_or_not(BCD2* v, int *result, void* params=nullptr) {
  int reachable = 0;
  if (v->color != UNKNOWN) {
    reachable = 1;
  } 
  *result = reachable;
}


void run_bc(char* filename, int v) {
  GraphMat::Graph<BCD2> G;
  GraphMat::edgelist_t<int> E;
  // Loads with <nvertices nvertices nedges> header
  GraphMat::load_edgelist(filename, &E, false, true, true);
  G.ReadEdgelist(E);
  E.clear();
  
  BC2 b;

  auto b_tmp = GraphMat::graph_program_init(b, G);

  G.setAllInactive();

  auto source = G.getVertexproperty(v);
  source.color = 0;
  G.setVertexproperty(v, source);
  G.setActive(v);

  struct timeval start, end;
  gettimeofday(&start, 0);

  zsim_roi_begin();
  GraphMat::run_graph_program(&b, G, GraphMat::UNTIL_CONVERGENCE, &b_tmp);
  zsim_roi_end();

  gettimeofday(&end, 0);
  printf("Time = %.3f ms \n", (end.tv_sec-start.tv_sec)*1e3+(end.tv_usec-start.tv_usec)*1e-3);
 
  GraphMat::graph_program_clear(b_tmp);

  int reachable_vertices = 0;
  G.applyReduceAllVertices(&reachable_vertices, reachable_or_not); //default reduction = sum
  if (GraphMat::get_global_myrank() == 0) printf("Reachable vertices = %d \n", reachable_vertices);

  for(int i = 1; i < std::min((uint64_t)25, (uint64_t)G.nvertices); i++) {
      if(G.vertexNodeOwner(i)) {
          std::cout << i << ": ";
          G.getVertexproperty(i).print();
      }

  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  if (argc < 3) {
    printf("Correct format: %s A.mtx source_vertex (1-based index)\n", argv[0]);
    return 0;
  }

  int source_vertex = atoi(argv[2]);
  run_bc(argv[1], source_vertex);

  MPI_Finalize();
  
}

