#include <vector>
#include <utility>
#include <algorithm>
#include "DataStructures/Heap.h"
#include "my_incl.h"
#include <stdlib.h>
#include "mpi.h"

const intT N = 1000000;

class pair_comp {
public:
  bool operator()(const std::pair<intT,intT>& lhs, 
             const std::pair<intT,intT>& rhs)
  {
    return lhs.second <= rhs.second;
  }
};// End Class
int main(intT argc, char** agrv)
{

  Heap<intT> heap(N);
  std::vector< std::pair<intT,intT> > v(N);
  srand(100);
  std::pair<intT,intT> p;
  for(intT i = 0; i < N; i++){
    intT n = rand();
    heap.set_metric(i,n);
    p = std::make_pair(i,n);
    v.push_back(p);
  }
  heap.heapify(N);
  
  Array1D<intT> result(N);
  for(intT i = 0; i< N; i++){
   intT j = heap.pop_off_top(N-i);
   result(i) = heap.get_metric(j);
  }

  std::sort(v.begin(), v.end(), pair_comp());
  // std::cout << result << std::endl;
  srand(NULL);
  
  return 0; 

}
