#include <vector>
#include "DataStructures/List2D.h"
#include "my_incl.h"
#include <mpi.h>

class Triangle
{
public:
  Triangle(const intT& node1, const intT& node2, const intT& node3)
  {
    nodes_[0] = node1;
    nodes_[1] = node2;
    nodes_[2] = node3;
  }
  Triangle(const intT (&nodes)[])
  {
    for(intT i = 0; i < 3; i++){
      nodes_[i] = nodes[i];
    }
  }
  //Triangle(const std::initializer_list<intT>& nodes);
  inline const intT& GetGlobalNode(const intT& local_node)
  {
    return nodes_[local_node];
  }
private:
  intT nodes_[3];
  intT edges_[3];
  realT volume;
  intT neighbor_elements_[3];
  intT elem_global_;
};

int main(intT argc, char* argv[])
{
  intT nelement = 0;
  std::cout << "Enter Number of elements: " ;
  std::cin >> nelement;
  std::cout << nelement << std::endl;
  List2D<intT> elem2node(nelement,3*nelement);
  for(intT i = 0; i < nelement; i++){
    elem2node.set_ncol(i,3);
  }
  
  double ts = MPI::Wtime();
  for(intT i = 0; i < nelement; i++){
    elem2node(i,0) = i;
    elem2node(i,1) = i + 1;
    elem2node(i,2) = i + 2;
  }
  double te = MPI::Wtime();
  std::cout << "Me: " << te - ts << std::endl;
  
  std::vector<Triangle> e2n;
  ts = MPI::Wtime();
  for(intT i = 0; i < nelement; i++){
    e2n.emplace_back(Triangle(i,i+1,i+2));
  } 
  te = MPI::Wtime();
  std::cout << "Vector: " << te - ts << std::endl;
  
  ts = MPI::Wtime();
  Array1D<realT> a(nelement);
  for(intT i = 0; i < nelement; i++){
    a(i) = (elem2node(i,0) + 
           elem2node(i,1) + 
            elem2node(i,2))/3.0;

  }
  te = MPI::Wtime();
  std::cout << "Me: " << te - ts << std::endl;
  ts = MPI::Wtime();
  std::vector<realT> avg(nelement);
  for(intT i = 0; i < nelement; i++){
    avg.push_back((e2n[i].GetGlobalNode(0) + 
                   e2n[i].GetGlobalNode(1) + 
                   e2n[i].GetGlobalNode(2))/3.0);
  } 
  te = MPI::Wtime();
  std::cout << "Vector: " << te - ts << std::endl;
}

