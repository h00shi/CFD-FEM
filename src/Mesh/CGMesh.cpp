/*
 * CGMesh.cpp
 *
 *  Created on: Oct 10, 2015
 *      Author: rabbit
 */
#include "Mesh/CGMesh.h"
//****************************************************************************80
CGMesh::CGMesh(UnstMeshReader& mesh_reader) : UnstMesh::UnstMesh(mesh_reader)
{
  UnstMesh::mesh_bc_faces_.FormBcFace2Element(UnstMesh::mesh_elements_);
  FormGraph();
}
//****************************************************************************80
CGMesh::~CGMesh()
{

  delete graph_;
}
void CGMesh::FormGraph()
{
  const List2D<intT>& node2element = mesh_elements_.get_node2element();
  const List2D<intT>& element2node = mesh_elements_.get_element2node();
  //---> Initialize size of list to zero
  intT nnz_adj = 0;
  //---> Loop over the nodes
  for(intT n = 0; n < UnstMesh::mesh_geometry_.get_nnode(); n++) {//Node_loop
    intT nn_tot = 0;
    //---> Loop over the elements that surround this node
    for(intT i = 0; i < node2element.get_ncol(n); i++){
      intT e = node2element(n,i);
      nn_tot += element2node.get_ncol(e);
    }
    //---> Temporary array for nodes attached to node n
    std::vector<int> temp (nn_tot,-1);
    //---> Initialize counter
    intT counter = 0;
    //---> Now loop over the elements containing that node
    for( intT i = 0; i < node2element.get_ncol(n); i++ ) {// element_loop
      //---> Get element containing node n
      intT e = node2element(n, i);
      //---> Loop over nodes attached to elem
      for ( intT j = 0; j < element2node.get_ncol(e); j++){ // Node_on_element
        //---> Get node index of jth node attached to elem
        intT node = element2node(e, j);
        //---> Store node in temp
        temp[counter] = node;
        //---> Increment counter
        counter += 1;
      } // Node_on_element
    } // End element_loop

    //---> Sort all possible nodes adjacent to node n in lexigraphical order
    std::sort (temp.begin(), temp.end());
    //---> Define variables to use unique function
    std::vector<int>::iterator itemp;
    /*---> Eleminate all duplicate nodes from list temp NOTE: Must call
       sort first for std::unique to work as expected.  */
    itemp = std::unique(temp.begin(), temp.end());
    //---> Resize temp to remove undefined entries
    temp.resize(std::distance( temp.begin(),itemp ) );
    //---> Return size of vector temp which is number of unique adjacent nodes
    intT c = temp.size();
    nnz_adj += c;

  } // End Node_loop
  List2D<intT> adj;
  //--------------------------------------------------------------------------
  adj.initialize(UnstMesh::mesh_geometry_.get_nnode(), nnz_adj);
  //--------------------------------------------------------------------------

  //---> Build Adjacency List
  //---> Loop over the nodes
  for (intT n = 0; n < UnstMesh::mesh_geometry_.get_nnode(); n++) {//Node_loop
    intT nn_tot = 0;
    //---> Loop over the elements that surround this node
    for(intT i = 0; i < node2element.get_ncol(n); i++){
      intT e = node2element(n,i);
      nn_tot += element2node.get_ncol(e);
    }
    //---> Temporary array for nodes attached to node n
    std::vector<int> temp (nn_tot,-1);
    //---> Initialize counter
    intT counter = 0;
    //---> Now loop over the elements containing that node
    for( intT i = 0; i < node2element.get_ncol(n); i++ ) {// element_loop
      //---> Get element containing node n
      intT e = node2element(n, i);
      //---> Loop over nodes attached to elem
      for ( intT j = 0; j < element2node.get_ncol(e); j++){ // Node_on_element
        //---> Get node index of jth node attached to elem
        intT node = element2node(e, j);
        //---> Store node in temp
        temp[counter] = node;
        //---> Increment counter
        counter += 1;
      } // Node_on_element
    } // End element_loop
    //---> Sort all possible nodes adjacent to node n in lexigraphical order
    std::sort (temp.begin(), temp.end());
    //---> Define variables to use unique function
    std::vector<int>::iterator itemp;
    /*---> Eleminate all duplicate nodes from list temp NOTE: Must call
       sort first for std::unique to work as expected.  */
    itemp = std::unique(temp.begin(), temp.end());
    //---> Resize temp to remove undefined entries
    temp.resize(std::distance( temp.begin(),itemp ) );
    //---> Return size of vector temp which is number of unique adjacent nodes
    intT c = temp.size();
    adj.set_ncol(n,c);
    for (intT i = 0; i < c; i++){adj(n,i) = temp[i];}
  } // End Node_loop

  graph_ = new Graph(adj);

} // End CGMesh::FormGraph();
