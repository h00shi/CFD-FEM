#include "Mesh/DummyMesh.h"
//****************************************************************************80
std::string DummyMesh::SetupCubicHexMesh(const realT& Lx, const realT& Ly, 
                                         const realT& Lz, const intT& Nx, 
                                         const intT& Ny, const intT& Nz)
{
  //---> Setup mesh stream
  std::ostringstream mesh_stream;
  mesh_stream.precision(std::numeric_limits<realT>::digits10);
    
  intT nnode = (Nx)*(Ny)*(Nz);
  intT ndim = 3;
  intT nelement = (Nx - 1)*(Ny - 1)*(Nz - 1);
  intT nnz_element2node = nelement*8;
  intT nbc_face = (Nx - 1)*(Ny - 1)*2 + (Nx - 1)*(Nz - 1)*2 + 
    (Ny - 1)*(Nz - 1)*2;
  intT nnz_bc_face2node = nbc_face*4;
    
  //---> Write out nodal header 
  mesh_stream << nnode << " " << ndim << std::endl;
  
  //---> Set dx, dy, dz
  realT dx = Lx/(double(Nx - 1));
  realT dy = Ly/(double(Ny - 1));
  realT dz = Lz/(double(Nz - 1));
  //---> Write coordinates as contiguous stream
  for(intT k = 0; k < Nz; k++){// k loop 
    for(intT j = 0; j < Ny; j++){// j loop 
      for(intT i = 0; i < Nx; i++){// i loop 
        realT x = i*dx;
        realT y = j*dy;
        realT z = k*dz;
        mesh_stream << x << " ";
        mesh_stream << y << " ";
        mesh_stream << z << " ";
      } // End i loop 
    } // End j loop 
  }// End k loop 
  mesh_stream << std::endl;
  
  //---> Write number of nodes per element as continuous stream
  mesh_stream << nelement << " " << nnz_element2node << std::endl;
  for(intT e = 0; e < nelement; e++){
    mesh_stream << 8 << " ";
  }
  mesh_stream << std::endl;
  //---> Write element2node;
  for(intT k = 0; k < Nz-1; k++){// k loop 
    for(intT j = 0; j < Ny-1; j++){// j loop 
      for(intT i = 0; i < Nx-1; i++){// i loop 
        intT node0 = i + j*Nx + k*Nx*Ny;
        intT node1 = node0 + 1;
        intT node2 = node0 + 1 + Nx;
        intT node3 = node0 + Nx;

        intT node4 = node0 + Nx*Ny;
        intT node5 = node0 + 1 + Nx*Ny;
        intT node6 = node0 + 1 + Nx + Nx*Ny;
        intT node7 = node0 + Nx + Nx*Ny;
        mesh_stream << node0 << " " << node1 << " " << node2 << " " << node3 
                    << " " << node4 << " " << node5 << " " << node6 << " " 
                    << " " << node7 << " ";
      } // End i loop 
    } // End j loop 
  }// End k loop  
  //---> Write boundary face header
  mesh_stream << nbc_face << " " << nnz_bc_face2node << std::endl;
  for(intT f = 0; f < nbc_face; f++){
    mesh_stream << 4 << " ";
  }
  
  intT num_faces[6];
  //---> Boundary Face2 node
  //--->Bottom face of cube
  num_faces[0] = (Ny-1)*(Nx-1);
  for(intT j = 0; j < Ny-1; j++){
    for(intT i = 0; i < Nx-1; i++){
      intT node0 = i + j*Nx; //first node of face
      mesh_stream << node0          << " ";
      mesh_stream << node0 + Nx     << " ";
      mesh_stream << node0 + 1 + Nx << " ";
      mesh_stream << node0 + 1      << " ";
    }
  }

  //z = 1 face
  num_faces[1] = (Ny-1)*(Nx-1);
  for(intT j = 0; j < Ny-1; j++){ 
    for(intT i = 0; i < Nx-1; i++){
        intT node0 = i + j*Nx+ (Nz -1)*Nx*Ny; //first node of face
      mesh_stream << node0          << " ";
      mesh_stream << node0 + 1      << " ";
      mesh_stream << node0 + 1 + Nx << " ";
      mesh_stream << node0 + Nx     << " ";
    } 
  }

  //x = 0 face of cube
  num_faces[2] = (Nz-1)*(Ny-1);
  for(intT k = 0; k < Nz-1; k++){
    for(intT j = 0; j < Ny-1; j++){
      intT node0 = j*Nx + k*Nx*Ny; //first node of face
      mesh_stream << node0              << " ";
      mesh_stream << node0 + Nx*Ny      << " ";
      mesh_stream << node0 + Nx + Nx*Ny << " ";
      mesh_stream << node0 + Nx         << " ";
    }
  }

  //x = 1 face of cube
  num_faces[3] = (Nz-1)*(Ny-1);
  for(intT k = 0; k < Nz-1; k++){
    for(intT j = 0; j < Ny-1; j++){
      intT node0 = Nx + j*Nx + k*Nx*Ny; //first node of face
      mesh_stream << node0              << " ";
      mesh_stream << node0 + Nx         << " ";
      mesh_stream << node0 + Nx + Nx*Ny << " ";
      mesh_stream << node0 + Nx*Ny      << " ";
    }
  }

  //y = 0 face of cube
  num_faces[4] = (Nz-1)*(Nx-1);
  for(intT k = 0; k < Nz-1; k++){
    for(intT i = 0; i < Nx-1; i++){
      intT node0 = i + (k)*(Nx+1)*(Ny+1); //first node of face
      mesh_stream << node0             << " ";
      mesh_stream << node0 + 1         << " ";
      mesh_stream << node0 + 1 + Nx*Ny << " ";
      mesh_stream << node0 + Nx*Ny     << " ";
    }
  }

  //---> Back face of cube
  num_faces[5] = (Nz-1)*(Nx-1);
  for(intT k = 0; k < Nz-1; k++){
    for(intT i = 0; i < Nx-1; i++){
      intT node0 = i + (Ny-1)*Nx + k*Nx*Ny; //first node of face
      mesh_stream << node0             << " ";
      mesh_stream << node0 + Nx*Ny     << " ";
      mesh_stream << node0 + 1 + Nx*Ny << " ";
      mesh_stream << node0 + 1         << " ";
    }
  }
  mesh_stream << std::endl;

  //---> Write boundary ID
  mesh_stream << nbc_face << "\n";
  for(intT bc = 0; bc < 6; bc++){ // Boundary ID
    for(intT i = 0; i < num_faces[bc]; i++){ // Faces on each bc-id
      mesh_stream << bc << " ";
    } // End Faces on each bc-id
  }// End Boundary ID
  
  mesh_stream << std::endl;
  return mesh_stream.str();
}// End DummyMesh::SetupDummyHexMesh
