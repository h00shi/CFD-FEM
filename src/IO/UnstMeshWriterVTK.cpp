#include "IO/UnstMeshWriterVTK.h"
//****************************************************************************80
UnstMeshWriterVTK::UnstMeshWriterVTK(const UnstMesh& mesh_ref) :
  UnstMeshWriter::UnstMeshWriter(mesh_ref)
{
  //---> Do some necessary initialization for the mesh
  vtk_mesh_ = vtkSmartPointer<vtkUnstructuredGrid>::New();

  intT nnode = UnstMeshWriter::mesh_.get_nnode();
  int ndim = UnstMeshWriter::mesh_.get_ndim();
  Array2D<realT> const & x = UnstMeshWriter::mesh_.get_x();
#ifdef DEV_DEBUG
  assert(ndim <= 3 and ndim > 0);
#endif  

  //---> Setup point data for insertion into vtk_mesh_
  vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
  vtk_points->SetNumberOfPoints(nnode);
  realT mesh_points[3] = {0.0, 0.0, 0.0};
  for(intT n = 0; n < nnode; n++){
    for(intT d = 0; d < ndim; d++){
      mesh_points[d] = x(n,d);
    }
    vtk_points->SetPoint(n,mesh_points);
  }
  vtk_mesh_->SetPoints(vtk_points);

  intT nelement = UnstMeshWriter::mesh_.get_nelement();
  List2D<intT> const & element2node = UnstMeshWriter::mesh_.get_element2node();
  vtkCell* vtk_cell = vtkTetra::New(); 

  for(intT e = 0; e < nelement; e++){ // Element loop 
    intT const vtk_cell_type = 
      UnstMeshWriter::
      mesh_.get_VTKType(UnstMeshWriter::mesh_.get_element_type()(e));
   
    if(vtk_cell->GetCellType() != vtk_cell_type){ //Check delete element type
      vtk_cell->Delete();
      switch(vtk_cell_type){// Switch cell type
      case VTKCellType::VTK_LINE:
	vtk_cell = vtkLine::New();
	break;
      case VTKCellType::VTK_TRIANGLE:
	vtk_cell = vtkTriangle::New();
	break;
      case VTKCellType::VTK_QUAD:
	vtk_cell = vtkQuad::New();
	break;
      case VTKCellType::VTK_TETRA:
	vtk_cell = vtkTetra::New();
	break;
      case VTKCellType::VTK_WEDGE:
	vtk_cell = vtkWedge::New();
	break;
      case VTKCellType::VTK_PYRAMID:
	vtk_cell = vtkPyramid::New();
	break;
      case VTKCellType::VTK_HEXAHEDRON:
	vtk_cell = vtkHexahedron::New();
	break;
      default:
	std::cerr << "Unsupported vtkcell_type found in WriteUnstGridVTK\n";
	exit(EXIT_FAILURE);
      }// Switch cell type
    }// Check delete current element type

    intT const nnode_per_element = vtk_cell->GetNumberOfPoints();
    for(intT j = 0; j < nnode_per_element; j++){
      vtk_cell->GetPointIds()->SetId(j,element2node(e,j));
    }
    vtk_mesh_->InsertNextCell(vtk_cell_type, vtk_cell->GetPointIds());
    
  }// End Element Loop 
  vtk_cell->Delete();
  
}// End UnstMeshWriterVTK::UnstMeshWriterVTK

//****************************************************************************80
UnstMeshWriterVTK::~UnstMeshWriterVTK()
{
}// End UnstMeshWriterVTK::~UnstMeshWriterVTK

//****************************************************************************80
void UnstMeshWriterVTK::Write(const std::string& fbase)
{
  SystemModule::cout << "Writing Mesh to VTK output to file base: " 
		     << fbase << std::endl;
  
  //---> Setup Writer;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> mesh_writer = 
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  //---> Setup Suffix;
  std::string suffix = "-" + 
    std::to_string(Communication::GetCommRank()) + ".vtu";
  
  //---> Register the mesh with the writer;
  mesh_writer->SetInputData(vtk_mesh_);
  
  //---> Register the filename with VTK object
  mesh_writer->SetFileName((fbase + suffix).c_str());
  
  //---> Tell it to write
  mesh_writer->Write();
  
  //---> Write the parallel meta data
  if (Communication::GetCommRank() == 0 and Communication::GetCommSize() > 1) {
    // Check parallel
    std::string pvtfilename = fbase + ".pvtu";
    std::ofstream outfile;
    outfile.open(pvtfilename.c_str());
    // output header information
    outfile << "<?xml version=\"1.0\"?>" << "\n";
    outfile << "<VTKFile type=\"PUnstructuredGrid\">" << "\n";
    outfile << "<PUnstructuredGrid GhostLevel=\"0\">" << "\n";

    // output point data fields
    outfile << "<PPointData>" << "\n";
    vtkSmartPointer<vtkFieldData> field_data = vtk_mesh_->GetPointData();
    for (int i = 0; field_data->GetArray(i) != NULL; i++) {
      vtkSmartPointer<vtkDataArray> array = field_data->GetArray(i);
      
      outfile << "<PDataArray type=\"Float64\" Name=\""
	      << array->GetName() << "\" NumberOfComponents=\""
	      << array->GetNumberOfComponents() <<"\"/>" << "\n";
    }
    outfile << "</PPointData>" << "\n";
    
    outfile << "<PCellData>" << "\n";
    field_data =vtk_mesh_->GetCellData();
    for (int i = 0; field_data->GetArray(i) != NULL; i++) {
      vtkSmartPointer<vtkDataArray> array = field_data->GetArray(i);
      
      outfile << "<PDataArray type=\"Float64\" Name=\""
	      << array->GetName() << "\" NumberOfComponents=\""
	      << array->GetNumberOfComponents() <<"\"/>" << "\n";
    }
    outfile << "</PCellData>" << "\n";
    
    // output points (coordinates)
    outfile << "<PPoints>" << "\n";
    outfile << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" 
	    << "\n";
    outfile << "</PPoints>" << "\n";

    // write out piece files
    std::string piece_filename;
    for (intT i = 0; i < Communication::GetCommSize(); i++) {
      // setup filename
      piece_filename = fbase + "-" + std::to_string(i) + ".vtu";
      outfile << "<Piece Source=\"" + piece_filename + "\"/>" << "\n";
    }

    // output footer
    outfile << "</PUnstructuredGrid>" << "\n";
    outfile << "</VTKFile>" << std::endl;
    outfile.close();
   
  } // End Check Parallel

  SystemModule::cout << "Done writing mesh to VTK Output" << std::endl 
		     << std::endl;
  return;
}// End UnstMeshWriterVTK::Write
