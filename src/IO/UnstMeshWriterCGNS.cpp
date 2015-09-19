#include "UnstMeshWriterCGNS.h"
//****************************************************************************80
UnstMeshWriterCGNS::UnstMeshWriterCGNS(const UnstMesh& mesh_ref) :
  UnstMeshWriter::UnstMeshWriter(mesh_ref)
{
   
  x_.initialize(UnstMeshWriter::mesh_.get_nnode());
  y_.initialize(UnstMeshWriter::mesh_.get_nnode());
  z_.initialize(UnstMeshWriter::mesh_.get_nnode());

  bar_con_.initialize( std::max(1,UnstMeshWriter::mesh_.get_nbar()),  2);
  tri_con_.initialize( std::max(1,UnstMeshWriter::mesh_.get_ntri()),  3);
  quad_con_.initialize( std::max(1,UnstMeshWriter::mesh_.get_nquad()),
                        4);
  tet_con_.initialize(  std::max(1,UnstMeshWriter::mesh_.get_ntet()), 4);
  prism_con_.initialize( std::max(1,UnstMeshWriter::mesh_.get_nprism()),
                         6);
  pyr_con_.initialize( std::max(1,UnstMeshWriter::mesh_.get_npyr()), 5);
  hex_con_.initialize( std::max(1,UnstMeshWriter::mesh_.get_nhex()), 8);

  bc_node_con_.initialize(
			  std::max(1, UnstMeshWriter::mesh_.get_nbc_node()), 1);
  bc_edge_con_.initialize(
			  std::max(1, UnstMeshWriter::mesh_.get_nbc_edge()), 2);
  bc_tri_con_.initialize(
			 std::max(1, UnstMeshWriter::mesh_.get_nbc_tri()),  3);
  bc_quad_con_.initialize(
			  std::max(1, UnstMeshWriter::mesh_.get_nbc_quad()), 4);
    
  bc_id_edge_.initialize(UnstMeshWriter::mesh_.get_nbc_id(),
			 UnstMeshWriter::mesh_.get_nbc_edge() );
  bc_id_tri_.initialize(UnstMeshWriter::mesh_.get_nbc_id(),
			UnstMeshWriter::mesh_.get_nbc_tri() );
  bc_id_quad_.initialize(UnstMeshWriter::mesh_.get_nbc_id(),
			 UnstMeshWriter::mesh_.get_nbc_quad() );

  //---> Using the mesh pointer setup the arrays for writing CGNS files
  //---> Loop over mesh and re-order coordinate values as separate 1-D Arrays.
  //     Unfortunately it appears that CGNS developers are equally STUPID and
  //     MORNIC to all other CS people, who don't seem to understand that I
  //     DON'T STORE MY DATA IN THEIR WAY.
  for (intT n = 0; n < UnstMeshWriter::mesh_.get_nnode(); n++) {
    x_(n) = UnstMeshWriter::mesh_.get_x()(n,0);

    switch (UnstMeshWriter::mesh_.get_ndim()) {
    case 1:
      y_(n) = 0.0;
      z_(n) = 0.0;
      break;
    case 2:
      y_(n) = UnstMeshWriter::mesh_.get_x()(n,1);
      z_(n) = 0.0;
      break;
    case 3:
      y_(n) = UnstMeshWriter::mesh_.get_x()(n,1);
      z_(n) = UnstMeshWriter::mesh_.get_x()(n,2);
      break;
    }// End switch (ndim)
  }

  intT itri = 0;
  intT iquad = 0;
  intT itet = 0;
  intT iprism = 0;
  intT ipyr = 0;
  intT ihex = 0;

  for (intT e = 0; e < UnstMeshWriter::mesh_.get_nelement(); e++) {
    // Element type connecitivty

    switch ( UnstMeshWriter::mesh_.get_element_type()(e) ) {
    case 0:
      bar_con_(e,0) =
	UnstMeshWriter::mesh_.get_element2node()(e,0) + 1;
      bar_con_(e,1) =
	UnstMeshWriter::mesh_.get_element2node()(e,1) + 1;
      break;
    case 1:
      tri_con_(itri,0) =
	UnstMeshWriter::mesh_.get_element2node()(e,0) + 1;
      tri_con_(itri,1) =
	UnstMeshWriter::mesh_.get_element2node()(e,1) + 1;
      tri_con_(itri,2) =
	UnstMeshWriter::mesh_.get_element2node()(e,2) + 1;
      itri += 1;
      break;
    case 2:
      quad_con_(iquad,0) =
	UnstMeshWriter::mesh_.get_element2node()(e,0) + 1;
      quad_con_(iquad,1) =
	UnstMeshWriter::mesh_.get_element2node()(e,1) + 1;
      quad_con_(iquad,2) =
	UnstMeshWriter::mesh_.get_element2node()(e,2) + 1;
      quad_con_(iquad,3) =
	UnstMeshWriter::mesh_.get_element2node()(e,3) + 1;
      iquad += 1;
      break;
    case 3:
      tet_con_(itet,0) =
	UnstMeshWriter::mesh_.get_element2node()(e,0) + 1;
      tet_con_(itet,1) =
	UnstMeshWriter::mesh_.get_element2node()(e,1) + 1;
      tet_con_(itet,2) =
	UnstMeshWriter::mesh_.get_element2node()(e,2) + 1;
      tet_con_(itet,3) =
	UnstMeshWriter::mesh_.get_element2node()(e,3) + 1;
      itet += 1;
      break;
    case 4:
      prism_con_(iprism,0) =
	UnstMeshWriter::mesh_.get_element2node()(e,0) + 1;
      prism_con_(iprism,1) =
	UnstMeshWriter::mesh_.get_element2node()(e,1) + 1;
      prism_con_(iprism,2) =
	UnstMeshWriter::mesh_.get_element2node()(e,2) + 1;
      prism_con_(iprism,3) =
	UnstMeshWriter::mesh_.get_element2node()(e,4) + 1;
      prism_con_(iprism,4) =
	UnstMeshWriter::mesh_.get_element2node()(e,5) + 1;
      prism_con_(iprism,5) =
	UnstMeshWriter::mesh_.get_element2node()(e,6) + 1;
      iprism += 1;
      break;
    case 5:
      pyr_con_(ipyr,0) =
	UnstMeshWriter::mesh_.get_element2node()(e,0) + 1;
      pyr_con_(ipyr,1) =
	UnstMeshWriter::mesh_.get_element2node()(e,1) + 1;
      pyr_con_(ipyr,2) =
	UnstMeshWriter::mesh_.get_element2node()(e,2) + 1;
      pyr_con_(ipyr,3) =
	UnstMeshWriter::mesh_.get_element2node()(e,3) + 1;
      pyr_con_(ipyr,4) =
	UnstMeshWriter::mesh_.get_element2node()(e,4) + 1;
      ipyr += 1;
      break;
    case 6:
      hex_con_(ihex,0) =
	UnstMeshWriter::mesh_.get_element2node()(e,1) + 1;
      hex_con_(ihex,1) =
	UnstMeshWriter::mesh_.get_element2node()(e,2) + 1;
      hex_con_(ihex,2) =
	UnstMeshWriter::mesh_.get_element2node()(e,3) + 1;
      hex_con_(ihex,3) =
	UnstMeshWriter::mesh_.get_element2node()(e,0) + 1;
      hex_con_(ihex,4) =
	UnstMeshWriter::mesh_.get_element2node()(e,5) + 1;
      hex_con_(ihex,5) =
	UnstMeshWriter::mesh_.get_element2node()(e,6) + 1;
      hex_con_(ihex,6) =
	UnstMeshWriter::mesh_.get_element2node()(e,7) + 1;
      hex_con_(ihex,7) =
	UnstMeshWriter::mesh_.get_element2node()(e,4) + 1;
      ihex += 1;
      break;
    }
  } // End Element type connecitivty

  Array1D<int> nEdgeId(UnstMeshWriter::mesh_.get_nbc_id());
  Array1D<int> nTriFaceId(UnstMeshWriter::mesh_.get_nbc_id());
  Array1D<int> nQuadFaceId(UnstMeshWriter::mesh_.get_nbc_id());
  nEdgeId.set_value(0);
  nTriFaceId.set_value(0);
  nQuadFaceId.set_value(0);

  //---> Count number of boundary edges, tris and quad that belong to each
  //     boundary id number separately.
  for(intT f = 0; f < UnstMeshWriter::mesh_.get_nbc_face(); f++) {
    intT id = UnstMeshWriter::mesh_.get_bc_face_id()(f);

    switch
      (UnstMeshWriter::mesh_.get_bc_face2node().get_ncol(f) ) {
    case 2:
      nEdgeId(id) +=1;
      break;
    case 3:
      nTriFaceId(id) += 1;
      break;
    case 4:
      nQuadFaceId(id) +=1;
      break;
    }// End swtich
  }// loop over boundary faces;

  bc_id_edge_.set_ncol(nEdgeId);
  bc_id_tri_.set_ncol(nTriFaceId);
  bc_id_quad_.set_ncol(nQuadFaceId);

  //---> Reset counters
  nEdgeId.set_value(0);
  nTriFaceId.set_value(0);
  nQuadFaceId.set_value(0);

  //---> For each boundary ID collect the face numbers of each type of face
  //     and group these faces by boundary id number
  for(intT f = 0; f < UnstMeshWriter::mesh_.get_nbc_face(); f++) {
    intT id = UnstMeshWriter::mesh_.get_bc_face_id()(f);
    switch
      (UnstMeshWriter::mesh_.get_bc_face2node().get_ncol(f) ) {
    case 2:
      bc_id_edge_(id, nEdgeId(id)) = f;
      nEdgeId(id) += 1;
      break;
    case 3:
      bc_id_tri_(id, nTriFaceId(id)) = f;
      nTriFaceId(id) += 1;
      break;
    case 4:
      bc_id_quad_(id, nQuadFaceId(id)) = f;
      nQuadFaceId(id) += 1;
      break;
    }// End swtich
  }// loop over boundary faces;

  //--------------------------- Boundary Connectivity ------------------------
  //---> The idea here is to list all boundary face connecitivities for each
  //     type by boundary id order.  So fa
  intT iBcEdge = 0;
  intT iBcTri = 0;
  intT iBcQuad = 0;
  //---> Loop over boundary ID's
  for(intT id = 0; id < UnstMeshWriter::mesh_.get_nbc_id(); id++) {

    //---> For all Edges on this Id get the boundary face number and
    // connectivity
    for(intT i = 0; i < bc_id_edge_.get_ncol(id); i++){
      //---> Get the boundary face number
      intT f = bc_id_edge_(id,i);
      //---> Copy over connectivity to bc_edge_con_ for this id
      bc_edge_con_(iBcEdge,0) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,0) + 1;
      bc_edge_con_(iBcEdge,1) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,1) + 1;
      //---> Increment boundary Edge counter;
      iBcEdge += 1;
    }

    //---> For all tris on this ID get the boundary face number and
    // connecitivity
    for(intT i = 0; i < bc_id_tri_.get_ncol(id); i++){
      //---> Get the boundary face number
      intT f = bc_id_tri_(id,i);
      //---> Copy over connectivity to bc_tri_con_ for this id
      bc_tri_con_(iBcTri,0) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,0) + 1;
      bc_tri_con_(iBcTri,1) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,1) + 1;
      bc_tri_con_(iBcTri,2) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,2) + 1;
      //---> Increment boundary Tri counter;
      iBcTri += 1;
    }

    //---> For all quadss on this ID get the boundary face number and
    // connecitivity
    for(intT i = 0; i < bc_id_quad_.get_ncol(id); i++){
      //---> Get the boundary face number
      intT f = bc_id_quad_(id,i);
      //---> Copy over connectivity to bc_quad_con_ for this id
      bc_quad_con_(iBcQuad,0) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,0) + 1;
      bc_quad_con_(iBcQuad,1) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,1) + 1;
      bc_quad_con_(iBcQuad,2) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,2) + 1;
      bc_quad_con_(iBcQuad,3) =
	UnstMeshWriter::mesh_.get_bc_face2node()(f,3) + 1;
      //---> Increment boundary Edge counter;
      iBcQuad += 1;
    }

  }// end boundary id loop

    
} // UnstMeshWriterCGNS


//****************************************************************************80
UnstMeshWriterCGNS::~UnstMeshWriterCGNS(){}


//****************************************************************************80
void UnstMeshWriterCGNS::Write(const std::string& fbase )
{

  intT funit;
  intT BaseIndex;
  intT ZoneIndex;
  intT CoordIndex;
  intT SectionIndex;

  cgsize_t isize[3][1];

  // char basename[33];
  // cgerr_t ier =
  cg_open( (fbase + std::string(".cgns") ).c_str(), CG_MODE_WRITE,
	   &funit  );

  std::cout << "Writing CGNS Mesh File: "
	    << (fbase + std::string(".cgns")) << std::endl;

  cg_base_write(funit, "Base",
		UnstMeshWriter::mesh_.get_ndim(),
		std::max(2,UnstMeshWriter::mesh_.get_ndim()),
		&BaseIndex);
  //---> Set size storage array
  //---> First number of nodes;
  isize[0][0] = UnstMeshWriter::mesh_.get_nnode();
  isize[1][0] = UnstMeshWriter::mesh_.get_nelement();
  isize[2][0] = 0;
  cg_zone_write(funit, BaseIndex, "VolumeZone", isize[0], Unstructured,
		&ZoneIndex);
  // ---> Loop over mesh and re-order coordinate values as separate 1-D
  //     Arrays.
  switch ( UnstMeshWriter::mesh_.get_ndim() ) {
  case 1:
    cg_coord_write(funit, BaseIndex, ZoneIndex, RealDouble, "CoordinateX",
		   x_.get_ptr(0), &CoordIndex);
    cg_coord_write(funit, BaseIndex, ZoneIndex, RealDouble, "CoordinateY",
		   y_.get_ptr(0), &CoordIndex);
    break;
  case 2:
    cg_coord_write(funit, BaseIndex, ZoneIndex, RealDouble, "CoordinateX",
		   x_.get_ptr(0), &CoordIndex);
    cg_coord_write(funit, BaseIndex, ZoneIndex, RealDouble, "CoordinateY",
		   y_.get_ptr(0), &CoordIndex);
    break;
  case 3:
    cg_coord_write(funit, BaseIndex, ZoneIndex, RealDouble, "CoordinateX",
		   x_.get_ptr(0), &CoordIndex);
    cg_coord_write(funit, BaseIndex, ZoneIndex, RealDouble, "CoordinateY",
		   y_.get_ptr(0), &CoordIndex);
    cg_coord_write(funit, BaseIndex, ZoneIndex, RealDouble, "CoordinateZ",
		   z_.get_ptr(0), &CoordIndex);
    break;
  }

  intT istart = 1;
  intT iend = 0;

  intT nBar =   UnstMeshWriter::mesh_.get_nbar();
  intT nTri =   UnstMeshWriter::mesh_.get_ntri();
  intT nQuad =  UnstMeshWriter::mesh_.get_nquad();
  intT nTet =   UnstMeshWriter::mesh_.get_ntet();
  intT nPrism = UnstMeshWriter::mesh_.get_nprism();
  intT nPyr =   UnstMeshWriter::mesh_.get_npyr();
  intT nHex =   UnstMeshWriter::mesh_.get_nhex();

  intT nBcTri = UnstMeshWriter::mesh_.get_nbc_tri();
  intT nBcQuad = UnstMeshWriter::mesh_.get_nbc_quad();


  intT fbegin = 0;
  switch ( UnstMeshWriter::mesh_.get_ndim() ) {

  case 1:
    iend += nBar - 1;
    cg_section_write(funit, BaseIndex, ZoneIndex, "Elem", BAR_2, istart,
		     iend, 0,
		     bar_con_.get_ptr(0,0),
		     &SectionIndex);
    istart = iend + 1;
    break;
  case 2:

    if( nTri > 0) {

      iend += nTri;
      cg_section_write(funit, BaseIndex, ZoneIndex, "Elem", TRI_3, istart,
		       iend, 0, tri_con_.get_ptr(0,0), &SectionIndex);
      istart = iend + 1;
    }
    if( nQuad > 0) {

      iend += nQuad;
      cg_section_write(funit, BaseIndex, ZoneIndex, "Elem", QUAD_4, istart,
		       iend, 0,
		       quad_con_.get_ptr(0,0),
		       &SectionIndex);
      istart = iend + 1;
    }

    fbegin = 0;
    for (intT id = 0; id < UnstMeshWriter::mesh_.get_nbc_id(); id++) {

      std::stringstream ss;
      ss << id;

      std::string BCID( "BC-" + ss.str() );
      std::cout << id << " " << BCID << std::endl;
      ss << id;
      iend += bc_id_edge_.get_ncol(id);
      cg_section_write(funit, BaseIndex, ZoneIndex, BCID.c_str(), BAR_2,
		       istart, iend, 0, bc_edge_con_.get_ptr(fbegin,0),
		       &SectionIndex);

      fbegin += bc_id_edge_.get_ncol(id);
      istart = iend + 1;
    }

    break;
  case 3:
    if( nTet > 0) {
      iend += nTet;
      cg_section_write(funit, BaseIndex, ZoneIndex, "Elem", TETRA_4, istart,
		       iend, 0,
		       tet_con_.get_ptr(0,0),
		       &SectionIndex);
      istart = iend + 1;
    }
    if(nPrism > 0 ) {
      iend += nPrism;
      cg_section_write(funit, BaseIndex, ZoneIndex, "Elem", PENTA_6, istart,
		       iend, 0,
		       prism_con_.get_ptr(0,0),
		       &SectionIndex);
      istart = iend + 1;
    }
    if( nPyr > 0 ) {
      iend += nPyr;
      cg_section_write(funit, BaseIndex, ZoneIndex, "Elem", PYRA_5, istart,
		       iend, 0,
		       pyr_con_.get_ptr(0,0),
		       &SectionIndex);
      istart = iend + 1;
    }
    if( nHex > 0 ) {
      iend += nHex;

      cg_section_write(funit, BaseIndex, ZoneIndex, "Elem", HEXA_8, istart,
		       iend, 0,
		       hex_con_.get_ptr(0,0),
		       &SectionIndex);
      istart = iend + 1;
    }

    if( nBcTri > 0 ) {
      fbegin = 0;
      for (intT id = 0; id < UnstMeshWriter::mesh_.get_nbc_id(); id++) {

	std::stringstream ss;
	ss << id;

	std::string BCID( "BC-" + ss.str() );
	std::cout << id << " " << BCID << std::endl;
	ss << id;
	iend += bc_id_tri_.get_ncol(id);
	cg_section_write(funit, BaseIndex, ZoneIndex, BCID.c_str(), TRI_3,
			 istart, iend, 0, bc_tri_con_.get_ptr(fbegin,0),
			 &SectionIndex);

	fbegin += bc_id_tri_.get_ncol(id);
	istart = iend + 1;
      }
    }

    if( nBcQuad > 0 ) {
      fbegin = 0;
      for (intT id = 0; id < UnstMeshWriter::mesh_.get_nbc_id(); id++) {

	std::stringstream ss;
	ss << id;

	std::string BCID( "BC-" + ss.str() );
	// std::cout << id << " " << BCID << std::endl;
	ss << id;
	iend += bc_id_quad_.get_ncol(id);
	cg_section_write(funit, BaseIndex, ZoneIndex, BCID.c_str(), QUAD_4,
			 istart, iend, 0, bc_quad_con_.get_ptr(fbegin,0),
			 &SectionIndex);

	fbegin += bc_id_quad_.get_ncol(id);
	istart = iend + 1;
      }
    }

    break;
  }


  cg_close(funit);

  std::cout << "Successfully wrote Unstructured Mesh to "
	    << (fbase + std::string(".cgns")) << std::endl;
}
