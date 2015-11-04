/*
 * PetscSolver.cpp
 *
 *  Created on: Oct 25, 2015
 *      Author: rabbit
 */
#include "LinearSolver/PetscSolver.h"

PetscSolver::PetscSolver(CSRMatrix<realT>& A,
    intT num_its) :
    LinearSolver(5.0e-16, 5.0e-16), A_(A),
    num_its_(num_its)
{

  InitializeGlobalArrays();
  nrow_ = A_.get_row_offset().get_size(0) - 1;
  InitializePetscObjects();

}

PetscSolver::~PetscSolver()
{
 if(destroy_objects_) {DestroyPetscObjects();}
 return;
}

bool PetscSolver::Factorize()
{
#ifdef PETSC

  for(intT i = 0; i < nrow_; i++){
    intT rstart = A_.get_row_offset()(i);
    intT nnz_row = A_.get_row_offset()(i + 1) - rstart;
    intT row = petsc_global_row_(i);

    MatSetValues(Amat_, 1, &row, nnz_row,
        colidx_global_.get_ptr(rstart),
        A_.get_data().get_ptr(rstart),
        INSERT_VALUES);
  }

   MatAssemblyBegin(Amat_, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(Amat_, MAT_FINAL_ASSEMBLY);

   //MatView(Amat_, PETSC_VIEWER_STDOUT_WORLD);

#endif

return(true);
}

void PetscSolver::Solve(const List2D<double>& b, List2D<double>& x)
{
#ifdef PETSC


  VecSetValues(bvec_, nrow_, petsc_global_row_.get_ptr(0),
         b.get_ptr(0,0), INSERT_VALUES);

  VecAssemblyBegin(bvec_);
  VecAssemblyEnd(bvec_);

  VecSetValues(xvec_, nrow_, petsc_global_row_.get_ptr(0),
      x.get_ptr(0,0), INSERT_VALUES);

  VecAssemblyBegin(xvec_);
  VecAssemblyEnd(xvec_);

  KSPSetOperators(ksp_, Amat_,Amat_);
  KSPSolve(ksp_,bvec_,xvec_);

  VecGetValues(xvec_, nrow_, petsc_global_row_.get_ptr(0),
      x.get_ptr(0,0));

#endif

  return;
}

//****************************************************************************80
void PetscSolver::InitializeGlobalArrays()
{

  row_start_.initialize(Communication::GetCommSize() + 1);
  petsc_global_row_.initialize(A_.get_nrow_per_node());
  intT nrow = A_.get_row_offset().get_size(0) - 1;

  Communication::AllGather(&nrow, 1, row_start_.get_ptr(1), 1);

  for(intT p = 0; p < Communication::GetCommSize(); p++){
    row_start_(p+1) += row_start_(p);
  }

  for(intT i = 0; i < nrow; i++){
    petsc_global_row_(i) = i + row_start_(Communication::GetCommRank());
  }

  Array1D<intT> ones(A_.get_nrow_per_node().get_size(0));

  ones.set_value(1);
  if( Communication::GetCommRank() > 1) {
//    Communicator::mpi_communicator->SendData(petsc_global_row_.get_ptr(0),
//        A_.get_nrow_per_node());
//    Communicator::mpi_communicator->RecvData(petsc_global_row_.get_ptr(0),
//        A_.get_nrow_per_node());
//    Communicator::mpi_communicator->Waitall();
  }

  colidx_global_.initialize_copy_pattern(A_.get_column_idx());

  for(intT i = 0; i < colidx_global_.get_total_size(); i++){
    intT local_col = A_.get_column_idx()(i);
    colidx_global_(i) = petsc_global_row_(local_col);
  }

}
//****************************************************************************80
void PetscSolver::InitializePetscObjects()
{
#ifdef PETSC

  PetscInitialize(PETSC_NULL,PETSC_NULL,PETSC_NULL, PETSC_NULL);

  MatCreateMPIAIJWithArrays(MPI_COMM_WORLD, nrow_, nrow_,
      PETSC_DETERMINE, PETSC_DETERMINE, A_.get_row_offset().get_ptr(0),
      colidx_global_.get_ptr(0), A_.get_data().get_ptr(0), &Amat_);

  MatAssemblyBegin(Amat_, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Amat_, MAT_FINAL_ASSEMBLY);

  VecCreate(MPI_COMM_WORLD, &xvec_);
  VecSetSizes(xvec_, nrow_, PETSC_DETERMINE);
  VecSetType(xvec_, VECMPI);

  VecCreate(MPI_COMM_WORLD, &bvec_);
  VecSetSizes(bvec_, nrow_, PETSC_DETERMINE);
  VecSetType(bvec_, VECMPI);

  KSPCreate(PETSC_COMM_WORLD, &ksp_);
  KSPSetType(ksp_, KSPPREONLY);
  KSPGetPC(ksp_, &pc_);
  PCSetType(pc_, PCHYPRE);
  PCHYPRESetType(pc_, "boomeramg");
  PetscOptionsSetValue("-pc_hypre_bommeramg_cycle_type","V");
  PetscOptionsSetValue("-pc_hypre_boomeramg_cycle_max_levels","25");
  std::string its = std::to_string(num_its_);
  PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", its.c_str());
  KSPSetFromOptions(ksp_);
  PCSetFromOptions(pc_);

#endif
}
//****************************************************************************80
void PetscSolver::DestroyPetscObjects()
{
#ifdef PETSC
  KSPDestroy(&ksp_);
  MatDestroy(&Amat_);
  VecDestroy(&bvec_);
  VecDestroy(&xvec_);
  PetscFinalize();
#endif
  destroy_objects_ = false;
}
