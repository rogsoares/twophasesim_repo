/*
 * KSP_solver.cpp
 *
 *  Created on: 6 de jun de 2016
 *      Author: rogerio
 */


#include "EBFV1_elliptic.h"

namespace PRS{

//

double EBFV1_elliptic::setMatrixFreeOperation(int num_free_nodes,int num_global_nodes, int dim){
	debug_msg("EBFV1_elliptic::setMatrixFreeOperation(): START");

	MatCreateShell(PETSC_COMM_WORLD,num_free_nodes,num_free_nodes,num_free_nodes,num_free_nodes,pDStruct,&matrix);
	MatSetFromOptions(matrix);
	MatShellSetOperation(matrix, MATOP_MULT,(void(*)(void))&EBFV1_elliptic::MatMultUser);

	// use last pressure solution as guess solution for iterative solver and then reduce number of iterations

	// call PETSc to solver system of equation
	PetscInt its;
	PetscBool guessNonZero = PETSC_FALSE;
	KSP ksp;
	PC preconditioner;

	static bool non_zero = false;
	if (non_zero){
		guessNonZero = PETSC_TRUE;
	}

	double startt = MPI_Wtime();
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators(ksp,matrix,pDStruct->G_freenodes);
	KSPSetType(ksp,KSPBCGS);
	KSPGetPC(ksp,&preconditioner);
//	PCSetType(preconditioner,PCHYPRE);
	PCSetType(preconditioner,PCILU);
//	PCHYPRESetType(preconditioner,"ams");
	KSPSetInitialGuessNonzero(ksp,guessNonZero);
	KSPSetTolerances(ksp,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	KSPSetFromOptions(ksp);
	KSPSolve(ksp,pDStruct->RHS,pDStruct->solution);
	KSPGetIterationNumber(ksp,&its);
	KSPDestroy(&ksp);
	MatDestroy(&matrix);
	double endt = MPI_Wtime();

	non_zero = true;

	cout << setprecision(8);
	cout << "\nKSP: solver[" << endt - startt;
	cout << "]\tits[" << its << "]";
	//printVectorToFile(pDStruct->solution,"solution.txt");//exit(1);

	debug_msg("EBFV1_elliptic::setMatrixFreeOperation(): END");
	return 0;
}
}
