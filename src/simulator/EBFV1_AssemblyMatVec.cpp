/*
 * EBFV1_AssemblyMatVec.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: rogerio
 */


#include "EBFV1_elliptic.h"

namespace PRS{

	double EBFV1_elliptic::assembly_EFG_RHS(MeshData *pMData, GeomData *pGCData, Physical *pPPData, SimulatorParameters *pSimPar){
		debug_msg("EBFV1_elliptic::assembly_EFG_RHS(): START");

		if (Perform_Assembling){
			allocate_matrices_vectors(pMData,pGCData->getNumDomains(),pGCData->getMeshDim());
			initialize_MAS(pGCData);
			build_global_matrices(pGCData,pSimPar);
			Perform_Assembling = false;
		}

		// this multiplication will be performed every new time step.
		update_edges_mobility(pGCData,pPPData);

		// multiply global matrices E and G by mobility
		matmult_mobility(pGCData);

		// assembly global G_tmp and E. F matrix has already been assembled
		assembly_global_matrices(pGCData->getNumDomains());

		// get sub matrix from global matrix for system of equations
		extract_freenodes_submatrices(pMData,pGCData->getNumDomains());

//		printMatrixToFile(pDStruct->G,"G_matrix.txt");
//		printMatrixToFile(pDStruct->G_freenodes,"G_freematrix.txt");
//		for (int i=0; i<pGCData->getNumDomains(); i++){
//			char fname1[128]; sprintf(fname1,"E_matrix_%d.txt",i);
//			char fname2[128]; sprintf(fname2,"E_freematrix_%d.txt",i);
//			char fname3[128]; sprintf(fname3,"F_matrix_%d.txt",i);
//			char fname4[128]; sprintf(fname4,"F_freematrix_%d.txt",i);
//			printMatrixToFile(pDStruct->E[i],fname1);
//			printMatrixToFile(pDStruct->E_freenodes[i],fname2);
//			printMatrixToFile(pDStruct->F[i],fname3);
//			printMatrixToFile(pDStruct->F_freenodes[i],fname4);
//		}

		// build RHS vector based on dirichlet nodes
		build_rhs_vector(pMData,pGCData->getNumDomains());

		MatDestroy(&pDStruct->G);
		for (int i=0; i<pGCData->getNumDomains(); i++){
			MatDestroy(&pDStruct->E[i]);
			MatDestroy(&pDStruct->F[i]);
		}

		// complete RHS building adding well flow rate (Neumann)
		wells(pGCData,pSimPar,pMData);

		debug_msg("EBFV1_elliptic::assembly_EFG_RHS(): END");
		return 0;
	}

	int EBFV1_elliptic::assembly_global_matrices(int ndom){
		debug_msg("EBFV1_elliptic::assembly_global_matrices(): START");

		MatAssemblyBegin(pDStruct->G,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(pDStruct->G,MAT_FINAL_ASSEMBLY);
		for (int dom=0; dom<ndom; dom++){
			MatAssemblyBegin(pDStruct->E[dom],MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(pDStruct->E[dom],MAT_FINAL_ASSEMBLY);
			MatAssemblyBegin(pDStruct->F[dom],MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(pDStruct->F[dom],MAT_FINAL_ASSEMBLY);
		}

		debug_msg("EBFV1_elliptic::assembly_global_matrices(): END");
		return 0;
	}

	int EBFV1_elliptic::allocate_matrices_vectors(MeshData* pMData, int ndom, int dim){
		debug_msg("EBFV1_elliptic::allocate_matrices_vectors: START");

		int num_global_nodes, num_free_nodes;
		pMData->get_num_global_nodes(num_global_nodes);
		pMData->get_num_free_nodes(num_free_nodes);
//		cout << "num_global_nodes : " << num_global_nodes << endl;
//		cout << "num_free_nodes   : " << num_free_nodes << endl;
//		cout << "ndom             : " << ndom << endl;
//		cout << "dim              : " << dim << endl;
//		exit(1);

		pDStruct->ndom = ndom;
		pDStruct->E = new Mat[ndom];
		pDStruct->F = new Mat[ndom];
		pDStruct->F_freenodes = new Mat[ndom];
		pDStruct->E_freenodes = new Mat[ndom];
		MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,num_global_nodes,num_global_nodes,100,PETSC_NULL,100,PETSC_NULL,&pDStruct->G);
		for (int dom=0; dom<ndom; dom++){
			MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,num_global_nodes,num_global_nodes*dim,100,PETSC_NULL,100,PETSC_NULL,&pDStruct->E[dom]);
			MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,num_global_nodes*dim,num_global_nodes,100,PETSC_NULL,100,PETSC_NULL,&pDStruct->F[dom]);
		}
		VecCreate(PETSC_COMM_WORLD,&pDStruct->RHS);
		VecSetSizes(pDStruct->RHS,PETSC_DECIDE,num_free_nodes);
		VecSetFromOptions(pDStruct->RHS);
		VecDuplicate(pDStruct->RHS,&pDStruct->solution);
		VecCreate(PETSC_COMM_WORLD,&pDStruct->z);
		VecSetSizes(pDStruct->z,PETSC_DECIDE,num_global_nodes*dim);
		VecSetFromOptions(pDStruct->z);

		debug_msg("EBFV1_elliptic::allocate_matrices_vectors: END");
		return 0;
	}

	void EBFV1_elliptic::update_edges_mobility(GeomData* pGCData, Physical *pPPData){
		debug_msg("EBFV1_elliptic::update_edges_mobility(): START");

		int i, idx_I, idx_J;
		double Sw_0, Sw_1;

		int nedges = 0;
		int ndom = pGCData->getNumDomains();
		for(i=0; i<ndom; i++){
			nedges += pGCData->getNumEdgesPerDomain(i);
		}

		for (i=0; i<nedges; i++){
			idx_I = pMAS->indices[i][0] - 1;					// it return global vertex ID
			idx_J = pMAS->indices[i][1] - 1;					// it return global vertex ID
			pPPData->getSaturation(idx_I,Sw_0);
			pPPData->getSaturation(idx_J,Sw_1);
			pMAS->edge_lambda[i] = 1.0;//0.5*(pPPData->getTotalMobility(Sw_0) + pPPData->getTotalMobility(Sw_1));
		}

		debug_msg("EBFV1_elliptic::update_edges_mobility(): END");
	}

	int EBFV1_elliptic::matmult_mobility(GeomData* pGCData){
		debug_msg("EBFV1_elliptic::matmult_mobility(): START");

		int i, j, k, counter;
		int ndom = pGCData->getNumDomains();
		int dim = pGCData->getMeshDim();
		int idxm[2], idxn[2*dim], pos1, pos2;
		int ncols = 2*dim;
		double Eij[4*dim];
		double Gij[4];

		counter = 0;
		for(i=0; i<ndom; i++){
			int nedges = pGCData->getNumEdgesPerDomain(i);
			for(j=0; j<nedges; j++){

				idxm[0] = pMAS->indices[counter][0] - 1;
				idxm[1] = pMAS->indices[counter][1] - 1;

				pos1 = dim*idxm[0];
				pos2 = dim*idxm[1];

				for (k=0; k<dim; k++){
					idxn[k] = pos1+k;
					idxn[dim+k] = pos2+k;
				}

				for(k=0; k<4; k++){
					Gij[k] = pMAS->Gij[counter][k]*pMAS->edge_lambda[counter];
				}

				//cout << "dom: " << i << "|\t" << idxm[0] << " " << idxm[1] << endl;
				for(k=0; k<4*dim; k++){
					Eij[k] = pMAS->Eij[counter][k]*pMAS->edge_lambda[counter];
					//cout << Eij[k] << " ";
				}
				//cout << endl;
				counter++;
				MatSetValues(pDStruct->G,2,idxm,2,idxm,Gij,ADD_VALUES);
				MatSetValues(pDStruct->E[i],2,idxm,ncols,idxn,Eij,ADD_VALUES);
				//break;
			}
		}

		//cout << "Counter: " << counter << endl;
		debug_msg("EBFV1_elliptic::matmult_mobility(): END");
		return 0;
	}

	int EBFV1_elliptic::extract_freenodes_submatrices(MeshData* pMData, int ndom){
		debug_msg("EBFV1_elliptic::extract_freenodes_submatrices(): START");

		IS pos_freenodes;		// indices for free nodes
		IS pos_extendednodes;	// indices for all nodes times dimension: pos = num_global_nodes*dim
		pMData->get_IS_free_nodes(&pos_freenodes);
		pMData->get_IS_extendednodes_nodes(&pos_extendednodes);

		MatGetSubMatrix(pDStruct->G,pos_freenodes,pos_freenodes,MAT_INITIAL_MATRIX,&pDStruct->G_freenodes);
		for (int i=0; i<ndom; i++){
			MatGetSubMatrix(pDStruct->F[i],pos_extendednodes,pos_freenodes,MAT_INITIAL_MATRIX,&pDStruct->F_freenodes[i]);
			MatGetSubMatrix(pDStruct->E[i],pos_freenodes,pos_extendednodes,MAT_INITIAL_MATRIX,&pDStruct->E_freenodes[i]);
		}

		debug_msg("EBFV1_elliptic::extract_freenodes_submatrices(): END");
		return 0;
	}

	int EBFV1_elliptic::build_rhs_vector(MeshData* pMData, int ndom){
		debug_msg("EBFV1_elliptic::build_rhs_vector(): START");

		VecZeroEntries(pDStruct->RHS);

		// Get contribution from matrices related to dirichlet (prescribed) nodes
		// get G matrix contribution:
		// -------------------------------------------------------------------------------------
		Mat A;										// auxiliary matrix
		Vec Vec_Dirichlet;
		IS pos_freenodes;							// indices for free nodes
		IS pos_dirichletnodes;						// indices for all nodes times dimension: pos = num_global_nodes*dim
		PetscScalar *a; a = 0;
		int i, num_free_nodes;

		pMData->get_num_free_nodes(num_free_nodes);
		pMData->get_IS_free_nodes(&pos_freenodes);
		pMData->get_IS_dirichlet_nodes(&pos_dirichletnodes);
		pMData->get_Vec_dirichlet(&Vec_Dirichlet);

		MatGetSubMatrix(pDStruct->G,pos_freenodes,pos_dirichletnodes,MAT_INITIAL_MATRIX,&A);
		MatMult(A,Vec_Dirichlet,pDStruct->RHS);
		MatDestroy(&A);
		VecGetArray(pDStruct->RHS,&a);

		for (i=0; i<num_free_nodes; i++){
			a[i] = -a[i];
		}

		VecRestoreArray(pDStruct->RHS,&a);
		a = 0;

		// get E*F matrix contribution per sub_domains:
		// -------------------------------------------------------------------------------------
		Mat EF;
		Vec RHS_tmp;
		VecDuplicate(pDStruct->RHS,&RHS_tmp);
		for (int dom=0; dom<ndom; dom++){
			MatMatMult(pDStruct->E[dom],pDStruct->F[dom],MAT_INITIAL_MATRIX,1.0,&EF);
			MatGetSubMatrix(EF,pos_freenodes,pos_dirichletnodes,MAT_INITIAL_MATRIX,&A);
			MatDestroy(&EF);
			MatMult(A,Vec_Dirichlet,RHS_tmp);
			MatDestroy(&A);
			VecAXPY(pDStruct->RHS,-1.0,RHS_tmp);
			VecZeroEntries(RHS_tmp);
		}
		VecDestroy(&RHS_tmp);

		debug_msg("EBFV1_elliptic::build_rhs_vector(): END");
		return 0;
	}
}
