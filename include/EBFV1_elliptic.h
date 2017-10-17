#ifndef _EBFV1_ELLIPTIC_EQUATION_H_
#define _EBFV1_ELLIPTIC_EQUATION_H_

#include "MeshData.h"
#include "Physical.h"


/*
 * --------------------------------------------------------------------------------------------
 * Data_struct: used to carry on user's defined matrices and vectors when Krylov solver is called.
 *              PETSc uses these matrices and vector to solve: A*x=b.
 *              Here, A = E*F + G, and it leads to high computational cost. To work around this
 *              issue, only matrix-vector products are conducted via Data_struct
 * --------------------------------------------------------------------------------------------
 */
struct Data_struct{
	Mat G;						// geometric coefficients matrices
	Mat* E;
	Mat* F;
	Mat G_freenodes;
	Mat* F_freenodes;
	Mat* E_freenodes;
	Vec solution;
	Vec RHS;
	Vec z;

	int ndom;
};

//        PRS: Petroleum Reservoir Simulator
namespace PRS{

	// MAS - Matrices Assembly Support
	// holds all Gij and Dij coefficients for entire simulation without mobility term auxiliary
	// Gij and Dij coefficients do not need to be recomputed every time elliptic solver is required
	struct MAS{
		double** Gij;			// matrix nedges by 4, where nedges means all mesh edges
		double** Eij;			// matrix nedges by 6, where nedges means all mesh edges
		int** indices;			// matrix indices used for assembling (global vertex ID)
		double* edge_lambda;	// average mobility for IJ edge: edge_lambda[ith_edge] = (lambda_I + lambda_J)/2
	};


	int MatMultUser(Mat mat, Vec u, Vec y);

	class EBFV1_elliptic{
	public:

		EBFV1_elliptic();
		~EBFV1_elliptic();
		double solver(MeshData *pMData, GeomData *pGDCata, SimulatorParameters *pSimPar, Physical *pPPData);

	private:

		double KSP_solver(Mat A, Mat pcMatrix, Vec y, Vec x, int Krylov_restart, PetscBool guessNonZero, KSPType ksptype, PCType pctype, PetscInt &its);

		// Geometric coefficient functions
		int divergence_E(const double *Cij, int edge, int dom, int dom_flag, int idx0_global, int idx1_global, int id0, int id1, int dim, int counter, double *versor, const double *K, bool is_K_Isotropic);
		int divergence_G(const double *Cij, int edge, int dom, int dom_flag, int idx0_global, int idx1_global, int id0, int id1, int dim, int,
				         double length, double *versor, const double *K, bool is_K_Isotropic);
		int gradient_F_edges(Mat F, const double *Cij, int dom, int idx0, int idx1, int id0, int id1, int dim, GeomData* pGCData);
		int gradient_F_bdry(Mat, int dom, GeomData* pGCData, SimulatorParameters* pSimPar);
		int F_bdryFaces( Mat, const int&);
		int F_bdryEdges(int,Mat);
		void initialize_MAS(GeomData* pGCData);
		void finalize_MAS(GeomData* pGCData);

		void update_edges_mobility(GeomData* pGCData, Physical *pPPData);
		int matmult_mobility(GeomData* pGCData);
		int G_assembly(Mat, GeomData* pGCData);
		int E_assembly(Mat, GeomData* pGCData, int, int&);
		int assembly_global_matrices(int ndom);
		int allocate_matrices_vectors(MeshData* pMData, int ndom, int dim);
		int build_global_matrices(GeomData* pGCData, SimulatorParameters* pSimPar);
		int build_rhs_vector(MeshData* pMData, int ndom);
		int extract_freenodes_submatrices(MeshData* pMData, int ndom);
		int free_matrices(GeomData* pGCData);

		// Associate to mesh nodes new pressure values computed
		double updatePressure(MeshData *pMData, GeomData *pGCData, Physical *pPPData, SimulatorParameters *pSimPar);

		// gradient functions
		void calculatePressureGradient(GeomData *pGCData, Physical *pPPData);
		void calc_p_grad_1(int dom, GeomData *pGCData, Physical *pPPData);
		void calc_p_grad_2(int dom, GeomData *pGCData, Physical *pPPData);
		void calc_p_grad_3(GeomData *pGCData, Physical *pPPData);
		void resetPressureGradient(GeomData *pGCData, Physical *pPPData);

		// PETSc matrix for matrix-free procedure
		Mat matrix;

		// Pointer to be used inside matrix-free procedure
		double setMatrixFreeOperation(int num_free_nodes,int num_global_nodes, int dim);

		/* Group set of calls to compute and assembly matrices E,F and G*/
		double assembly_EFG_RHS(MeshData *pMData, GeomData *pGDCata, Physical *pPPData, SimulatorParameters *pSimPar);

		//set well contribution to right hand side
		int wells(GeomData* pGCData, SimulatorParameters* pSimPar, MeshData* pMData);

		Data_struct *pDStruct;
		MAS* pMAS;					// pointer for Matrix Assembly Support (MAS) struct
		bool Perform_Assembling;		// allows matrix assembly just once

		/*
		 * --------------------------------------------------------------------------------------------
		 * MatMultUser is called several time by PETSc during KSP solver until convergence be reached.
		 * --------------------------------------------------------------------------------------------
		 */

		static int MatMultUser(Mat mat, Vec u, Vec y){
			void *ctx;
			MatShellGetContext(mat,&ctx);
			Data_struct* mats = (Data_struct*)(ctx);

			// the desired multiplication => A*u = y, where A = EF+G
			MatMult(mats->G_freenodes,u,y);							// y = G*u;

			// step 2:  y = y + (E*F*u)_dom1 + (E*F*u)_dom2 + ... + (E*F*u)_domN
			for (int i=0; i<mats->ndom; i++){
				VecZeroEntries(mats->z);
				MatMult(mats->F_freenodes[i],u,mats->z);			// z = F*u
				MatMultAdd(mats->E_freenodes[i],mats->z,y,y);		// y = y + E*z
			}
			return 0;
		}
	};
}
#endif
