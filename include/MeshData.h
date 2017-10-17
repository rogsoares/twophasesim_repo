#ifndef MESHDATA_H_
#define MESHDATA_H_

#include "SimulatorParameters.h"

namespace PRS{

	class MeshData{
	public:

		MeshData();
		~MeshData();

		void initialize(ParMeshAdapt* pMesh, SimulatorParameters* pSimPar, GeomData *);
		int set_dirichlet_free_nodes_arrays(ParMeshAdapt* pMesh, SimulatorParameters* pSimPar);
		bool isPrescribed(int ID);
		void get_dirichlet_value(int ID, double& val);
		void get_num_global_nodes(int &ngn);
		void get_num_free_nodes(int &nfn);
		void get_num_dirichlet_nodes(int &ndn);
		void get_free_nodes_ptr(int** nfn_ptr);
		void get_dirichlet_nodes_ptr(int **ndn_ptr);
		void get_IS_free_nodes(IS *free_rows);
		void get_IS_dirichlet_nodes(IS *dir_rows);
		void get_IS_extendednodes_nodes(IS *extendednodes_nodes);
		void get_Vec_dirichlet(Vec *v);

		void get_node_condition(int ID, int &pos);	// IN: node's ID	OUT: index position for free or dirichlet
													// Ex.: ID = 278 has production well (Neumann, not Dirichlet).
													// Which RHS row the well flow rate must be inserted into?

	private:

		IS IS_free_rows;
		IS IS_dirichlet_rows;
		IS IS_extendednodes_nodes;
		Vec Vec_Dirichlet;

		int num_free_nodes;
		int num_global_nodes;
		int num_dirichlet_nodes;
		int *free_nodes_ptr;					// free_nodes_ptr[i] = ID - 1
		int *dirichlet_nodes_ptr;			// dirichlet_nodes_ptr[i] = ID - 1
		int *dirichlet_flags;
		double *dirichlet_values_ptr;		// dirichlet_nodes_ptr[i] = ID - 1
		int *free_dirichlet_global;			// a vector of num_global_nodes size used to map node IDs
											// it helps to know if node is free or dirichlet based on its ID.
	};
}
#endif /*MESHDATA_H_*/
