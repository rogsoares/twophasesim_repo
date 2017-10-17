#include "GeomData.h"

namespace PRS{

	void GeomData::allocatePointers(int nnodes_global, int dim){
		debug_msg("allocatePointers: START");

		try{

			int ndom = getNumDomains();
			Cij = new Matrix<double>[ndom];
			Cij_norm = new Matrix<double>[ndom];
			volume = new Matrix<double>[ndom];
			edges = new Matrix<int>[ndom];
			elem = new Matrix<int>[ndom];
			nodes = new Matrix<int>[ndom];
			ID = new Matrix<int>[ndom];
			edge_length = new Matrix<double>[ndom];
			edge_versor = new Matrix<double>[ndom];
			volume_bdry = new Matrix<double>[ndom];
			ID_bdry = new Matrix<int>[ndom];
			volume_global = new Matrix<double>[1];
			volume_global[0].allocateMemory(nnodes_global);
			volume_global[0].initialize(.0);


			if (dim==3){
				faces_bdry = new Matrix<int>[ndom];
			}

			if (dim==2){
				edges_bdry = new Matrix<int>[ndom];
				for (int k=0; k<ndom; k++){
					edges_bdry[k].allocateMemory(numDomBDRYEdges[k],6);
					edges_bdry[k].initialize(0);
				}
				Dij = new Matrix<double>[ndom];
			}


			external_bdry_elem = new Matrix<int>[1];

			// 4 = 2 indices + 2 flags // 6 = 3 indices + 3 flags
			external_bdry_elem[0].allocateMemory((dim==2)?numExtBdryEdges:numExtBdryFaces,2*dim);

			versor_ExtBdryElem = new Matrix<double>[1];
			versor_ExtBdryElem[0].allocateMemory((dim==2)?numExtBdryEdges:numExtBdryFaces,3);
			pConnectivities = new Matrix<int>[1];
			pConnectivities->allocateMemory(numElem,elemtype);
			pCoords = new Matrix<double>[1];
			pCoords->allocateMemory(numNodes,3);



			for (int k=0; k<ndom; k++){

				int nedges = numDomEdges[k];
				int nnodes = numNodesPerDomain[k];
				int nbnodes = numBdryNodesPerDomain[k];
				int nelem = numDomElem[k];

//				cerr << "sub-dom :  " << k <<
//						" nedges :  " << nedges <<
//						" nnodes :  " << nnodes <<
//						" nbnodes:  " << nbnodes <<
//						" nelem  :  " << nelem << "\n";

				if (dim==2){
					Cij[k].allocateMemory(nedges,3);
					Cij[k].initialize(.0);
				}
				else{
					//Cij[k].allocate_array_pointers(nedges,3);
					Cij[k].allocateMemory(nedges,3);
					Cij[k].initialize(.0);
				}

				int size = (dim==2)?numDomBDRYEdges[k]:numDomBDRYFaces[k];
				//cerr << "SIZE: " << size << endl;

				if (dim==2){
					Dij[k].allocateMemory(size,3);
					Dij[k].initialize(.0);
				}


				edges[k].allocateMemory(nedges,6);
				edges[k].initialize(0);

				// stores local and global indices for triangles (2-D) or tetrahedra (3-D)
				elem[k].allocateMemory(nelem,2*(dim+1));
				elem[k].initialize(0);

				if (dim==3){
					faces_bdry[k].allocateMemory(size,9);
					faces_bdry[k].initialize(0);
				}

				nodes[k].allocateMemory(nnodes);
				volume[k].allocateMemory(nnodes);
				volume[k].initialize(.0);
				ID[k].allocateMemory(nnodes);
				ID[k].initialize(0);

				volume_bdry[k].allocateMemory(nbnodes);
				volume_bdry[k].initialize(.0);
				ID_bdry[k].allocateMemory(nbnodes);
				ID_bdry[k].initialize(0);

				Cij_norm[k].allocateMemory(nedges);
				Cij_norm[k].initialize(.0);
				edge_length[k].allocateMemory(nedges);
				edge_length[k].initialize(.0);
				edge_versor[k].allocateMemory(nedges,3);
				edge_versor[k].initialize(.0);
			}
		}
		catch  (std::bad_alloc& ba){
			throw Exception(__LINE__,__FILE__,"Memory allocation failed. Bad allocation.");
		}

		//exit(1);
		debug_msg("allocatePointers: END");
	}

	void GeomData::deallocatePointers(int _dim){
		debug_msg("deallocatePointers: START");


		_dim = dim;
		int ndom = getNumDomains();
	//	volume_global[0].freeMemory();
	//	delete[] volume_global; volume_global = 0;
//		if (dim==2){
//			external_bdry_elem[0].freeMemory();
//			for (int k=0; k<ndom; k++){
//				edges_bdry[k].freeMemory();
//			}
//			delete[] external_bdry_elem; external_bdry_elem=0;
//			delete[] edges_bdry; edges_bdry = 0;
//			delete[] numDomBDRYEdges; numDomBDRYEdges = 0;
//		}

		dealloc_INT_vector(pNodeID);

//		versor_ExtBdryElem[0].freeMemory();
//		delete[] versor_ExtBdryElem; versor_ExtBdryElem=0;

		pConnectivities[0].freeMemory();
		delete[] pConnectivities; pConnectivities = 0;

//		pCoords[0].freeMemory();
//		delete[] pCoords; pCoords = 0;


		for (int k=0; k<ndom; k++){
			Cij[k].freeMemory();
			Cij_norm[k].freeMemory();
			edge_length[k].freeMemory();
			Dij[k].freeMemory();
//			//volume[k].freeMemory();
			edges[k].freeMemory();
			elem[k].freeMemory();
			nodes[k].freeMemory();
			ID[k].freeMemory();
			ID_bdry[k].freeMemory();
//			volume_bdry[k].freeMemory();
			edge_versor[k].freeMemory();

			if (dim==3){
				faces_bdry[k].freeMemory();
			}
		}
		delete[] Cij; Cij = 0;
		delete[] numDomEdges; numDomEdges = 0;
		delete[] Cij_norm; Cij_norm = 0;
		delete[] Dij; Dij = 0;
//		//delete[] volume; volume = 0;
		delete[] edges; edges = 0;
		delete[] elem; elem = 0;
		delete[] nodes; nodes = 0;
		delete[] ID; ID = 0;
//		delete[] volume_bdry; volume_bdry = 0;
		delete[] ID_bdry; ID_bdry = 0;
		delete[] edge_length; edge_length = 0;
		delete[] edge_versor; edge_versor = 0;
		if (dim==3){
			delete[] faces_bdry; faces_bdry = 0;
		}

		debug_msg("deallocatePointers: END");
	}
}
