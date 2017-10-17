#ifndef GEOMETRICCOEFFICIENTSDATA_H_
#define GEOMETRICCOEFFICIENTSDATA_H_

#include "Matrix.h"
#include "ParMeshAdapt.h"

// where edges' coefficients will be stored
struct DataCoefficients{
	int dom1;	// domains to where Cij/Dij belongs
	int dom2;	// domains to where Cij/Dij belongs
	double* Cij;
	double* Dij;

	int num_subdomains;	// number of sub-domains sharing this edge
	int* sub_domain_list;  // C-style index sub-domains flag

	std::map<int,double*> C_ij;
	double length;

	std::list<int> face_IDs;
};

namespace PRS{

	/*! \class GeomData GeomData.h
	 *  \brief GeomData is designed to set/get set of data associated to mesh entities.
	 *  The set of data is defined by the structure Coefficients above.
	 */
	class GeomData{
	public:

		GeomData();
		~GeomData();

		int ndom;
		int *domains;

		void initialize(ParMeshAdapt* pMesh);
		void calculate_1(ParMeshAdapt* pMesh);
		void calculate_1__3D(ParMeshAdapt* pMesh);
		void calculate_2(ParMeshAdapt* pMesh);
		void calculate_3(ParMeshAdapt* pMesh);
		void calculate_3__3D(ParMeshAdapt* pMesh);

		void transfer_Cij(ParMeshAdapt* pMesh);
		void transfer_Cij_Dij(ParMeshAdapt* pMesh);
		void transfer_Volume(ParMeshAdapt* pMesh);
		void transfer_Mesh(ParMeshAdapt* pMesh);
		void transfer_mesh_3D(ParMeshAdapt* pMesh);

		// mapping functions
		// ---------------------------------------------------------------------------------------------------------
		void mapping(ParMeshAdapt* pMesh);
		void set_sequential_numbering(ParMeshAdapt* pMesh, map_int_int& mapIDtoIndex_global);
		void set_nodeIDs_subdomain(ParMeshAdapt* pMesh, set_int *dom_ids);
		void set_bdry_data(ParMeshAdapt* pMesh, set_int *dom_bdry_ids, int k, map_int_int& mapBdryIDtoIndex);
		void set_bdrynodeIDs_subdomain(ParMeshAdapt* pMesh, set_int *dom_bdry_ids, map_int_int& mapIDtoIndex_global);
		void set_bdrynodeIDs_subdomain_3D(ParMeshAdapt* pMesh, set_int *dom_bdry_ids, map_int_int& mapIDtoIndex_global);
		void set_global_data(ParMeshAdapt* pMesh, set_int *dom_ids, int k, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global);
		void mapping_edges(ParMeshAdapt* pMesh, int k, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global, map_int_int& mapBdryIDtoIndex);
		void mapping_edges_2D(ParMeshAdapt* pMesh, int k, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global, map_int_int& mapBdryIDtoIndex);
		void mapping_edges_3D(ParMeshAdapt* pMesh, int k, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global);
		void mapping_elements(ParMeshAdapt* pMesh, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global);

		void setTotalReservoirVolume(double );
		double getReservoirVolume() const;


		void setMeshDim(int mdim) { dim = mdim; }
		int getMeshDim() const { return dim; }

		void setNumGEdges(int nge) { numGEdges = nge; };
		int getNumGEdges() const { return numGEdges; }

		void setSmallestEdgeLength(double sel) { smallestEdge = sel; }
		double getSmallestEdgeLength() const { return smallestEdge; }

		double getReservoirHeight() const { return reservoirHeight; }
		void setReservoirHeight(double h) { reservoirHeight = h; }

		void getCij(int dom, int row, double* cij);
		void getCij(int dom, int row, const double* &cij);
		void setCij(int dom, int row, double* cij);
		void setCij(int dom, int row, double** cij);
		void getDij(int dom, int row, double* dij);
		void getDij(int dom, int row, const double* &dij);
		void setDij(int dom, int row, double* dij);
		void getDij(int dom, int row, double** dij);

		// calculate: edges length, edges versor (domain and boundary)
		void calculateEdgeProperties(ParMeshAdapt* pMesh);
		void allocatePointers(int,int);
		void deallocatePointers(int);
		int getNumEdgesPerDomain(int i) const;
		int getNumElemPerDomain(int i) const;
		int getNumBDRYEdgesPerDomain(int i) const;
		int getNumBdryFacesPerDomain(int i) const;
		int getNumNodesPerDomain(int i) const;
		int getNumNodes() const;
		void cleanData(ParMeshAdapt* pMesh);
		// needed to create an indexed data structured and though get fast access (direct access, O(1)) to any data.

		// creates a mapping for saturation and pressure solution based in a new mesh (the adapted mesh)
		void mappingNodesIds_Tmp( );

		void getEdge(int dom, int row, const int* &ptr);
		void getEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global);
		void getEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global, int &flag1, int &flag2);

		void getNodeIdx_Global(int dom, int i, int &idx);
		void setVolume(int idx, double v);
		void getVolume(int idx,double &v);
		void getVolume(int dom, int idx, double& vol);
		void getVolume(int dom, int idx_0, int idx_1, double& volumeI, double& volumeJ);

		// return the i_th global node ID
		int getNodeID(int i_th) const;

		void getID(int dom, int idx_0, int idx_1, int& id0, int &id1);
		void getID(int dom, int idx, int& id);
		void getBdryEdge(int dom, int row, int &idx_0, int &idx_1);
		void getBdryEdge(int dom, int row, int &idx_0, int &idx_1, int &idx0_global, int &idx1_global);
		// return
		void getBdryFace(int dom, int row, int &idx_0, int &idx_1, int &idx_2, int &idx0_global, int &idx1_global, int &idx2_global);
		// get control volume for one vertex
		void getBdryVolume(int dom, int idx, double& vol);
		// get control volume for trhee vertices
		void getBdryVolume(int dom, int idx0, int idx1, int idx2, double *vol);
		// get control volume for three vertices
		void getBdryVolume(int dom, int idx0, int idx1, int idx2, double& volumeI, double& volumeJ, double& volumeK);
		// get control volume for two vertices
		void getBdryVolume(int dom, int idx_0, int idx_1, double& volumeI, double& volumeJ);
		// get ID for two vertices on boundary
		void getBdryID(int dom, int idx_0, int idx_1, int& id0, int &id1);
		// get ID for three vertices on boundary
		void getBdryID(int dom, int idx_0, int idx_1, int idx_2, int& id0, int &id1, int& id2);
		// get ID for vertex on boundary
		void getBdryID(int dom, int idx, int& id);
		// get flag number defined by user on .geo file for domain i (i=0,1,2,...,n-1), where n is the number of domains
		int getDomFlag(int i) const;
		// return edge length
		void getLength(int dom, int idx, double &length) const;
		// return versor built over edge pointing from node I to node J, where Node ID I is ALWAYS less than node ID J.
		void getVersor(int dom, int idx, double **v);
		void setCij_norm(int dom, int idx, double val);
		void getCij_norm(int dom, int idx, double &val);
		void getVersor_ExternalBdryElement(int idx, double* versor);

		// Used by Saturation Gradient
		void getExternalBdryEdges(int idx, int &idx0_global, int &idx1_global, int &flag1, int &flag2);
		void getExternalBdryFaces(int idx, int &idx0_global, int &idx1_global, int &idx2_global, int &flag1, int &flag2, int &flag3);


		int getNumExternalBdryEdges() const;
		int getNumExternalBdryFaces() const;

		void calculate_extFaceVersor( );

		void setTotalNumberOfEdges(int n);
		void getTotalNumberOfEdges(int &n) const;
		void setMeshNodes(int n);
		void getMeshNodes(int &n) const;
		void setNumDomains(int n);
		void setDomainList(const int* domlist);
		int getNumDomains() const;
		const int* getDomainList() const;
		int getNumElements() const;
		void getConnectivities(int row, const int **connectivities);
		void getCoordinates(int row, const double** coords);


		// idx is an array of size 2*(dim+1) which stores local and global indices for an element
		// if triangle: 3 locals and 3 globals: total = 2*(2+1) = 6 indices
		// if tetra: 4 locals and 4 globals:    total = 2*(3+1) = 8 indices
		void getElement(int dom, int row, int* idx);

		void getElement(int dom, int row, const int* &idx);

		// calculate number of elements sharing the same vertex for all mesh vertices
		void calculateNumElemSharingVertex();

		// return number of elements sharing a vertex
		int getNumFacesSharingVertex(int elem) const;

		void read_geometry(const char* filename);

		int num_surfaces;
		std::map<int, std::set<int> > surface_map;	// face geometry flag -> volume(s) geometry flags list

		int _ndom;
		int* domainList;
		std::map<int,int> map_domain_flags;
		double reservoirHeight;
		double reservoirVolume;
		int dim;
		double smallestEdge;
		int numGEdges;					// number or global edges
		int* numDomEdges;				// number of edges per sub-domain
		int* numDomElem;				// number of elements per sub-domain
		int* numDomFaces_tmp;			// number of face per sub-domain adapted mesh
		int* numDomBDRYEdges;			// number of bdry edges per sub-domain
		int* numDomBDRYFaces;			// number of bdry faces per sub-domain
		int* numNodesPerDomain;			// number of nodes per sub-domain
		int* numBdryNodesPerDomain;		// number of nodes per sub-domain
		int* numElemSharingVertex;		// number of elements sharing a vertex for all mesh vertices
		int numExtBdryEdges;
		int numExtBdryFaces;
		int numNodes;					// total number of mesh nodes
		double* elem_CDL;				// characteristic dimension length of an element
		double* node_CDL;
		double* elem_HR;					// store element height ration: h_new/h_old

		int* pNodeID;

		int elemtype;					// elemtype = 3 (2-D triangle: 3 nodes), elemtype = 4 (3-D tetrahedron: 4 nodes)
		int numElem;						// number of mesh elements (2D/3D)
		Matrix<int> *pConnectivities;	// matrix for element connectivities
		Matrix<double>* pCoords;			// Mesh node's coordinates: x, y, z

		Matrix<int> *ID;					// node ID per domain
		Matrix<int> *edges;				// edges id0-id1 per domain, where id0 and id1 are array indices and not node IDs.
		Matrix<int> *elem;				// stores local (per domain) and global indices for face node's IDs
		Matrix<int> *faces_tmp;			// (for interpolation):stores local (per domain) and global indices for face node's IDs
		Matrix<int> *nodes;				// same as edges
		Matrix<double>* volume;			//node volumes per domain
		Matrix<double>* volume_global;

		Matrix<int> *ID_bdry;			// node ID per domain
		Matrix<int> *edges_bdry;			// edges id0-id1 per domain, where id0 and id1 are array indices and not node IDs.
		Matrix<int> *faces_bdry;
		Matrix<double>* volume_bdry;		//node volumes per domain

		Matrix<double>* versor_ExtBdryElem;	// versor for external elements (2-D: edges, 3-D: triangles)
		Matrix<int>* external_bdry_elem;		// external_bdry_elem: idx0_global, idx1_global, idx2_global, flag1, flag2, flag3

		Matrix<double>* Cij;				// Cij vector
		Matrix<double>* Dij;

		Matrix<double>* edge_versor;
		Matrix<double>* edge_length;
		Matrix<double>* Cij_norm;
	};
}

#endif /*GEOMETRICCOEFFICIENTSDATA_H_*/

