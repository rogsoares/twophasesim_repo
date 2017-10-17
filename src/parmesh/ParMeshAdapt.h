/*
 * ParMeshAdapt.h
 *
 *  Created on: 29 de fev de 2016
 *      Author: rogerio
 */

#ifndef SRC_PARMESHADAPT_H_
#define SRC_PARMESHADAPT_H_

//namespace PMA{

#include "PMA_Element.h"

static int call_func_getedge = 0;
static int call_func_findedge = 0;

typedef std::list<ElementFromFile> list_ElmFromFile;
typedef std::list<ElementFromFile>::iterator list_ElmFromFile_Iter;


class ParMeshAdapt{

public:

	ParMeshAdapt();
	ParMeshAdapt(int argc, char** argv);
	~ParMeshAdapt();

	// I/O functions
	// -------------------------------------------------------------------------------------------
	void read(const char* filename);
	void write();
	void write(const char* filename);
	void statistics();

	// vertex functions
	// -------------------------------------------------------------------------------------------
	void createVertex(int ID, VertexInfo* vinfo);
	void createVertex(int ID, double x, double y, double z, int physical=0, int geom=0);
	void setVertex(int ID, int physical, int geom);
	void getVertex(int ID, VertexInfo** vinfo);
	void updateNumVertices();
	void getNumVertices(int &n);

	// edge functions
	// -------------------------------------------------------------------------------------------
	void createEdgeDataStructure(int level);
	void createEdgeDataStructure();
	bool createEdge(int id0, int id1, int physical, int geom, bool bdry, bool toBeSplit, int MPV_id, IDX_str& idx);
	void create_hangnode_edge(int id1, int id2,TEDGEDB &edgeDB, int* index);

	void getEdge(int id0, int id1, EdgeInfo** einfo);

	void createEdgeDataStruct_3D();
	void createEdge(int id0, int id1, int physical, int geometry, bool bdry, IDX_str* idx);
	void getEdge(int id0, int id1, EDGE_CREATION& edge_creation, EdgeIter* iter);

	void findEdge(int id0, int id1, bool& found);
	void findEdge(int id0, int id1, bool& found0, bool& found1, TEDGEDB::iterator& iter_out);
	void getNumEdges(int &n);
	void printEdges(TEDGEDB &edgeDB);
	void find_and_remove_edge(int id0, int id1, int level);
	void update_greenedge_neighbors();
	void update_greenedge_neighbors2(int id1, int id2, IDX_str& idx);
	void kill_greenEdge_neighbors(int id1, int id2, int level);
	void cleanEdgeDataStructure();

	// element functions
	// -------------------------------------------------------------------------------------------
	void getElement(IDX_str &idx, Element **pElm);
	void getNumElementsOfComputationalMesh(int &n) const;
	void printElm(std::list<IDX_str> &list);


	void OneLevelSmoothing();
	void OneLevelSmoothing_part1();
	void OneLevelSmoothing_part2();


	void getNumElements(int level, int& nElem);
	void getVertices(const int*, int, VertexInfo**);
	void allocateMemory(int elem_base, int level);
	void smooth_split();

	// it initiates the process of element subdivisions
	void refine();
	void splitElement(Element* pElm, int base, int Level, int leaf_parent, int &max_vertex_ID);
	void splitGreen(int* greenElmIDs, int base, int level, int leaf_parent);
	void RefinementSmoothing();
//	void RefinementSmoothing(int &count);
//	void RefinementSmoothing_part2(int &count);
	void RefinementSmoothing_part3();


	// Unrefine functions:
	// ----------------------------------------------------------
	void unrefine();
	void unrefine(int level);
	//void removeGreenLeaves();
	void setElementsTobeRemoved();

	/*
	 * red-refined triangles will produce extra points along an edge of neighboring triangles that are not refined to the same level.
	 * These triangles must be refined irregularly, which are labeled as green triangles and are of lower quality (shown in Figure 2.15).
	 * Green-refinement is only performed if a single point is inserted. If more points are inserted or higher refinement of this triangle
	 * is wanted, the green triangles are removed and the parent triangle is red-refined, too.
	 */
	void buildGreenElements();

	void calculateMemoryUsage();

	// free memory deleting all pointers to unused mesh elements (triangles and edges).
	// It's supposed that after removing some elements (unrefining process) from mesh, the
	// memory allocated for them still remains. It's means that those pointers are "hidden"
	// from element loop.
	void cleanupMemory();

	// transfer data from temporarily list to the main element structure
	void initialize();
	void finalize();

	void getParent(int &index);
	void printEdge(int id1, int id2);
	int dim() const { return _dim; }

	void create_main_datastruct(list_ElmFromFile *pElmList);
	void create_triangle_datastruct(list_ElmFromFile *pElmList);
	void getBdryFace_list(list_ElmFromFile **pElmList);

	// holds the indices of the main data structure to those elements which form the computational mesh.
	// Always loop this list to reach the computational mesh, which means the top leaves.
	std::list<IDX_str> idxlist;

	ElementBase* pBase_EDS; // pointer: Element Data Structure
							// pEdgeDB[0] -> edge data structure (level: 0) of the original mesh
							// pEdgeDB[1] -> edge data structure (level: 1) for leaves of elements from level 0
							// pEdgeDB[2] -> edge data structure (level: 2) for leaves of elements from level 1
							// pEdgeDB[3] -> edge data structure (level: 3) for leaves of elements from level 2
	TVertexDB VertexDB;		// mesh vertices
	TVertexDB CS_VerticesDB;// Control Surface Vertices Database
	TEDGEDB *pEdgeDB;		// pointer to an array of edge data structure





	// while reading mesh file, store all elements into tmpElementList
	// They will all be transfered to pEDS structure
	std::list<EdgeFromFile> tmpEdgeList;
	list_ElmFromFile tmpTrianglesList;
	list_ElmFromFile tmpTetraList;

	int ** IndexMat;	// a matrix of indices used a leaf access its parent


	// Comput. Mesh Edge DataBase
	TEDGEDB CMEDB;
	int maxVertex_ID;
	int nElem;
	int nVertices;
	int nBdryEdges;
	int nBdryFaces;
	double t_start;
	double t_end;
	char* output_filename;
	ETYPE etype;

	int _dim;
};

#endif /* SRC_PARMESHADAPT_H_ */
