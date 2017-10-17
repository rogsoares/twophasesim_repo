/*
 * Element.h
 *
 *  Created on: 29 de fev de 2016
 *      Author: rogerio
 */

#ifndef SRC_ELEMENT_H_
#define SRC_ELEMENT_H_

#include "include.h"

typedef double Coords;

enum ETYPE {TRIANGLE, QUAD, TETRA};

struct IDX_str{
	int base;
	int level;
	int leaf;
};

struct VertexInfo{
	double* coords;
	int physical;
	int geom;
	bool ghost;
	std::map<int,double> volume;
	bool toCompMesh;
	int ID_mapped;		// this is not the real vertex ID. It's used to output.
};

enum EDGE_CREATION {EDGE_EXIST, EDGE_NOT_EXIST, EDGES_NODE_NOT_EXIST};

struct EdgeInfo{
	bool markedToSplit; 	// true: a new vertex has been assigned to edge's middle point
						// It means it's ready to be split.
	int MPV_id;			// Middle Point Vertex ID: for refinement purposes
	int physical;		// physical flag: dirichlet, neumann, etc
	int geom;			// geometry identification
	int numRCopies;		// parallel computing only
	int numFaces;		// number of faces around edge
	bool bdry;			// true: boundary edge, false: internal
	bool toCompMesh; 	// if true, says to element to be included into the computational mesh again
	bool green;

	/*
	 for 2-D meshes, an edge is shared by two elements. Each edge must know who are its
	 neighbor elements. pNeighbor_1 e pNeighbor_2 points to an array of integers that
	 say where those elements are located into the mesh data structured.

	 pNeighbor_1 = {base,level,jthLeave}
	 Means:
				base		: i_th root element
				level	: 0 (root), 1, 2 or 3
				jthleave	: first index of the j_th element required

	 Ex.: base = 100, level = 2, jthleave=6

	 for level 2, there are 2^(2*level) elements (leaves): 16
	 for jhleave=6, we want the third element (leaf):

     leave		indices
       1         0,1,2
       2			 3,4,5
       3			 6,7,8
      ...
	   16		12,13,14


	   To get the connectivities (ID1, ID2, ID3) of that leaf:

	   jthleave   -> ID1
	   jthleave+1 -> ID2
	   jthleave+2 -> ID3

	   base = 100 means that the required element is a leaf of the 100th element of the
	   original mesh.
	*/

	int* pNeighbor_1;
	int* pNeighbor_2;

	// set_of_Neighbos: store a set of pointers that point to all tetrahedrals sharing the same edge
	// it's not applied to 2-D mesh where the number of elements are always one (external boundary edges) or two (internal edges).
	std::set<IDX_str*> set_of_Neighbors;


	void* pEdgeStuffes;
};

struct EdgeFromFile{
	int id1;
	int id2;
	int physical;
	int geom;
};

struct ElementFromFile{
	int id1;
	int id2;
	int id3;
	int id4;
	int physical;
	int geom;
	ETYPE etype;
};

struct Element{
	int id1;
	int id2;
	int id3;
	int id4;
	int ID;
	int numrefine;		// how many times element must be split: max=3 (generates (2^(2n) = 64 leaves)
						// which level its leaves will reach

	int numunrefine;	// how many times element must be split: max=3 (generates (2^(2n) = 64 leaves)
	bool toRefine;
	bool toUnRefine;
	bool toCompMesh; // if true, says to element to be included into the computational mesh again
	bool levelAboveMe;
	bool green;
	int count;
	ETYPE etype;

	// combinations are: id1-id3, id2-id3
	int green_id1;
	int green_id2;
	int childLeaf_0;
	int childLeaf_1;

	double* pCentroid;
	double** pFaceCentroid;
	double volume;
};

struct ElementBase{

	Element** pLeaves;	// a matrix of elements
					// each row contains the number of expected elements for that level
					// pEDB[0][0]	 first row : one element only (original mesh)
					// pEDB[1][0:3]	 second row: four elements
					// pEDB[2][0:15]	 third row : sixteen elements
					// pEDB[3][0:63] fourth row:	 sixty four elements
	int topLevel;

	// all leaves should inherit the element base properties
	int physical;
	int geom;
};



void middlePoint(const double* p1, const double* p2, double* mp);
void centroid(const double** pointers, double* pCentroid, int number_of_pointer);
void calculate_area(double* I, double* J, double* K, double &area);


// debugging tools
void printElement(Element* elm);
void printElement(Element* elm, int withVertexID);


// edge data base
typedef std::map<int, std::map<int, EdgeInfo*> > TEDGEDB;
typedef std::map<int, std::map<int, EdgeInfo*> >::iterator EdgeIter;
typedef std::map<int, EdgeInfo*>::iterator EdgeIter2;


typedef std::map<int, VertexInfo*>::iterator VIter;
typedef std::map<int, VertexInfo*> TVertexDB;
typedef std::map<int, VertexInfo*>::iterator TVertexDBIter;
typedef std::list<IDX_str>::iterator IdxIter;



#endif /* SRC_ELEMENT_H_ */
