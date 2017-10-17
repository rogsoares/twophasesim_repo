/*
 * createEdgeDataStruct_3D.cpp
 *
 *  Created on: 20 de jun de 2016
 *      Author: rogerio
 */


#include "ParMeshAdapt.h"

void ParMeshAdapt::createEdgeDataStruct_3D(){
	cout << "Create edge data structure: START\n";
	double t1 = MPI_Wtime();

	// loop over computational mesh: take only elements at level 'level'
	for (IdxIter it = idxlist.begin(); it != idxlist.end(); it++){
		IDX_str* idx = &(*it);

		Element* elm = &pBase_EDS[idx->base].pLeaves[0][idx->leaf];
		createEdge(elm->id1,elm->id2,0,0,0,idx);
		createEdge(elm->id1,elm->id3,0,0,0,idx);
		createEdge(elm->id1,elm->id4,0,0,0,idx);
		createEdge(elm->id2,elm->id3,0,0,0,idx);
		createEdge(elm->id2,elm->id4,0,0,0,idx);
		createEdge(elm->id3,elm->id4,0,0,0,idx);
	}

	cout << "Create edge data structure: END\n";
	double t2 = MPI_Wtime();

	// ------------------------------------------------------------------------------------------------------------------------
	double frac,h,m,s;
	double t = t2-t1;
	frac = modf(t/3600.,&h);
	frac = modf(frac*60.,&m);
	frac = modf(frac*60.,&s);
	cout << setprecision(0) << "Elapsed Time: " << h << "h " << m << "min " << s << "s\n";
}

void ParMeshAdapt::createEdge(int id0, int id1, int physical, int geometry, bool bdry, IDX_str* idx){

	if (id0 > id1){
		swap(id0,id1);			// edge's IDs: id0 > id1 ALWAYS!!!!
	}

	EdgeIter iter;
	EdgeInfo* edge; edge = 0;
	EDGE_CREATION edge_creation;

	// look for the desired edge into database
	getEdge(id0,id1,edge_creation,&iter);

	// if any ID has been found, then create a new edge
	if (edge_creation==EDGE_EXIST){

		std::map<int, EdgeInfo*> map = iter->second;
		edge = map[id1];
		edge->set_of_Neighbors.insert(idx);
	}
	else{
		edge = new EdgeInfo;

		// the three next lines are filled for edges presented into gmsh file as boundary condition elements
		edge->physical = physical;
		edge->geom = geometry;
		edge->bdry = bdry;

		if (!bdry){
			edge->set_of_Neighbors.insert(idx);
		}

		// neither id0 non id1 belongs to edge data struct
		if (edge_creation==EDGE_NOT_EXIST){
			std::map<int,EdgeInfo*> secondPart;
			secondPart.insert( std::pair<int,EdgeInfo*>(id1,edge) );
			(*pEdgeDB)[id0] = secondPart;
		}
		// id0 was found into edge data base, but id1.
		else{
			iter->second.insert( std::pair<int,EdgeInfo*>(id1,edge) );
		}
	}
}

void ParMeshAdapt::getEdge(int id0, int id1, EDGE_CREATION& edge_creation, EdgeIter* iter){

	// we start assuming edge exits
	edge_creation = EDGE_EXIST;

	*iter = pEdgeDB->find(id0);
	if (*iter==pEdgeDB->end()){
		edge_creation = EDGE_NOT_EXIST;
//		*iter = 0;
		return;
	}

	std::map<int,EdgeInfo*> map2 = (*iter)->second;
	EdgeIter2 iter2 = map2.find(id1);
	if (iter2 == map2.end()){
		edge_creation = EDGES_NODE_NOT_EXIST;
	}
}
