/*
 * PMA_init.cpp
 *
 *  Created on: 11 de mar de 2016
 *      Author: rogerio
 */

#include "ParMeshAdapt.h"

ParMeshAdapt::ParMeshAdapt(){
}

ParMeshAdapt::ParMeshAdapt(int argc, char** argv){
	MPI_Init(&argc, &argv);
}

ParMeshAdapt::~ParMeshAdapt(){
	finalize();
}

void ParMeshAdapt::initialize(){
#ifdef __DEBUG_STEPS__
	cout << "initialize: START\n";
#endif

	t_start = MPI_Wtime();

	list_ElmFromFile *pElmList; pElmList=0;
	if (_dim==2){
		pElmList = &tmpTrianglesList;
	}
	else{
		pElmList = &tmpTetraList;
	}

	// number of original element mesh. This will be the base structure
	// ----------------------------------------------------------------------------------------------------------------
	nElem = (int)pElmList->size();

	// Allocate memory for pBase_EDS: That's the first step to build the Main data structure
	// ----------------------------------------------------------------------------------------------------------------
	create_main_datastruct(pElmList);
	pElmList->clear();
	pElmList = 0;

	// if 3-D mesh, a data struct for boundary faces are required
	// ----------------------------------------------------------------------------------------------------------------
	if (_dim==2){
		tmpTrianglesList.clear();
	}

	// the edge data struct will be create
	// ----------------------------------------------------------------------------------------------------------------
	pEdgeDB = new TEDGEDB;
	IDX_str idx;
	idx.base = 0;
	idx.level = 0;
	idx.leaf = 0;
	// ----------------------------


	// create edge data structure for the base mesh
	// start edge data structure passing edge from geometry which may hold physical properties
	std::list<EdgeFromFile>::iterator edge_it;
	if (_dim==2){
		for(edge_it = tmpEdgeList.begin(); edge_it!=tmpEdgeList.end();edge_it++){
			createEdge((*edge_it).id1,(*edge_it).id2,(*edge_it).physical,(*edge_it).geom,true,false,-1,idx);
		}
		createEdgeDataStructure(0);
	}
	else{
		for(edge_it = tmpEdgeList.begin(); edge_it!=tmpEdgeList.end();edge_it++){
			createEdge((*edge_it).id1,(*edge_it).id2,(*edge_it).physical,(*edge_it).geom,true,&idx);
		}
		createEdgeDataStruct_3D();
	}
	tmpEdgeList.clear();

#ifdef __DEBUG_STEPS__
	cout << "initialize: END\n";
#endif
}

void ParMeshAdapt::cleanupMemory(){
	cout << "ParMeshAdapt::cleanupMemory(): START\n";

	for (IdxIter it = idxlist.begin(); it != idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element* elm = &pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		delete elm; elm = 0;
		delete[] pBase_EDS[idx->base].pLeaves;
		pBase_EDS[idx->base].pLeaves = 0;
	}
	delete[] pBase_EDS; pBase_EDS = 0;


	// edge data base where bdry data must be transfered from
	for (auto iter1 = pEdgeDB[0].begin(); iter1 != pEdgeDB[0].end(); iter1++){
		for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
			EdgeInfo* einfo = iter2->second;
			if (_dim==2){
				if (einfo->pNeighbor_1){
					delete[] einfo->pNeighbor_1;
					einfo->pNeighbor_1 = 0;
				}
				if (einfo->pNeighbor_2){
					delete[] einfo->pNeighbor_2;
					einfo->pNeighbor_2 = 0;
				}
			}
			delete einfo; einfo=0;
		}
		iter1->second.clear();
	}
	pEdgeDB[0].clear();

	VertexInfo *vertex;
	std::map<int,double> volume;
	for(VIter vit = VertexDB.begin(); vit!=VertexDB.end(); vit++){
		vertex = vit->second;
		delete[] vertex->coords; vertex->coords = 0;
		vertex->volume.clear();
	}
	VertexDB.clear();

	cout << "ParMeshAdapt::cleanupMemory(): END\n";
}

void ParMeshAdapt::create_main_datastruct(list_ElmFromFile *pElmList){
	int base_idx = 0;
	pBase_EDS = new ElementBase[nElem];
	for (list_ElmFromFile_Iter it=pElmList->begin(); it!=pElmList->end(); it++){

		ElementFromFile *elmff = &(*it);

		// allocate the level array (max: 3 levels of refinement (0,1,2 and 3)
		// At the moment, only the elements of level 0 will be allocated. In this case, only one element!
		pBase_EDS[base_idx].pLeaves = new Element*[4];

		// where the element of the original mesh will be stored.
		// the other levels, initialize as NULL pointers
		pBase_EDS[base_idx].pLeaves[0] = new Element[1];
		pBase_EDS[base_idx].pLeaves[1] = 0;
		pBase_EDS[base_idx].pLeaves[2] = 0;
		pBase_EDS[base_idx].pLeaves[3] = 0;

		// set base element with properties. All leaves, including root element, inherit them
		pBase_EDS[base_idx].geom = elmff->geom;
		pBase_EDS[base_idx].physical = elmff->physical;
		pBase_EDS[base_idx].topLevel = 0;

		// in level=root, there is only one element (it's obvious!)
		int level = 0;
		int leaf = 0;

		// create the list to structure indices. The computational mesh will be accessed through it.
		IDX_str idx;
		idx.base = base_idx;
		idx.level = 0;
		idx.leaf = 0;
		idxlist.push_back(idx);

		// the root element now is created!
		Element *elm = &pBase_EDS[base_idx].pLeaves[level][leaf];
		elm->id1 = (*it).id1;
		elm->id2 = (*it).id2;
		elm->id3 = (*it).id3;
		if (_dim==3){
			elm->id4 = (*it).id4;
		}

		//		elm->toRefine = false;
		//		elm->toUnRefine = false;
		//		elm->toCompMesh = true;
		//		elm->numrefine = 0;
		//		elm->levelAboveMe = false;
		//		elm->green = false;

		base_idx++;
	}
}

void ParMeshAdapt::create_triangle_datastruct(list_ElmFromFile *pElmList){

//	typedef std::map<int,Element> map_third;
//	typedef std::map<int,map_id3> map_second;
//	typedef std::map<int,map_id2> map_main;
//
//
//	int id1, id2, id3;
//	map_main theMap;
//	map_second map2;
//	map_third map3;
//	for(list_ElmFromFile_Iter it = pElmList->begin(); it!= pElmList->end(); it++){
//
//		int ID[3] = {(*it).id1, (*it).id2, (*it).id3};
//		sort(ID,ID+3);
//
//		map3.insert( std::pair<int,Element> >([ID[2]],*it) );
//		map2.insert( std::pair<int,map_third> >([ID[2]],map3) );
//		map2[ID[1]] = map2;
//
//
//	}



}

void ParMeshAdapt::finalize(){
	t_end = MPI_Wtime();

	cout << "Total time: " << t_end - t_start << endl;

	MPI_Finalize();
}
