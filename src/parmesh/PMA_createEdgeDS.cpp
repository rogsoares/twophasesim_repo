/*
 * PMA_createEdgeDS.cpp
 *
 *  Created on: 1 de mar de 2016
 *      Author: rogerio
 */


#include "ParMeshAdapt.h"

// Create the edge data structure based on elements from the original mesh
void ParMeshAdapt::createEdgeDataStructure(int level){
	// create edge data structure for level 'level'
	bool created[3];
	int count = 0;

	// loop over computational mesh: take only elements at level 'level'
	for (IdxIter it = idxlist.begin(); it != idxlist.end(); it++){
		IDX_str* idx = &(*it);
		if (idx->level == level){
			Element* elm = &pBase_EDS[idx->base].pLeaves[level][idx->leaf];
			created[0] = createEdge(elm->id1,elm->id2,-1,-1,false,false,-1,*idx);
			created[1] = createEdge(elm->id1,elm->id3,-1,-1,false,false,-1,*idx);
			created[2] = createEdge(elm->id2,elm->id3,-1,-1,false,false,-1,*idx);
			//printElement(elm);
			for (int i=0; i<3; i++){
				if (created[i]) count++;
			}
		}
	}

	//cout << "Number of edges created: " << count << endl;

	// transfer bdry edge data from level-1 to level:
	if (level){

		TEDGEDB::iterator iter;
		EdgeInfo* edge; edge = 0;

		// edge data base where bdry data must be transfered to
		TEDGEDB edgeDB = pEdgeDB[0];

		// edge data base where bdry data must be transfered from
		for (EdgeIter iter1 = pEdgeDB[0].begin(); iter1 != pEdgeDB[0].end(); iter1++){
			for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* einfo = iter2->second;

				// take only split bdry edges
				if (einfo->bdry && einfo->markedToSplit){

					int id1 = iter1->first;
					int id2 = iter2->first;
					int id3 = einfo->MPV_id;

					// edge from level:  level-1 =>  id1 - id2
					// edge from level:  level   =>  id1 - id3 AND id2 - id3

					// get desired edge into database at level 'level'
					getEdge(id1,id3,&edge);

					edge->bdry = true;
					edge->physical = einfo->physical;
					edge->geom = einfo->geom;

					// get desired edge into database at level 'level'
					getEdge(id2,id3,&edge);

					edge->bdry = true;
					edge->physical = einfo->physical;
					edge->geom = einfo->geom;
				}
			}
		}
	}
	//cout << ". Done\n";
}

bool ParMeshAdapt::createEdge(int id0, int id1, int physical, int geom, bool bdry,bool toBeSplit, int MPV_id, IDX_str& idx){

	// edge's IDs: id0 > id1 ALWAYS!!!!
	if (id0 > id1){
		swap(id0,id1);
	}


	bool found0, found1;
	TEDGEDB::iterator iter;
	EdgeInfo* einfo; einfo = 0;

	// look for the desired edge into database
	findEdge(id0,id1,found0,found1,iter);


	//	if (id1==4099 && id0==3590){
	//		cout << "Creating: " << id0 << "  "  << id1 << endl;
	//		cout << "found0: " << found0 << "  found1: "  << found1 << endl;
	//	}

	// if the pair of ID already exits, update the number of elements sharing edge and then return.
	if (found0 && found1){
		getEdge(id0,id1,&einfo);

//		if (id0==201 && id1==784){
//			Element* elm = &pBase_EDS[idx.base].pLeaves[idx.level][idx.leaf];
//			cout << "Creating 2: " << id0 << "  "  << id1 << endl;
//			cout << "found0: " << found0 << "  found1: "  << found1 << endl;
//			cout << "2: idx.base: " << idx.base << "  idx.level: "  << idx.level <<  "  idx.leaf: " << idx.leaf << endl;
//			cout << "elem IDs: " << elm->id1 << "  "  << elm->id2 <<  "  " << elm->id3 << endl;
//			cout << "einfo->pNeighbor_1: " << einfo->pNeighbor_1 << "  einfo->pNeighbor_2 " << einfo->pNeighbor_2 << endl;
//		}

		if (einfo->bdry && !einfo->pNeighbor_1){

			// the another element sharing edge
			einfo->pNeighbor_1 = new int[3];
			einfo->pNeighbor_1[0] = idx.base;
			einfo->pNeighbor_1[1] = idx.level;
			einfo->pNeighbor_1[2] = idx.leaf;
		}
		else{
			einfo->numFaces++;

			// the another element sharing edge
			einfo->pNeighbor_2 = new int[3];
			einfo->pNeighbor_2[0] = idx.base;
			einfo->pNeighbor_2[1] = idx.level;
			einfo->pNeighbor_2[2] = idx.leaf;
		}
		return false;
	}

	// if found0 or found1 are false, the edge does not exit.
	einfo = new EdgeInfo;

	einfo->geom = geom;
	einfo->physical = (bdry)?physical:-1;
	einfo->bdry = bdry;
	einfo->green = false;
	einfo->toCompMesh = true;

	// the line below works when an edge over geometry (boundary) asks to be split into two new ones.
	// In general, edges are created after elements had been split.
	// the variable "toBeSplit" is true when a loop over edges over geometry is performed
	einfo->markedToSplit = toBeSplit;

	// if edge is new, it counts one element.
	einfo->numFaces = 1;

	// element that share the edge
	if (einfo->bdry){
		einfo->pNeighbor_1 = 0;
	}
	else{
		einfo->pNeighbor_1 = new int[3];
		einfo->pNeighbor_1[0] = idx.base;
		einfo->pNeighbor_1[1] = idx.level;
		einfo->pNeighbor_1[2] = idx.leaf;
	}
	einfo->pNeighbor_2 = 0;

//	if (id0==201 && id1==784){
//		Element* elm = &pBase_EDS[idx.base].pLeaves[idx.level][idx.leaf];
//		cout << "Creating: " << id0 << "  "  << id1 << endl;
//		cout << "found0: " << found0 << "  found1: "  << found1 << endl;
//		cout << "2: idx.base: " << idx.base << "  idx.level: "  << idx.level <<  "  idx.leaf: " << idx.leaf << endl;
//		cout << "elem IDs: " << elm->id1 << "  "  << elm->id2 <<  "  " << elm->id3 << endl;
//		cout << "einfo->pNeighbor_1: " << einfo->pNeighbor_1 << "  einfo->pNeighbor_2 " << einfo->pNeighbor_2 << endl;
//	}

	std::map<int,EdgeInfo*> secondPart;
	// first if: id0 already exist into edge DB, but id1.
	// second if: both id0 and id1 do not exist.
	if (found0 && !found1){
		iter->second.insert( std::pair<int,EdgeInfo*>(id1,einfo) );
	}
	else if (!found0 && !found1){
		secondPart.insert( std::pair<int,EdgeInfo*>(id1,einfo) );
		(*pEdgeDB)[id0] = secondPart;
	}

	return true;
}

void ParMeshAdapt::getEdge(int id0, int id1, EdgeInfo** einfo){
	if (id0>id1){
		swap(id0,id1);
	}

	TEDGEDB::iterator iter1;
	iter1 = pEdgeDB->find(id0);
	if (iter1==pEdgeDB->end()){
		cout << "Edge (" << id0 << "," << id1 << ") not found!\n";
		exit(1);
	}

	std::map<int,EdgeInfo*> map2 = iter1->second;
	std::map<int,EdgeInfo*>::iterator iter2 = map2.find(id1);
	if (iter2 == map2.end()){
		cout << "Edge (" << id0 << "," << id1 << ") not found!\n";
		exit(1);
	}
	*einfo = iter2->second;
}

void ParMeshAdapt::findEdge(int id0, int id1, bool& found){
	found = false;
	if (id0>id1){
		swap(id0,id1);
	}
	EdgeIter iter1;
	iter1 = pEdgeDB->find(id0);			// find first edge vertex ID
	if (iter1 != pEdgeDB->end()){			// if found, find second edge vertex ID
		std::map<int,EdgeInfo*>::iterator iter2 = iter1->second.find(id1);
		if (iter2!=iter1->second.end()){
			found = true;
		}
	}
}

void ParMeshAdapt::findEdge(int id0, int id1, bool& found0, bool& found1, TEDGEDB::iterator& iter_out){
	iter_out = pEdgeDB->find(id0);			// find first edge vertex ID
	if (iter_out != pEdgeDB->end()){		// if found, find second edge vertex ID
		found0 = true;
		std::map<int,EdgeInfo*>::iterator iter = iter_out->second.find(id1);
		found1 = (iter!=iter_out->second.end())?true:false;
	}
	else{
		found0 = false;
		found1 = false;
	}
}

void ParMeshAdapt::getNumEdges(int &n){
	n = 0;
	for(EdgeIter it1 = pEdgeDB->begin(); it1!=pEdgeDB->end(); it1++){
		n += (int)it1->second.size();
	}
}

void ParMeshAdapt::update_greenedge_neighbors(){
#ifdef __DEBUG_STEPS__
	cout << "update_greenedge_neighbors(): START\n";
#endif

//	Element* elm; elm = 0;
//	for (EdgeIter iter1 = pEdgeDB[0].begin(); iter1 != pEdgeDB[0].end(); iter1++){
//		for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
//			EdgeInfo* einfo = iter2->second;
//			int id1 = iter1->first;
//			int id2 = iter2->first;
//
//			if (!einfo->bdry && einfo->toCompMesh){
//
//				if (einfo->pNeighbor_1){
//					// elm does not belong to the computational mesh. It has been replaced by its two children
//					elm = &pBase_EDS[einfo->pNeighbor_1[0]].pLeaves[einfo->pNeighbor_1[1]][einfo->pNeighbor_1[2]];
//
//					// if element is flagged as green it means the edge's neighbor element must point to one of it children
//					if (elm->green){
//						//cout << "aresta: " << id1 << "  " << id2 << ":: ";
//						// elm->green_id1 and elm->green_id2 represent the IDs of one of the green element which is not the
//						// say to edge to which parent's child it should point.
//						// say to edge it points to an element one level above
//						einfo->pNeighbor_1[1] = einfo->pNeighbor_1[1] + 1;
//						einfo->pNeighbor_1[2] = (id1==elm->green_id1 && id2==elm->green_id2)?elm->childLeaf_0:elm->childLeaf_1;
////						cout << "childLeaf_0: " << elm->childLeaf_0 << ", childLeaf_1: " << elm->childLeaf_1 << ", ";
////						cout << "pNeighbor_1: " << " " << einfo->pNeighbor_1[0] << " "
////								                       << einfo->pNeighbor_1[1] << " "
////													   << einfo->pNeighbor_1[2] << "\t ";
////						cout << endl;
////						printElement(elm);
//					}
//				}
//				else{
//					cout << "Error!\n Edge "  << id1 << " "  << id2 << " as pNeighbor_1 NULL \n";
//					exit(1);
//				}
//
//				if (einfo->pNeighbor_2){
//					elm = &pBase_EDS[einfo->pNeighbor_2[0]].pLeaves[einfo->pNeighbor_2[1]][einfo->pNeighbor_2[2]];
//					if (elm->green){
//						cout << "aresta: " << id1 << "  " << id2 << ":: ";
//						einfo->pNeighbor_2[1] = einfo->pNeighbor_2[1] + 1;
//						einfo->pNeighbor_2[2] = (id1==elm->green_id1 && id2==elm->green_id2)?elm->childLeaf_0:elm->childLeaf_1;
//						cout << "childLeaf_0: " << elm->childLeaf_0 << ", childLeaf_1: " << elm->childLeaf_1 << ", ";
//						cout << "pNeighbor_2: " << " " << einfo->pNeighbor_2[0] << " "
//								                       << einfo->pNeighbor_2[1] << " "
//													   << einfo->pNeighbor_2[2];
//						cout << endl;
//						printElement(elm);
//					}
//				}
//				else{
//					cout << "Error!\n Edge "  << id1 << " "  << id2 << " as pNeighbor_2 NULL"; exit(1);
//				}
//			}
//		}
//	}

#ifdef __DEBUG_STEPS__
	cout << "update_greenedge_neighbors(): END\n";
#endif
}

void ParMeshAdapt::create_hangnode_edge(int id1, int id2, TEDGEDB &edgeDB, int* index){
//	if (id1 > id2){
//		swap(id1,id2);
//	}
//
//	EdgeInfo* einfo = new EdgeInfo;
//
//	einfo->geom = -1;
//	einfo->physical = -1;
//	einfo->bdry = false;
//	einfo->green = false;
//	einfo->toCompMesh = true;
//	einfo->markedToSplit = false;
//	einfo->numFaces = 2;
//
//	einfo->pNeighbor_1 = new int[3];
//	einfo->pNeighbor_1[0] = index[0];
//	einfo->pNeighbor_1[1] = index[1]+1;
//	einfo->pNeighbor_1[2] = 4*index[2];
//	einfo->pNeighbor_2 = new int[3];
//	einfo->pNeighbor_2[0] = index[0];
//	einfo->pNeighbor_2[1] = index[1]+1;
//	einfo->pNeighbor_2[2] = 4*index[2]+1;
//
//	bool found0, found1;
//	TEDGEDB::iterator iter;
//	findEdge(id1,id2,found0,found1,iter,edgeDB);
//	std::map<int,EdgeInfo*> secondPart;
//	// first if: id0 already exist into edge DB, but id1.
//	// second if: both id0 and id1 do not exist.
//	if (found0 && !found1){
//		iter->second.insert( std::pair<int,EdgeInfo*>(id2,einfo) );
//	}
//	else if (!found0 && !found1){
//		secondPart.insert( std::pair<int,EdgeInfo*>(id2,einfo) );
//		edgeDB[id1] = secondPart;
//	}

}

void ParMeshAdapt::cleanEdgeDataStructure(){
	for (EdgeIter iter1 = pEdgeDB[0].begin(); iter1 != pEdgeDB[0].end(); iter1++){
		for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
			EdgeInfo* einfo = iter2->second;

			if (einfo->pNeighbor_1){
				delete[] einfo->pNeighbor_1; einfo->pNeighbor_1 =0;
			}
			if (einfo->pNeighbor_2){
				delete[] einfo->pNeighbor_2; einfo->pNeighbor_2 =0;
			}
		}
		iter1->second.clear();
	}
	pEdgeDB[0].clear();
}
