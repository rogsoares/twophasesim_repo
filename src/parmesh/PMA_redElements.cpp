/*
 * PMA_split.cpp
 *
 *  Created on: 29 de fev de 2016
 *      Author: rogerio
 */


#include "ParMeshAdapt.h"

void ParMeshAdapt::refine(){

//	double t1 = MPI_Wtime();

	int numLeaves[4] = {1,4,16,64};

	// before split elements, let me know what is the largest vertex ID.
	int max_vertex_ID;
	VIter vit;
	for(vit = VertexDB.begin(); vit!=VertexDB.end(); vit++){
		max_vertex_ID = vit->first;
	}
	//cout << "Max vertex ID: " << max_vertex_ID << endl;

	// start element splitting process
	IDX_str* idx; idx=0;


	// -----------------------------------------------------------------------------------------------
	// loop through elements of the computational mesh of size 'size'
	// As element leaves are being created, they are added to the end of idxlist, so idxlist increases
	// varying its size. We must go through only over the first 'size' elements of the list
	// ------------------------------------------------------------------------------------------------


	std::list<IDX_str> tmp_idxlist_1;	// store elements which were not marked to be refined (they belong to the computational mesh)
	std::list<IDX_str> tmp_idxlist_2;	// store elements recently created from level 'level' (so, they belong to level 'level+1' and
										// may be required to be split to level 'level+2'

	int max_level = 1;

	// loop over all computational elements mesh level by level, i.e, from lowest (level=0)
	// to highest (level=max_level) searching for elements marked to be refined.
	// ---------------------------------------------------------------------------------------------------------
	for (int i=0; i<=3; i++){
		int count = 0;
		IdxIter it = idxlist.begin();
		while ( it != idxlist.end() ){
			idx = &(*it);
			Element* elm = &this->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
			if (idx->level == i){
				if (elm->toRefine){
					count++;

					// update max_level
					max_level = std::max(elm->numrefine,max_level);

					// allocate memory for leaves to level 'level+1' if necessary
					if (!pBase_EDS[idx->base].pLeaves[idx->level+1]){
						pBase_EDS[idx->base].pLeaves[idx->level+1] = new Element[numLeaves[idx->level+1]];
						//cout << "Parent at level["<<idx->level<<"] has " << numLeaves[idx->level+1] << " children\n" << endl;
					}
					splitElement(elm,idx->base,idx->level,idx->leaf,max_vertex_ID);
					elm->levelAboveMe = true;

					// Add leaves to tmp_idxlist_2. They will be analyzed next.
					for (int j=0; j<4; j++){

						IDX_str idx_leaf;
						idx_leaf.base = idx->base;
						idx_leaf.level = idx->level+1;
						idx_leaf.leaf = 4*idx->leaf + j;
						tmp_idxlist_2.push_back(idx_leaf);
					}
					elm->toCompMesh = false;
					elm->toRefine = false;
					elm->numrefine = 0;
				}
				else{
					// at the end, this list will hold all elements of the computational mesh
					tmp_idxlist_1.push_back(*idx);
					elm->levelAboveMe = false;
				}
			}
			it++;
		}

		// let's update idxlist. First, clean it up.
		if (count){
			if ( tmp_idxlist_2.size() ){
				idxlist.clear();
				// assign tmp_idxlist_2 elements to idxlist
				idxlist.assign(tmp_idxlist_2.begin(),tmp_idxlist_2.end());
				// clean tmp_idxlist_2
				tmp_idxlist_2.clear();
				if (i+1<=3){
					createEdgeDataStructure(i+1);
				}
			}
		}
	}

	idxlist.clear();
	idxlist.assign(tmp_idxlist_1.begin(),tmp_idxlist_1.end());
	tmp_idxlist_1.clear();

	//double t2 = MPI_Wtime();
}


/*
 * splitElement always split an element into four new elements. There is any special treatment for any element.
 */
void ParMeshAdapt::splitElement(Element* pElm, int base, int level, int leaf_parent, int &max_vertex_ID){

	if (!pEdgeDB[0].size()){
		cerr << "ERROR: Edge data structure has not been created! Exiting...\n";
		exit(1);
	}

	// get triangle edges:
	EdgeInfo* einfo; einfo = 0;
//	getEdge(pElm->id1,pElm->id2,&einfo,pEdgeDB[0]);

	VertexInfo* V1=0;
	VertexInfo* V2=0;
	VertexInfo* V_new=0;

	if (!einfo->markedToSplit){
		++max_vertex_ID;
		getVertex(pElm->id1,&V1);
		getVertex(pElm->id2,&V2);

		V_new = new VertexInfo;
		V_new->coords = new Coords[3];
		middlePoint(V1->coords,V2->coords,V_new->coords);
		createVertex(max_vertex_ID,V_new);
		einfo->MPV_id = max_vertex_ID;
		einfo->markedToSplit = true;
		einfo->toCompMesh = false;
	}
	int R = einfo->MPV_id;

	//getEdge(pElm->id1,pElm->id3,&einfo,pEdgeDB[0]);
	if (!einfo->markedToSplit){
		++max_vertex_ID;
		getVertex(pElm->id1,&V1);
		getVertex(pElm->id3,&V2);

		V_new = new VertexInfo;
		V_new->coords = new Coords[3];
		middlePoint(V1->coords,V2->coords,V_new->coords);
		createVertex(max_vertex_ID,V_new);
		einfo->MPV_id = max_vertex_ID;
		einfo->markedToSplit = true;
		einfo->toCompMesh = false;
	}
	int T = einfo->MPV_id;

	//getEdge(pElm->id2,pElm->id3,&einfo,pEdgeDB[0]);
	if (!einfo->markedToSplit){
		++max_vertex_ID;
		getVertex(pElm->id2,&V1);
		getVertex(pElm->id3,&V2);
		V_new = new VertexInfo;
		V_new->coords = new Coords[3];
		middlePoint(V1->coords,V2->coords,V_new->coords);
		createVertex(max_vertex_ID,V_new);
		einfo->MPV_id = max_vertex_ID;
		einfo->markedToSplit = true;
		einfo->toCompMesh = false;
	}
	int S = einfo->MPV_id;

	int leaf = 4*leaf_parent;

	// create element for the next level: pElm->id1-R-T, S-R-T, R-S-pElm->id2, S-T-pElm->id3
	Element *elm_1 = &this->pBase_EDS[base].pLeaves[level+1][leaf];
	Element *elm_2 = &this->pBase_EDS[base].pLeaves[level+1][leaf+1];
	Element *elm_3 = &this->pBase_EDS[base].pLeaves[level+1][leaf+2];
	Element *elm_4 = &this->pBase_EDS[base].pLeaves[level+1][leaf+3];

	elm_1->id1 = pElm->id1;
	elm_1->id2 = R;
	elm_1->id3 = T;
	elm_1->toUnRefine = false;
	elm_1->toCompMesh = true;
	elm_1->numrefine = pElm->numrefine-1;
	elm_1->toRefine = (elm_1->numrefine > 0) && (elm_1->numrefine <= 3);
	elm_1->green = false;

	elm_2->id1 = S;
	elm_2->id2 = R;
	elm_2->id3 = T;
	elm_2->toUnRefine = false;
	elm_2->toCompMesh = true;
	elm_2->numrefine = pElm->numrefine-1;
	elm_2->toRefine = (elm_2->numrefine > 0) && (elm_2->numrefine <= 3);
	elm_2->green = false;

	elm_3->id1 = R;
	elm_3->id2 = S;
	elm_3->id3 = pElm->id2;
	elm_3->toUnRefine = false;
	elm_3->toCompMesh = true;
	elm_3->numrefine = pElm->numrefine-1;
	elm_3->toRefine = (elm_3->numrefine > 0) && (elm_3->numrefine <= 3);
	elm_3->green = false;

	elm_4->id1 = S;
	elm_4->id2 = T;
	elm_4->id3 = pElm->id3;
	elm_4->toUnRefine = false;
	elm_4->toCompMesh = true;
	elm_4->numrefine = pElm->numrefine-1;
	elm_4->toRefine = (elm_4->numrefine > 0) && (elm_4->numrefine <= 3);
	elm_4->green = false;
}
