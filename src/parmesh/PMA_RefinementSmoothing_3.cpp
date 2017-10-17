/*
 * PMA_RefinementSmoothing3.cpp
 *
 *  Created on: 3 de abr de 2016
 *      Author: rogerio
 */

#include "ParMeshAdapt.h"

void ParMeshAdapt::RefinementSmoothing_part3(){

	int base1, level1, leaf1;
	int base2, level2, leaf2;
	Element *elm2, *elm1, *elm;
	//int count = 0;

	/*
	 * Identify which elements are surrounded by two other neighbor elements refined one level above.
	 * These elements are not identified during RefinementSmoothing steps (1 and 2) because they don't exist.
	 * ------------------------------------------------------------------------------------------------------------
	 */

	for (EdgeIter iter1 = this->pEdgeDB[0].begin(); iter1 != this->pEdgeDB[0].end(); iter1++){
		for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
			EdgeInfo* einfo = iter2->second;
//			int id1 = iter1->first;
//			int id2 = iter2->first;

			// get element: 1
			base1 = einfo->pNeighbor_1[0];
			level1 = einfo->pNeighbor_1[1];
			leaf1 = einfo->pNeighbor_1[2];
			elm1 = &this->pBase_EDS[base1].pLeaves[level1][leaf1];

			if (!einfo->bdry && einfo->pNeighbor_2){

				// get element: 2
				base2 = einfo->pNeighbor_2[0];
				level2 = einfo->pNeighbor_2[1];
				leaf2 = einfo->pNeighbor_2[2];

				elm2 = &this->pBase_EDS[base2].pLeaves[level2][leaf2];
				if ( !elm2->numrefine ){
					if ( elm1->levelAboveMe ){
						elm2->count++;
					}
				}

				if ( !elm1->numrefine ){
					if ( elm2->levelAboveMe ){
						elm1->count++;
					}
				}
			}
		}
	}

	/*
	 * Refine elements identified at previous loop as surrounded by higher refined elements.
	 * ------------------------------------------------------------------------------------------------------------
	 */

	// before split elements, let me know what is the largest vertex ID.
	int max_vertex_ID;
	VIter vit;
	for(vit = VertexDB.begin(); vit!=VertexDB.end(); vit++){
		max_vertex_ID = vit->first;
	}

	int numLeaves[4] = {1,4,16,64};
	std::list<IDX_str> tmp_idxlist;	// store elements which were not marked to be refined (they belong to the computational mesh)
	for (IdxIter it = idxlist.begin(); it!=idxlist.end(); it++){

		IDX_str* idx = &(*it);
		elm = &pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		if (elm->count==2){
			// allocate memory for leaves to level 'level+1' if necessary
			if (!pBase_EDS[idx->base].pLeaves[idx->level+1]){
				pBase_EDS[idx->base].pLeaves[idx->level+1] = new Element[numLeaves[idx->level+1]];
			}
			splitElement(elm,idx->base,idx->level,idx->leaf,max_vertex_ID);
			elm->levelAboveMe = true;

			// Add leaves to tmp_idxlist_1. They will be analyzed next.
			for (int j=0; j<4; j++){
				IDX_str idx_leaf;
				idx_leaf.base = idx->base;
				idx_leaf.level = idx->level+1;
				idx_leaf.leaf = 4*idx->leaf + j;
				tmp_idxlist.push_back(idx_leaf);
			}
			elm->toCompMesh = false;
			elm->toRefine = false;
			elm->numrefine = 0;
			elm->count = 0;
		}
		else{
			tmp_idxlist.push_back(*it);
		}
	}

	idxlist.clear();
	idxlist.assign(tmp_idxlist.begin(),tmp_idxlist.end());
	tmp_idxlist.clear();

	for (int i=0; i<=3; i++){
		createEdgeDataStructure(i);
	}
}
