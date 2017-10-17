/*
 * PMA_greenElements.cpp
 *
 *  Created on: 31 de mar de 2016
 *      Author: rogerio
 */

#include "ParMeshAdapt.h"

/*
 * Strategy:
 *
 * Loop over all edges from level: 0 to highest-1.
 *
 * For each edge: check if it has been split
 *
 * 		if true:
 *
 * 			look for its neighbor elements Elm1 and Elm2.
 *
 * 			if (Elm1 has leaves above) AND (Elm2 doesn't have leaves above)
 *
 * 				split Elm2 into two new elements
 *
 * 			end
 *
 * 		end
 * end
 *
 *
 */
void ParMeshAdapt::buildGreenElements(){
#ifdef __DEBUG_STEPS__
	cout << "buildGreenElements(): START\n";
#endif

	int max_vertex_ID;
	VIter vit;
	for(vit = VertexDB.begin(); vit!=VertexDB.end(); vit++){
		max_vertex_ID = vit->first;
	}

	std::list<int> greenEdges_List;

	// loop over CM edge DB
	int hang_node, id1,id2,id3;
	Element *elm, *elm1, *elm2;
	int* index;
	int greenEdge_1[2];
	int greenEdge_2[2];
	int count = 0;
	int numLeaves[4] = {1,4,16,64};

	//	for (level = 0; level<=3; level++){
	for (EdgeIter iter1 = pEdgeDB[0].begin(); iter1 != pEdgeDB[0].end(); iter1++){
		for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){

			id1 = iter1->first;
			id2 = iter2->first;

			// get edge
			EdgeInfo* einfo = iter2->second;

			//				if (id1==595 && id2==3267){
			//					cout << "einfo->markedToSplit: " << einfo->markedToSplit << endl;
			//					cout << "einfo->bdry: " << einfo->bdry << endl;
			//					cout << "Sim, eu existo!\n";
			//					cout << __LINE__<<endl;
			//					//exit(1);
			//				}

			if (einfo->markedToSplit){

				hang_node = einfo->MPV_id;

				// get element: 1
				elm1 = &pBase_EDS[einfo->pNeighbor_1[0]].pLeaves[einfo->pNeighbor_1[1]][einfo->pNeighbor_1[2]];

				if (!einfo->bdry){
					// get element: 2
					elm2 = &pBase_EDS[einfo->pNeighbor_2[0]].pLeaves[einfo->pNeighbor_2[1]][einfo->pNeighbor_2[2]];

					elm = 0;
					index = 0;

					// here we suppose elm2 is about to be split into two new elements
					if (elm1->levelAboveMe && !elm2->levelAboveMe){
						elm = elm2;
						index = einfo->pNeighbor_2;
					}

					// here we suppose elm2 is about to be split into two new elements
					if (!elm1->levelAboveMe && elm2->levelAboveMe){
						elm = elm1;
						index = einfo->pNeighbor_1;
					}

					// elm is the green element to be split into two new leaves
					if (elm){
						// before create new elements, check if there is space in memory
						if (!pBase_EDS[index[0]].pLeaves[index[1]+1]){
							pBase_EDS[index[0]].pLeaves[index[1]+1] = new Element[numLeaves[index[1]+1]];
						}
						//
						// precisamos saber entre quais IDs do elemento verde o hang node se localiza
						int greenElmIDs[6];
						int ID[3] = {elm->id1, elm->id2, elm->id3};
						sort(ID,ID+3);

						if (ID[0]==id1 && ID[1]==id2){
							greenElmIDs[0] = ID[0];
							greenElmIDs[1] = ID[2];
							greenElmIDs[2] = hang_node;

							greenElmIDs[3] = ID[1];
							greenElmIDs[4] = ID[2];
							greenElmIDs[5] = hang_node;

							greenEdge_1[0] = ID[0];
							greenEdge_1[1] = ID[2];
							greenEdge_2[0] = ID[1];
							greenEdge_2[1] = ID[2];

							id3 = ID[2];

						}
						else if (ID[1]==id1 && ID[2]==id2){
							greenElmIDs[0] = ID[0];
							greenElmIDs[1] = ID[1];
							greenElmIDs[2] = hang_node;

							greenElmIDs[3] = ID[0];
							greenElmIDs[4] = ID[2];
							greenElmIDs[5] = hang_node;

							greenEdge_1[0] = ID[0];
							greenEdge_1[1] = ID[1];
							greenEdge_2[0] = ID[0];
							greenEdge_2[1] = ID[2];

							id3 = ID[0];
						}
						else if (ID[0]==id1 && ID[2]==id2){
							greenElmIDs[0] = ID[0];
							greenElmIDs[1] = ID[1];
							greenElmIDs[2] = hang_node;

							greenElmIDs[3] = ID[1];
							greenElmIDs[4] = ID[2];
							greenElmIDs[5] = hang_node;

							greenEdge_1[0] = ID[0];
							greenEdge_1[1] = ID[1];
							greenEdge_2[0] = ID[1];
							greenEdge_2[1] = ID[2];

							id3 = ID[1];
						}

						splitGreen(greenElmIDs,index[0],index[1],index[2]);

						/*
						 * Elemento pai: set as green
						 *
						 * A green parent has two children.
						 * An edge pointing to a green parent element could be informed to
						 * not point to this parent anymore but to its children.
						 * Green parent knows who are its children: 0 and 1.
						 * Green parent knows which edge's ID belong to each child.
						 *
						 * These ideas could solve the terrible cpu performance raised
						 * by using "kill_greenEdge_neighbors"and "update_greenedge_neighbors"
						 * functions.
						 */
						elm->levelAboveMe = true;
						elm->toCompMesh = false;
						elm->green = true;
						elm->numrefine = 0;
						elm->green_id1 = id1;
						elm->green_id2 = id3;
						elm->childLeaf_0 = 4*index[2];
						elm->childLeaf_1 = 4*index[2]+1;
//						cout << "aresta: " << id1 << "  " << id2 << ":: ";
//						cout << "childLeaf: " << elm->childLeaf_0 << " " << elm->childLeaf_1 << endl;
						//printElement(elm);
						count++;

						//create_hangnode_edge(id3,hang_node,pEdgeDB[0],index);
					}
				}
			}
		}
	}

	cout << "Step 2:\n";

	count = 0;
	std::list<IDX_str> greenlist;
	std::list<IDX_str> cpmlist;

	cout << "idxlist size: " << idxlist.size() << endl;
	for (IdxIter it = idxlist.begin(); it!=idxlist.end(); it++){
		IDX_str* idx = &(*it);
		elm = &pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		// identify elements on idx list set to not be included into the computational list.
		if (elm->green){
			for (int i=0; i<2; i++){
				IDX_str idx_leaf;
				idx_leaf.base = idx->base;
				idx_leaf.level = idx->level+1;
				idx_leaf.leaf = 4*idx->leaf + i;
				greenlist.push_back(idx_leaf);
			}
			count++;
		}

		if (elm->toCompMesh){
			elm->numrefine = 0;
			cpmlist.push_back(*it);
		}
	}

	idxlist.clear();

	for (IdxIter it = greenlist.begin(); it!=greenlist.end(); it++){
		idxlist.push_back(*it);
	}
	greenlist.clear();

	for (int i=0; i<=3; i++){
		createEdgeDataStructure(i);
	}

	for (IdxIter it = cpmlist.begin(); it!=cpmlist.end(); it++){
		idxlist.push_back(*it);
	}
	cpmlist.clear();

	this->cleanEdgeDataStructure();
	for (int i=0; i<=3; i++){
		createEdgeDataStructure(i);
	}

	// looking for NULL pointers
	// -----------------------------------------------------------------------------
//	for (EdgeIter iter1 = pEdgeDB[0].begin(); iter1 != pEdgeDB[0].end(); iter1++){
//		for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
//			EdgeInfo* einfo = iter2->second;
//			id1 = iter1->first;
//			id2 = iter2->first;
//
//			if (!einfo->bdry && einfo->toCompMesh){
//				if (!einfo->pNeighbor_1 || !einfo->pNeighbor_2){
//					cout << "edge: " << id1 << " " << id2 << " has: ";
//					cout << "N1: " << einfo->pNeighbor_1 << " N2: " << einfo->pNeighbor_2 << "\n";
//					exit(1);
//				}
//				else{
//					cout << "pNeighbor_1: " << einfo->pNeighbor_1[0] << " " << einfo->pNeighbor_1[1] << " " << einfo->pNeighbor_1[2] << "\t";
//					cout << "pNeighbor_2: " << einfo->pNeighbor_2[0] << " " << einfo->pNeighbor_2[1] << " " << einfo->pNeighbor_2[2] << "\n";
//				}
//			}
//
//		}
//	}
//	update_greenedge_neighbors();
//
//	for (int i=0; i<=3; i++){
//		createEdgeDataStructure(i);
//	}

//	exit(1);
#ifdef __DEBUG_STEPS__
	cout << "buildGreenElements(): END\n";
#endif
}

void ParMeshAdapt::splitGreen(int* greenElmIDs,int base, int level, int leaf_parent){

	//cout << "green: ";
	int k = 0;
	int leaf = 4*leaf_parent;
	for (int i=0; i<2; i++){

		//cout << base << " " << level+1 << " " << leaf+i << "  \t";
		Element *elm = &pBase_EDS[base].pLeaves[level+1][leaf+i];
		elm->id1 = greenElmIDs[k++];
		elm->id2 = greenElmIDs[k++];
		elm->id3 = greenElmIDs[k++];
		elm->toUnRefine = false;
		elm->toCompMesh = true;
		elm->numrefine = 0;
		elm->toRefine = 0;
		elm->levelAboveMe = false;
		elm->green = false;
	}
	//cout << endl;
}
