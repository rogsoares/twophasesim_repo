/*
 * PMA_unrefine.cpp
 *
 *  Created on: 17 de mar de 2016
 *      Author: rogerio
 */

#include "ParMeshAdapt.h"


/*
 * PMA_UnrefinementSmoothing_1.cpp
 *
 *  Created on: 7 de abr de 2016
 *      Author: rogerio
 */

#include "PMA_refineTests.h"


/*
 * Verifique quais elementos devem ser removidos.
 * Esta verificacao é feita a partir dos elementos da malha computacional olhando para o pai. Este recebe a média das
 * marcacoes atribuida a seus filhos:
 *
 * 			marcacao
 * 		media        pai   filhos
 *       -3          -2      E
 * 	 -3 < m <=-2     -1      E
 * 	 -2 < m <=-1      0      E
 * 	 -1 < m <=0       1      P
 * 	  0 < m <=1       2      P
 * 	  1 < m <=2       3      P
 *
 *  E: eliminado
 *  P: permanece
 */

void ParMeshAdapt::setElementsTobeRemoved(){
#ifdef __DEBUG_STEPS__
	cout << "setElementsTobeRemoved: START\n";
#endif

	Element *child[4];
	Element *elm, *parent;
	int i, child_leaf, parent_leaf;

	std::list<IDX_str> list;

	for (IdxIter it = idxlist.begin() ;it!=idxlist.end(); it++){
		IDX_str* idx = &(*it);
		if (!idx->level){
			list.push_back(*idx);
		}
	}

	// elimine os elementos de cima (level=3) para baixo (level=1). Level = 0 não são alterados
	for (int level=3; level>0; level--){
		for (IdxIter it = idxlist.begin() ;it!=idxlist.end(); it++){
			IDX_str* idx = &(*it);

			if (idx->level == level){
				elm = &pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

				//cout << "elm->toCompMesh: " << elm->toCompMesh << ",  ";

				// as long as the list is run, some elements are set to be excluded from the computational mesh
				// indirectly by one of its brothers. So, when its turn come, calculations MUST not be repeated.
				if (elm->toCompMesh){
					parent_leaf = IndexMat[0][idx->leaf];
					parent = &pBase_EDS[idx->base].pLeaves[idx->level-1][parent_leaf];

					//					cout << "parent_leaf: " << parent_leaf << ",  ";
					//					cout << "parent->green: " << parent->green << ",  ";



					double mean = .0;
					// loop over all its children and calculate mean value

					int div = (parent->green)?2:4;


					for (i=0; i<div; i++){
						child_leaf = 4*parent_leaf + i;
						child[i] = &pBase_EDS[idx->base].pLeaves[idx->level][child_leaf];
						mean += child[i]->numrefine;
					}
					mean = (double)(mean/div);


					if ( mean<=-2.5 ){
						parent->numrefine = -2;
					}
					else if ( mean>-2.5 && mean<=-1.5 ){
						parent->numrefine = -1;
					}
					if ( mean>-1.5 && mean<=-0.5 ){
						parent->numrefine = 0;
					}

					//cout << "mean: " << mean << endl;

					// if arithmetic mean is less than -0.5, it means that all parent's children must
					// be put out of the computational mesh.
					if (mean <= -0.5){
						for (i=0; i<div; i++){
							child[i]->toCompMesh = false;
							//child[i] = 0;
						}

						IDX_str idx_parent;
						idx_parent.base = idx->base;
						idx_parent.level = idx->level-1;
						idx_parent.leaf = parent_leaf;
						list.push_back(idx_parent);
						parent->toCompMesh = true;
					}
					else{
						list.push_back(*idx);
					}
					//}
				}
			}
		}
	}

	idxlist.clear();
	for (IdxIter it = list.begin(); it!=list.end(); it++){
		idxlist.push_back(*it);
	}
	list.clear();

#ifdef __DEBUG_STEPS__
	cout << "setElementsTobeRemoved: END\n";
#endif
}



//
///*
//
// Unrefining criteria:
//
// 	 	   idx :   0   1   2   3
// 	 	 leaves:  -3  -3  -3  -3
// 	 	 parent:        -2
//
// 	 	 leaves:  -3  -2  -2  -2
// 	 	 parent:        -1
//
// 	 	 leaves:  -3  -2  -2  -2
// 	 	 parent:        -1
//
//
//
//
// */
//
//void ParMeshAdapt::unrefine(){
//	unrefine(3);
//}
//
//void ParMeshAdapt::unrefine(int level){
//
//	double t1 = MPI_Wtime();
//
//
//		Element *parent; parent=0;
//		std::list<IDX_str> tmp_idxlist;
//		std::list<IDX_str> tmp_idxlist1;
//
//
//		for (int i=level; i>0; i--){
//			cout << "Unrefinig for level: " << i << " to " << i-1 << endl;
//			IdxIter it = idxlist.begin();
//			int count=0;
//			while ( it != idxlist.end() ){
//				IDX_str *idx = &(*it);
//
//				// unrefine only elements located at level: 'level'
//				if (idx->level == i){
//
//					Element *elm = &pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
//				//	cout << "i = " << i << " idx->level = "  << idx->level << " elm->toUnRefine = "  << elm->toUnRefine;
//
//					if (elm->toUnRefine){
//						elm->toCompMesh = false;
//
//						// indices:       0                  1                  2                  3
//						// level 1:     elem_0      |     elem_0      |       elem_0      |     elem_0         |
//						// level 2: e_0 e_1 e_2 e_3 | e_4 e_5 e_6 e_7 | e_8 e_9 e_10 e_11 | e_0  e_1  e_2  e_3 |
//						// indices:  0   1   2   3  |  4   5   6   7  |  8   9   10   11  |  12   13   14   15 |
//
//						// Ex.: leaf e_6 -> index=6, parent has index=1
//						// Ex.: leaf e_12 -> index=6, parent has index=3
//
//						// do not include 'parent' back to the computational mesh more than once.
//						// it could happen caused by its leaves
//						parent = &this->pBase_EDS[idx->base].pLeaves[idx->level-1][ this->IndexMat[0][idx->leaf] ];
//
////						cout << "\tleaf " << idx->leaf << "  has parent " << this->IndexMat[0][idx->leaf] <<
////								" parent->toCompMesh: "  << parent->toCompMesh << endl;
//
//						if (!parent->toCompMesh){
//							parent->toCompMesh = true;
//
//							IDX_str idx_parent;
//							idx_parent.base = idx->base;
//							idx_parent.level = idx->level-1;
//							idx_parent.leaf = this->IndexMat[0][idx->leaf];
//
//							parent->toUnRefine = (idx_parent.level>0);
//							tmp_idxlist.push_back(idx_parent);
//						}
//					}
//					else{
//						tmp_idxlist1.push_back(*idx);
//					}
//				}
//				else{
//					// if element does not belong to level 'level' hold it into another list
//					tmp_idxlist.push_back(*idx);
//				}
//				it++;
//			}
//
//			cout << "idxlist: "  << idxlist.size() << endl;
//			cout << "tmp_idxlist: "  << tmp_idxlist.size() << endl;
//			cout << "tmp_idxlist1: "  << tmp_idxlist1.size() << endl;
//			cout << "\n\n\n";
//
//
//			idxlist.clear();
//			idxlist.assign(tmp_idxlist.begin(),tmp_idxlist.end());
//			tmp_idxlist.clear();
//
//
//			cout << "Removing vertices over some edges...";
//			// let's remove some vertices. These vertices are located in the middle of edges of level: 'level-1'
//			EdgeIter it1 = pEdgeDB[i-1].begin();
//			std::map<int, EdgeInfo*>::iterator it2;
//			for(;it1!=pEdgeDB[i-1].end();it1++){
//				for(it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
//					EdgeInfo* edge = it2->second;
//					if (edge->markedToSplit){
//						std::map<int, VertexInfo*>::iterator itmap = VertexDB.find(edge->MPV_id);
//						if (itmap==VertexDB.end()){
//							cout << "Error: vertex " << edge->MPV_id <<  "  not found\n";
//							exit(1);
//						}
//						VertexInfo *v = itmap->second;
//						delete[] v->coords; v->coords = 0;
//						VertexDB.erase(itmap);
//
//						edge->MPV_id = 0;
//					}
//				}
//			}
//			cout << "done!\n";
//		}
//		cout << "Done!\n";
//}
//
///*
//
//void ParMeshAdapt::unrefine(int level){
//
//
//	double t1 = MPI_Wtime();
//	Element *parent; parent=0;
//	std::list<IDX_str> tmp_idxlist;
//	std::list<IDX_str> tmp_idxlist1;
//
//	for (int i=level; i>0; i--){
//		cout << "Unrefinig for level: " << i << " to " << i-1 << endl;
//		IdxIter it = idxlist.begin();
//		int count=0;
//		while ( it != idxlist.end() ){
//			IDX_str *idx = &(*it);
//
//			// unrefine only elements located at level: 'level'
//			if (idx->level == i){
//
//				Element *elm = &pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
//				cout << "i = " << i << " idx->level = "  << idx->level << " elm->toUnRefine = "  << elm->toUnRefine << endl;
//
//				if (elm->toUnRefine){
//					elm->toCompMesh = false;
//
//					// indices:       0                  1                  2                  3
//					// level 1:     elem_0      |     elem_0      |       elem_0      |     elem_0         |
//					// level 2: e_0 e_1 e_2 e_3 | e_4 e_5 e_6 e_7 | e_8 e_9 e_10 e_11 | e_0  e_1  e_2  e_3 |
//					// indices:  0   1   2   3  |  4   5   6   7  |  8   9   10   11  |  12   13   14   15 |
//
//					// Ex.: leaf e_6 -> index=6, parent has index=1
//					// Ex.: leaf e_12 -> index=6, parent has index=3
//
//					// do not include 'parent' back to the computational mesh more than once.
//					// it could happen caused by its leaves
//					parent = &this->pBase_EDS[idx->base].pLeaves[idx->level-1][ this->IndexMat[0][idx->leaf] ];
//					if (!parent->toCompMesh){
//						parent->toCompMesh = true;
//
//						IDX_str idx_parent;
//						idx_parent.base = idx->base;
//						idx_parent.level = idx->level-1;
//						idx_parent.leaf = this->IndexMat[0][idx->leaf];
//
//						parent->toUnRefine = (idx_parent.level>0);
//						tmp_idxlist.push_back(idx_parent);
//					}
//				}
//				else{
//					tmp_idxlist1.push_back(*idx);
//				}
//			}
//			else{
//				// if element does not belong to level 'level' hold it into another list
//				tmp_idxlist.push_back(*idx);
//			}
//			it++;
//		}
//		idxlist.clear();
//		idxlist.assign(tmp_idxlist.begin(),tmp_idxlist.end());
//		tmp_idxlist.clear();
//
//
//		cout << "Removing vertices over some edges...";
//		// let's remove some vertices. These vertices are located in the middle of edges of level: 'level
//		EdgeIter it1 = pEdgeDB[i-1].begin();
//		std::map<int, EdgeInfo*>::iterator it2;
//		for(;it1!=pEdgeDB[i-1].end();it1++){
//			for(it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
//				EdgeInfo* edge = it2->second;
//				if (edge->markedToSplit){
//					std::map<int, VertexInfo*>::iterator it = VertexDB.find(edge->MPV_id);
//					VertexInfo *v = it->second;
//					delete[] v->coords; v->coords = 0;
//					VertexDB.erase(it);
//				}
//			}
//		}
//		cout << "done!\n";
//	}
//
//	idxlist.clear();
//	idxlist.assign(tmp_idxlist1.begin(),tmp_idxlist1.end());
//	tmp_idxlist1.clear();
//	cout << "Done!\n";
//
//
//
//
//	double t2 = MPI_Wtime();
//	cout << "UnRefining end.\tElapsed time: " << t2-t1 << endl;;
//}
//
//*/
//
////
////			/*
////			 get elm's parent (one level below) and allow it (parent)
////			 to be included into the computational mesh again
////
////			 Element to be unrefined ("removed" from computational mesh):
////			 	 	 	 	 leaf_elm = this->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf]
////
////			                   parent
////			                     |
////			       ----------------------------
////			       |        |        |        |
////			     leaf_1   leaf_2   leaf_3   leaf_4
////
////			     'leaf_elm' can be one of the four parent's leaves
////			 */
