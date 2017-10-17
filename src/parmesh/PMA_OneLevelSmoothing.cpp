/*
 * PMA_OneLevelSmoothing.cpp
 *
 *  Created on: 8 de abr de 2016
 *      Author: rogerio
 */

#include "ParMeshAdapt.h"

void ParMeshAdapt::OneLevelSmoothing(){
	OneLevelSmoothing_part1();
	OneLevelSmoothing_part2();
}

void ParMeshAdapt::OneLevelSmoothing_part1(){

#ifdef __DEBUG_STEPS__
	cout << "OneLevelSmoothing_part1(): START\n";
#endif

	int count;
	do{
		count = 0;
		for (EdgeIter iter1 = this->pEdgeDB[0].begin(); iter1 != this->pEdgeDB[0].end(); iter1++){
			for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* einfo = iter2->second;

				// check only edges on the computational mesh and sharing two elements
				// external boundary edges are excluded, but edge between to domains are check as well.
				if (einfo->toCompMesh && einfo->numFaces==2){

//					cout << "\n\n" <<  iter1->first << " " << iter2->first << endl;
//					if (!einfo->pNeighbor_1){
//						cout << "einfo->pNeighbor_1 NULL: " << einfo->pNeighbor_1 << endl;
//						exit(1);
//					}
//					if (!einfo->pNeighbor_2){
//						cout << "einfo->pNeighbor_2 NULL: " << einfo->pNeighbor_2 << endl;
//						exit(1);
//					}

					Element *elm1 = &this->pBase_EDS[einfo->pNeighbor_1[0]].pLeaves[einfo->pNeighbor_1[1]][einfo->pNeighbor_1[2]];
					Element *elm2 = &this->pBase_EDS[einfo->pNeighbor_2[0]].pLeaves[einfo->pNeighbor_2[1]][einfo->pNeighbor_2[2]];

//					if (!elm1){
//						cout << "elm1 NULL: " << elm1 << endl;
//						exit(1);
//					}
//
//					if (!elm2){
//						cout << "elm2 NULL: " << elm2 << endl;
//						exit(1);
//					}
//
//					cout << "green: " << einfo->green << "\n";
//					cout << "N1: " << einfo->pNeighbor_1[0] << " " << einfo->pNeighbor_1[1] << " " << einfo->pNeighbor_1[2] << "\n";
//					cout << "N2: " <<einfo->pNeighbor_2[0] << " " << einfo->pNeighbor_2[1] << " " << einfo->pNeighbor_2[2] << "\n";
//					cout << "ref:" << elm1->numrefine << ", " << elm2->numrefine << endl;
//					printElement(elm1);
//					printElement(elm2);

					//}
					// 2 --3 = 5
					int diff = elm1->numrefine - elm2->numrefine;
					if ( diff > 1 || diff < -1){
						if (elm2->numrefine < elm1->numrefine){
							elm2->numrefine = elm1->numrefine - 1;
							elm2->toRefine = true;
						}
						else{
							elm1->numrefine = elm2->numrefine - 1;
							elm1->toRefine = true;
						}
						count++;
					}
				}
			}
		}
	}while(count);

#ifdef __DEBUG_STEPS__
	cout << "OneLevelSmoothing_part1(): END\n";
#endif
}

void ParMeshAdapt::OneLevelSmoothing_part2(){
#ifdef __DEBUG_STEPS__
	cout << "OneLevelSmoothing_part2(): START\n";
#endif
	int count;
	do{
		count = 0;
		//		for (int i=0; i<=3; i++){
		for (EdgeIter iter1 = this->pEdgeDB[0].begin(); iter1 != this->pEdgeDB[0].end(); iter1++){
			for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* einfo = iter2->second;

				// check only edges on the computational mesh and sharing two elements
				// external boundary edges are excluded, but edge between to domains are check as well.
				if (einfo->toCompMesh && einfo->numFaces==2){
					Element *elm1 = &this->pBase_EDS[einfo->pNeighbor_1[0]].pLeaves[einfo->pNeighbor_1[1]][einfo->pNeighbor_1[2]];
					Element *elm2 = &this->pBase_EDS[einfo->pNeighbor_2[0]].pLeaves[einfo->pNeighbor_2[1]][einfo->pNeighbor_2[2]];

					if ( !elm2->numrefine && elm1->numrefine==1){
						elm2->count++;
					}

					if ( !elm1->numrefine && elm2->numrefine==1){
						elm1->count++;
					}
				}
			}
		}
		//		}

		for (IdxIter it = idxlist.begin(); it!=idxlist.end(); it++){
			IDX_str* idx = &(*it);
			Element *elm = &pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

			if (!elm->numrefine && elm->count==2){
				elm->toRefine = true;
				elm->numrefine = 1;
				count++;
			}
			elm->count = 0;
		}
		cout << "OneLevelSmoothing_part2: count = " << count << endl;
	}while(count);
#ifdef __DEBUG_STEPS__
	cout << "OneLevelSmoothing_part2(): END\n";
#endif
}
