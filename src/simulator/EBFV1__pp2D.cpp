/*
 * EBFV1__pp2D.cpp
 *
 *  Created on: 28 de abr de 2016
 *      Author: rogerio
 */

#include "EBFV1__pre-processors.h"

void pp_2D(ParMeshAdapt* pMesh, GeomData *pGCData){
	debug_msg("Pre-processor 2-D: START");

	double t1 = MPI_Wtime();

	// initialize volumes of control volumes per domain
	init_volumes(pMesh,pGCData->ndom);

	// calculate triangles centroids
	calculate_triagles_centroids_and_areas(pMesh);

	int physical;
	double mp[2], v[2], IJ[2], inner_prod;
	VertexInfo *I, *J;
	Element *elm1, *elm2;
	DataCoefficients* pCoeff; pCoeff = 0;

	// calculate Cij and volume of control volumes per domain
	for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
		for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){

			int id0 = iter1->first;
			int id1 = iter2->first;
			EdgeInfo* edge = iter2->second;
			edge->pEdgeStuffes = 0;

			initCoefficients_2D(&edge);
			DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
			pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
			pMesh->getVertex(id0,&I);
			pMesh->getVertex(id1,&J);

			elm1 = &pMesh->pBase_EDS[edge->pNeighbor_1[0]].pLeaves[edge->pNeighbor_1[1]][edge->pNeighbor_1[2]];
			physical = pMesh->pBase_EDS[edge->pNeighbor_1[0]].physical;

			// to which domain Cij must be associated
			pCoeff->dom1 = pGCData->map_domain_flags[physical];
			pCoeff->dom2 = pCoeff->dom1;

			for (int i=0; i<2; i++){
				mp[i] = 0.5*(I->coords[i]+J->coords[i]);
				IJ[i] = J->coords[i]-I->coords[i];
				v[i] = mp[i] - elm1->pCentroid[i];
			}

			inner_prod = v[1]*IJ[0] + (-v[0])*IJ[1];
			if ( inner_prod <= .0 ){
				v[0] = -v[0];
				v[1] = -v[1];
			}
			pCoeff->Cij[0] += v[1];
			pCoeff->Cij[1] += -v[0];

			if (!edge->bdry){
				elm2 = &pMesh->pBase_EDS[edge->pNeighbor_2[0]].pLeaves[edge->pNeighbor_2[1]][edge->pNeighbor_2[2]];

				// mp and IJ were already computed
				for (int i=0; i<2; i++){
					v[i] = mp[i] - elm2->pCentroid[i];
				}

				inner_prod = v[1]*IJ[0] + (-v[0])*IJ[1];
				if ( inner_prod <= .0 ){
					v[0] = -v[0];
					v[1] = -v[1];
				}
				pCoeff->Cij[0] += v[1];
				pCoeff->Cij[1] += -v[0];
			}
			else{

				pCoeff->Dij[0] = -0.5*(J->coords[1]-I->coords[1]);
				pCoeff->Dij[1] =  0.5*(J->coords[0]-I->coords[0]);
				inner_prod = pCoeff->Dij[0]*(mp[0] - elm1->pCentroid[0]) + pCoeff->Dij[1]*(mp[1] - elm1->pCentroid[1]);
				if (  inner_prod <= .0 ){
					pCoeff->Dij[0] = -pCoeff->Dij[0];
					pCoeff->Dij[1] = -pCoeff->Dij[1];
				}

				if (edge->numFaces==2){
					pCoeff->Dij[2] = -pCoeff->Dij[0];
					pCoeff->Dij[3] = -pCoeff->Dij[1];

					physical = pMesh->pBase_EDS[edge->pNeighbor_2[0]].physical;
					pCoeff->dom2 = pGCData->map_domain_flags[physical];
					elm2 = &pMesh->pBase_EDS[edge->pNeighbor_2[0]].pLeaves[edge->pNeighbor_2[1]][edge->pNeighbor_2[2]];

					// mp and IJ were already computed
					for (int i=0; i<2; i++){
						v[i] = mp[i] - elm2->pCentroid[i];
					}

					inner_prod = v[1]*IJ[0] + (-v[0])*IJ[1];
					if ( inner_prod <= .0 ){
						v[0] = -v[0];
						v[1] = -v[1];
					}

					pCoeff->Cij[2] += v[1];
					pCoeff->Cij[3] += -v[0];
				//	cout << physical << "   "  << pCoeff->Cij[2] <<  " " << pCoeff->Cij[3] << endl;
					elm2 = 0;
				}
			}
			elm1 = 0;
		}
	}

	debug_msg("Pre-processor 2-D: END");
	double t2 = MPI_Wtime();
	cout << "Pre-processor time: " << t2 - t1 << " seconds\n\n";

//	print_coefficients(pMesh);
	validate_coefficients(pMesh,pGCData);
	//exit(1);
}

