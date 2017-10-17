/*
 * validate_coefficients.cpp
 *
 *  Created on: 1 de jul de 2016
 *      Author: rogerio
 */

#include "EBFV1__pre-processors.h"

void validate_coefficients(ParMeshAdapt* pMesh, GeomData* pGCData){
	debug_msg("validate_coefficients 2-D: START\n");

	// Soma dos Cij e Dij da superficie de controle de um nó deve ser nula. Isso vale para qualquer nó da malha
	double t1 = MPI_Wtime();

	// initialize
	int n;
	pMesh->getNumVertices(n);
	n += 2;
	double** sum_mat = new double*[n];
	for (int i=0; i<n; i++){
		sum_mat[i] = new double[2];
	}

	VertexInfo* I; I = 0;
	VertexInfo* J; J = 0;
	EdgeInfo* edge; edge=0;
	DataCoefficients* pCoeff; pCoeff = 0;

	cout << setprecision(18) << scientific;
	int dom, id0, id1, physical1, physical2;
	int ndom = pGCData->getNumDomains();
	for (dom=0; dom<ndom; dom++){

		for (int i=0; i<n; i++){
			sum_mat[i][0] = 0;
			sum_mat[i][1] = 0;
		}

		cout << "Checking domain: " << dom + 1 << " of " << ndom << ":";
		for (EdgeIter iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){

				id0 = iter1->first - 1;
				id1 = iter2->first - 1;

				edge = iter2->second;
				pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
				physical1 = pMesh->pBase_EDS[edge->pNeighbor_1[0]].physical;

				if (dom==pGCData->map_domain_flags[physical1]){

					// Cij points from the small to the high edge vertex ID. id0 < id1 ALWAYS!!!
					sum_mat[id0][0] += pCoeff->Cij[0];
					sum_mat[id0][1] += pCoeff->Cij[1];
					sum_mat[id1][0] += -pCoeff->Cij[0];
					sum_mat[id1][1] += -pCoeff->Cij[1];

					// if edge is on boundary, include Dij
					if (edge->bdry){
						sum_mat[id0][0] += pCoeff->Dij[0];
						sum_mat[id0][1] += pCoeff->Dij[1];
						sum_mat[id1][0] += pCoeff->Dij[0];
						sum_mat[id1][1] += pCoeff->Dij[1];
					}

				}

				if (edge->numFaces==2){
					physical2 = pMesh->pBase_EDS[edge->pNeighbor_2[0]].physical;
					if (physical1 != physical2){
						if (dom==pGCData->map_domain_flags[physical2]){

							// Cij points from the small to the high edge vertex ID. id0 < id1 ALWAYS!!!
							sum_mat[id0][0] += pCoeff->Cij[2];
							sum_mat[id0][1] += pCoeff->Cij[3];
							sum_mat[id1][0] += -pCoeff->Cij[2];
							sum_mat[id1][1] += -pCoeff->Cij[3];

							// if edge is on boundary, include Dij
							if (edge->bdry){
								sum_mat[id0][0] += pCoeff->Dij[2];
								sum_mat[id0][1] += pCoeff->Dij[3];

								sum_mat[id1][0] += pCoeff->Dij[2];
								sum_mat[id1][1] += pCoeff->Dij[3];
							}

						}
					}
				}
			}
		}


		// check summation
		double total_sum[3] = {.0,.0,.0};
		for (int i=0; i<n; i++){
			total_sum[0] += sum_mat[i][0];
			total_sum[0] += sum_mat[i][1];
		}

		if ( fabs(total_sum[0])<1e-15 && fabs(total_sum[1])<1e-15 && fabs(total_sum[2])<1e-15){
			cout << "\nCoefficient validation: PASSED!\n";
		}
		else{
			cerr << "\nCoefficient validation: FAILED!\n";
			exit(1);
		}
		cout << "Total Summation       : " << total_sum[0] << " " << total_sum[1] << " " << total_sum[2] << endl;
		cout << endl;
	}

	for (int i=0; i<n; i++){
		delete[] sum_mat[i]; sum_mat[i] = 0;
	}
	sum_mat = 0;


	/*
	 * Soma dos volumes dos volumes de controle deve ser igual ao volume total do dominio
	 */

	double* p;
	map_int_pDouble volume_sum;
	for(VIter vit = pMesh->VertexDB.begin(); vit!=pMesh->VertexDB.end(); vit++){
		int id = vit->first;
		p = new double[ndom];
		for(int i=0; i<ndom; i++){
			p[i] = .0;
		}
		volume_sum[id] = p;
		p = 0;
	}


	double vol;
	int index, physical;
	VertexInfo *vertex;
	map_int_double volume;
	for(VIter vit = pMesh->VertexDB.begin(); vit!=pMesh->VertexDB.end(); vit++){
		int id = vit->first;
		vertex = vit->second;
		volume = vertex->volume;
		for(map_int_double_Iter it = volume.begin(); it!=volume.end(); it++){
			physical = it->first;
			vol = it->second;
			index = pGCData->map_domain_flags[physical];
			p = volume_sum[id];
			p[index] += vol;
			p = 0;
		}
	}

	//check volume summation:
	double total_volume_sum[ndom];
	for(int i=0; i<ndom; i++){
		total_volume_sum[i] = .0;
	}

	map_int_pDouble_Iter vid_Iter;
	for(vid_Iter = volume_sum.begin(); vid_Iter!=volume_sum.end(); vid_Iter++){
		for(int i=0; i<ndom; i++){
			p = vid_Iter->second;
			total_volume_sum[i] += p[i];

		}
		delete[] p; p = 0;
	}
	volume_sum.clear();

	cout << setprecision(8) << scientific;
	double total = 0;
	cout << "Volume of control volume summation per domain:\n";
	for(int i=0; i<ndom; i++){
		cout << "Domain[" << i << "]   : " << total_volume_sum[i] << endl;
		total += total_volume_sum[i];
	}
	cout << "Total volume: " << total << "\n\n";


	double t2 = MPI_Wtime();
	cout << "Validation time: " << t2 - t1 << endl;

	debug_msg("validate_coefficients 2-D: END");
}


void validate_3D(ParMeshAdapt* pMesh, GeomData* pGCData){
	debug_msg("validate_coefficients 3-D: START\n");

	// Soma dos Cij e Dij da superficie de controle de um nó deve ser nula. Isso vale para qualquer nó da malha
	// ------------------------------------------------------------------------------------------------------------------------
	double t1 = MPI_Wtime();

	// initialize
	int n;
	pMesh->getNumVertices(n);
	n += 2;
	double** sum_mat = new double*[n];
	for (int i=0; i<n; i++){
		sum_mat[i] = new double[3];
	}

	VertexInfo* I; I = 0;
	VertexInfo* J; J = 0;
	EdgeInfo* edge; edge=0;
	ElementFromFile face;
	DataCoefficients* pCoeff; pCoeff = 0;

	cout << setprecision(18) << scientific;
	int i, j, dom, id0, id1, id2, row;
	int ndom = pGCData->getNumDomains();


	/*
	 * Soma dos volumes dos volumes de controle deve ser igual ao volume total do dominio
	 */

	double* p;
	map_int_pDouble volume_sum;
	for(VIter vit = pMesh->VertexDB.begin(); vit!=pMesh->VertexDB.end(); vit++){
		int id = vit->first;
		p = new double[ndom];
		for(int i=0; i<ndom; i++){
			p[i] = .0;
		}
		volume_sum[id] = p;
		p = 0;
	}


	double vol;
	int index, physical;
	VertexInfo *vertex;
	map_int_double volume;
	for(VIter vit = pMesh->VertexDB.begin(); vit!=pMesh->VertexDB.end(); vit++){
		int id = vit->first;
		vertex = vit->second;
		volume = vertex->volume;
		for(map_int_double_Iter it = volume.begin(); it!=volume.end(); it++){
			physical = it->first;
			vol = it->second;
			//cout << "volume: " << vol << endl;
			index = pGCData->map_domain_flags[physical];
			p = volume_sum[id];
			p[index] += vol;
			p = 0;
		}
	}

	//check volume summation:
	double total_volume_sum[ndom];
	for(int i=0; i<ndom; i++){
		total_volume_sum[i] = .0;
	}

	map_int_pDouble_Iter vid_Iter;
	for(vid_Iter = volume_sum.begin(); vid_Iter!=volume_sum.end(); vid_Iter++){
		for(int i=0; i<ndom; i++){
			p = vid_Iter->second;
			total_volume_sum[i] += p[i];

		}
		delete[] p; p = 0;
	}
	volume_sum.clear();

	cout << setprecision(8) << scientific;
	double total = 0;
	cout << "Volume of control volume summation per domain:\n";
	for(int i=0; i<ndom; i++){
		cout << "Domain[" << i << "]   : " << total_volume_sum[i] << endl;
		total += total_volume_sum[i];
	}
	cout << "Total volume: " << total << "\n\n";



	// perform coefficients summation sub-domains by sub-domain
	// ------------------------------------------------------------------------------------------------------------------------
	for (dom=0; dom<ndom; dom++){

		// initialize matrices of summations
		for (i=0; i<n; i++){
			for (j=0; j<3; j++){
				sum_mat[i][j] = 0;
			}
		}

		// loop over all edges of sub-domain 'dom'
		// ------------------------------------------------------------------------------------------------------------------------
		cout << "Checking domain: " << dom + 1 << " of " << ndom << ":";
//		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
//			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
//
//				id0 = iter1->first-1;
//				id1 = iter2->first-1;
//				edge = iter2->second;
//				pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
//
//				// loop over all sub-domains sharing edge: take only that equal to 'dom'
//				for (i=0; i<pCoeff->num_subdomains; i++){
//					if (dom==pCoeff->sub_domain_list[i]){
//
//						double* Cij = pCoeff->C_ij[dom+1];
//						// Cij points from the small to the high edge vertex ID. id0 < id1 ALWAYS!!!
//						sum_mat[id0][0] += Cij[0];
//						sum_mat[id0][1] += Cij[1];
//						sum_mat[id0][2] += Cij[2];
//
//						sum_mat[id1][0] += -Cij[0];
//						sum_mat[id1][1] += -Cij[1];
//						sum_mat[id1][2] += -Cij[2];
//
//						//cout << "\nsum_mat: " << sum_mat[id0][0] << " " << sum_mat[id0][1] << " " << sum_mat[id0][2] << "  ,  " << sum_mat[id1][0] << " " << sum_mat[id1][1] << " " << sum_mat[id1][2];
//					}
//				}
//			}
//		}


		// loop over all faces belonging to sub-domain 'dom'
		// ------------------------------------------------------------------------------------------------------------------------
		row = 0;
		double *Dij; Dij=0;

		list_ElmFromFile *face_list; face_list=0;
		pMesh->getBdryFace_list(&face_list);
		for(auto it=face_list->begin(); it!=face_list->end(); it++){
			face = *it;
			auto it1 = pGCData->surface_map.find(face.geom);
			for (auto it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
				if ( (*it2-1)==dom){

					pGCData->getDij(dom,row,&Dij);
					cout << "\ndom:" << dom << "  row:" << row << "\tFace: " << face.id1 << " " << face.id2 << " " << face.id3 << " Dij: " << Dij[0] << " " << Dij[1] << " " << Dij[2] << "";
					row++;


					id0 = face.id1 - 1;
					id1 = face.id2 - 1;
					id2 = face.id3 - 1;

					for(i=0; i<3; i++){
						sum_mat[id0][i] += Dij[i];
						sum_mat[id1][i] += Dij[i];
						sum_mat[id2][i] += Dij[i];
					}
					Dij = 0;
				}
			}
		}
		face_list=0;

		//		for (i=0; i<n; i++){
		//			cout << "Face[" << id0 << "," << id1 << "," << id2 << "]\t" << "sum_mat: " << sum_mat[i][0] << " " << sum_mat[i][1] << " " << sum_mat[i][2] << endl;
		//		}

		// check summation of Cij and Dij
		// ------------------------------------------------------------------------------------------------------------------------
		double total_sum[3] = {.0,.0,.0};
		for (int i=0; i<n; i++){
			total_sum[0] += sum_mat[i][0];
			total_sum[1] += sum_mat[i][1];
			total_sum[2] += sum_mat[i][2];
		}

		if ( fabs(total_sum[0])<1e-15 && fabs(total_sum[1])<1e-15 && fabs(total_sum[2])<1e-15){
			cout << "\nCoefficient validation: PASSED!\n";
		}
		else{
			cout << "\nTotal Summation       : " << total_sum[0] << " " << total_sum[1] << " " << total_sum[2] << endl;
			//throw Exception(__LINE__,__FILE__,"\nCoefficient validation: FAILED!\n");
			cerr << "Coefficient validation: FAILED!\n";
			cerr << "Remember that face's connectivities must be listed in counterclockwise!\n";
			exit(1);
		}
		cout << "Total Summation       : " << total_sum[0] << " " << total_sum[1] << " " << total_sum[2] << endl;
		cout << endl;
	}

	for (int i=0; i<n; i++){
		delete[] sum_mat[i]; sum_mat[i] = 0;
	}
	sum_mat = 0;





	double t2 = MPI_Wtime();
	cout << "Validation time: " << t2 - t1 << endl;
	//exit(1);
	debug_msg("validate_coefficients 3-D: END");
}

