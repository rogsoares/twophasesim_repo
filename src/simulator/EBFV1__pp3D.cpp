/*
 * EBFV1__pp3D.cpp
 *
 *  Created on: 28 de abr de 2016
 *      Author: rogerio
 */

#include "EBFV1__pre-processors.h"

void pp_3D(ParMeshAdapt* pMesh, GeomData *pGCData){
	debug_msg("Pre-processor 3-D: START");
	double t1 = MPI_Wtime();

	initCoefficients_3D(pMesh);

	//exit(1);
	init_volumes(pMesh,pGCData->ndom);			// initialize volumes of control volumes per domain
	calculate_volume(pMesh);					// calculate tetrahedral's volumes
	calculate_centroids(pMesh);					// calculate tetrahedral's centroids
	calculate_face_centroids(pMesh);
	calculate_Cij(pMesh);

	calculate_NumFacesPerSubDomain(pMesh,pGCData);
	calculate_Dij(pMesh,pGCData);
//
//	cout << setprecision(8) << fixed;
//	for (EdgeIter iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
//		for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
//			EdgeInfo* edge = iter2->second;
//			DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
//
//			cout << "edge[" << iter1->first << "," << iter2->first << "] ";
//			std::map<int,double*>::iterator it;
//			for(it=pCoeff->C_ij.begin(); it!=pCoeff->C_ij.end(); it++){
//				double *cij = it->second;
//				cout << "Cij[" << it->first << "]: " << cij[0] << " " << cij[1] << " " << cij[2] << "\t";
//			}
//			cout << endl;
//			pCoeff = 0;
//		}
//	}
//	exit(1);

	cleanup_auxiliary_data(pMesh);

	debug_msg("Pre-processor 3-D: END");
	double t2 = MPI_Wtime();
	cout << "Pre-processor ";
	convertSecToTime(t2-t1);
	//
	//print_coefficients(pMesh);
	validate_3D(pMesh,pGCData);
	//exit(1);
}

void compute_Cij(ParMeshAdapt* pMesh, int id0, int id1, int geom, const double* I_coords, const double* J_coords, const double* face1_centroid, const double* face2_centroid, const double* tetra_centroid){
	int i;
	double dot, sign, normal1[3], normal2[3], normal[3], *Cij, IJ[3], edge_tetra_vec[3], edge_face_vec1[3], edge_face_vec2[3];


	create_vector(I_coords,J_coords,IJ);
	double edge_center[3] = { 0.5*(J_coords[0]+I_coords[0]),0.5*(J_coords[1]+I_coords[1]),0.5*(J_coords[2]+I_coords[2])};
	create_vector(edge_center,tetra_centroid,edge_tetra_vec);
	create_vector(edge_center,face1_centroid,edge_face_vec1);
	create_vector(edge_center,face2_centroid,edge_face_vec2);

	// dessa forma, normal1 e normal2 tem a mesma direcao
	cross_product(edge_tetra_vec,edge_face_vec1,normal1);
	cross_product(edge_face_vec2,edge_tetra_vec,normal2);

	for (i=0; i<3; i++){
		normal[i] = 0.5*(normal1[i]+normal2[i]);
	}

	dot_product(IJ,normal,dot);
	sign = (dot<0)?-1.0:1.0;

	// pegue Cij da aresta referente ao sub-dominio corrente e atualize-o
	EdgeInfo* edge; edge=0;
	DataCoefficients* pCoeff; pCoeff = 0;

	pMesh->getEdge(id0,id1,&edge);
	pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
	if (!pCoeff->C_ij[geom]){
		alloc_DOUBLE_vector(__LINE__,__FILE__,Cij,3);
	}
	else{
		Cij = pCoeff->C_ij[geom];
	}

	for (i=0; i<3; i++){
		Cij[i] += sign*normal[i];
	}

	pCoeff->C_ij[geom] = Cij;
	edge->pEdgeStuffes = (DataCoefficients*)pCoeff;
	pCoeff = 0;

}

void calculate_Cij(ParMeshAdapt* pMesh){
	debug_msg("\tcalculate_Cij 3D: START");

	VertexInfo *vertex[4];
	int id1, id2, id3, id4;
	double *coords_1, *coords_2, *coords_3, *coords_4;
	double *cntr_1, *cntr_2, *cntr_3, *cntr_4;
	coords_1 = 0; cntr_1 = 0;
	coords_2 = 0; cntr_2 = 0;
	coords_3 = 0; cntr_3 = 0;
	coords_4 = 0; cntr_4 = 0;

	// loop over all tetra: take each one (of six) and build Cij
	// ---------------------------------------------------------------------------------------------------------------------------------
	for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);															// index to element data structure
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];		// element taken from data structure
		int geom = pMesh->pBase_EDS[idx->base].geom;							// to which domains element belongs


		id1 = elm->id1;
		id2 = elm->id2;
		id3 = elm->id3;
		id4 = elm->id4;

		// element's vertices
		pMesh->getVertex(id1,&vertex[0]);
		pMesh->getVertex(id2,&vertex[1]);
		pMesh->getVertex(id3,&vertex[2]);
		pMesh->getVertex(id4,&vertex[3]);

		coords_1 = vertex[0]->coords;
		coords_2 = vertex[1]->coords;
		coords_3 = vertex[2]->coords;
		coords_4 = vertex[3]->coords;

		cntr_1 = elm->pFaceCentroid[0];
		cntr_2 = elm->pFaceCentroid[1];
		cntr_3 = elm->pFaceCentroid[2];
		cntr_4 = elm->pFaceCentroid[3];

		// compute Cij
		// 0-1-2 : pFaceCentroid[0]
		// 0-1-3 : pFaceCentroid[1]
		// 1-2-3 : pFaceCentroid[2]
		// 0-2-3 : pFaceCentroid[3]
		// edges: 0-1  0-2  0-3  1-2  1-3  2-3

		// edge: id1 - id2
		compute_Cij(pMesh,id1,id2,geom,coords_1,coords_2,cntr_1,cntr_2,elm->pCentroid);

		// edge: id1 - id3
		compute_Cij(pMesh,id1,id3,geom,coords_1,coords_3,cntr_1,cntr_4,elm->pCentroid);

		// edge: id1 - id4
		compute_Cij(pMesh,id1,id4,geom,coords_1,coords_4,cntr_2,cntr_4,elm->pCentroid);

		// edge: id2 - id3
		compute_Cij(pMesh,id2,id3,geom,coords_2,coords_3,cntr_1,cntr_3,elm->pCentroid);

		// edge: id2 - id4
		compute_Cij(pMesh,id2,id4,geom,coords_2,coords_4,cntr_2,cntr_3,elm->pCentroid);

		// edge: id3 - id4
		compute_Cij(pMesh,id3,id4,geom,coords_3,coords_4,cntr_3,cntr_4,elm->pCentroid);
	}

	debug_msg("\tcalculate_Cij: END");
}

void calculate_NumFacesPerSubDomain(ParMeshAdapt* pMesh, GeomData *pGCData){
	debug_msg("\tcalculate_NumFacesPerSubDomain: START");
	ElementFromFile face;
	list_ElmFromFile *face_list; face_list=0;
	std::set<int> sub_domains_set;
	std::map<int /*sub-domain geometry number (NOT physical!!!)*/,int /*counter triangles belonging to that geometry number*/> sub_domain_count;

	pMesh->getBdryFace_list(&face_list);
	alloc_INT_vector(__LINE__,__FILE__,pGCData->numDomBDRYFaces,pGCData->getNumDomains());

	// initialize counter
	for(auto it=face_list->begin(); it!=face_list->end(); it++){
		face = *it;
		auto it1 = pGCData->surface_map.find(face.geom);
		for (auto it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
			sub_domain_count[*it2] = 0;
		}
	}

	for(auto it=face_list->begin(); it!=face_list->end(); it++){
		face = *it;
		auto it1 = pGCData->surface_map.find(face.geom);
		for (auto it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
			sub_domain_count[*it2]++;
		}
	}
	face_list=0;

	int i = 0;
	//cout << "sub_domain_count.size: " << sub_domain_count.size() << endl;
	for(auto it=sub_domain_count.begin(); it!=sub_domain_count.end(); it++){
		//cout << "sub-domain " << it->first << " has " << it->second << " triangles\n";

		if (i>=pGCData->getNumDomains()){
			throw Exception(__LINE__,__FILE__,"Out of range!");
		}
		pGCData->numDomBDRYFaces[i] = it->second;
		//cout << "pGCData->numDomBDRYFaces[i]: " << pGCData->numDomBDRYFaces[i] << endl;
		i++;
	}
	//exit(1);
	debug_msg("\tcalculate_NumFacesPerSubDomain: END");
}

void calculate_Dij(ParMeshAdapt* pMesh, GeomData *pGCData){
	debug_msg("\tcalculate_Dij 3D: START");

	// allocate Dij Matrix
	// ---------------------------------------------------------------------------------------------------------------------------------
	int rows, cols = 3;
	int ndom = pGCData->getNumDomains();
	pGCData->Dij = new Matrix<double>[ndom];

	for(int i=0; i<ndom; i++){
		rows = pGCData->numDomBDRYFaces[i];
		pGCData->Dij[i].allocateMemory(rows,cols);
		pGCData->Dij[i].initialize(0);
	}

	// loop over boundary faces (internal/external) and build Dij
	// ---------------------------------------------------------------------------------------------------------------------------------
	int dom;
	int* row_counter;
	double* Dij; Dij=0;
	double v1[3], v2[3];
	bool key;

	VertexInfo *vertex[3];
	ElementFromFile face;
	list_ElmFromFile *face_list; face_list=0;
	pMesh->getBdryFace_list(&face_list);

	alloc_INT_vector(__LINE__,__FILE__,row_counter,ndom);

	for(auto it=face_list->begin(); it!=face_list->end(); it++){

		face = *it;
		auto it1 = pGCData->surface_map.find(face.geom);
		key = true;

		// the following loop has only two elements:
		for (auto it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
			dom = *it2 - 1;
			int &row = row_counter[dom];

			pGCData->getDij(dom,row,&Dij);
			if (key){

				pMesh->getVertex(face.id1,&vertex[0]);
				pMesh->getVertex(face.id2,&vertex[1]);
				pMesh->getVertex(face.id3,&vertex[2]);

				create_vector(vertex[0]->coords,vertex[1]->coords,v1);
				create_vector(vertex[1]->coords,vertex[2]->coords,v2);
				cross_product(v1,v2,Dij);

				Dij[0] *= 0.166666666666666666666666667;
				Dij[1] *= 0.166666666666666666666666667;
				Dij[2] *= 0.166666666666666666666666667;

				cout << __LINE__ << "  dom: " << dom << "   Dij: " << face.id1 << " " << face.id2 << " " << face.id3 << "   Dij: " << Dij[0] << " " << Dij[1] << " " << Dij[2] << endl;
				key = false;
			}
			else{
				pMesh->getVertex(face.id1,&vertex[0]);
				pMesh->getVertex(face.id2,&vertex[1]);
				pMesh->getVertex(face.id3,&vertex[2]);

				create_vector(vertex[0]->coords,vertex[1]->coords,v1);
				create_vector(vertex[1]->coords,vertex[2]->coords,v2);
				cross_product(v1,v2,Dij);

				Dij[0] *= -0.166666666666666666666666667;
				Dij[1] *= -0.166666666666666666666666667;
				Dij[2] *= -0.166666666666666666666666667;
			}

			row++;
		}
		Dij = 0;
	}

	double sum[3] = {0,0,0};
	double* dij; dij=0;
	for(dom=0; dom<pGCData->getNumDomains(); dom++){
		for(int i=0; i<pGCData->numDomBDRYFaces[dom]; i++){
			pGCData->getDij(dom,i,&dij);
			//cout << "dom: " << dom << " Dij = {"  << dij[0] << " " << dij[1] << " " << dij[2] << "}\n";
			sum[0] += dij[0];
			sum[1] += dij[1];
			sum[2] += dij[2];
			cout << sum[0] << " " << sum[1] << " " << sum[2] << "\n";
			dij = 0;
		}
	}

	cout << setprecision(16) << scientific;
	cout << sum[0] << " " << sum[1] << " " << sum[2] << "\n"; exit(1);

	debug_msg("\tcalculate_Dij: END");
	exit(1);
}
























