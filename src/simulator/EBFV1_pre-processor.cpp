/*
 * EBFV1_pre-processor.cpp
 *
 *  Created on: 26 de abr de 2016
 *      Author: rogerio
 */

#include "EBFV1__pre-processors.h"

void pre_processor(ParMeshAdapt* pMesh, GeomData *pGCData){

	getMeshDomains(pMesh,pGCData);
	if (pMesh->dim()==2){
		pp_2D(pMesh,pGCData);
	}
	else{
		pp_3D(pMesh,pGCData);
	}
}

void sort_tetra(int& id1,int& id2,int& id3,int& id4){
	int ID[4] = {id1,id2,id3,id4};
	sort(ID,ID+4);
	id1 = ID[0];
	id2 = ID[1];
	id3 = ID[2];
	id4 = ID[3];
}

void initCoefficients_3D(ParMeshAdapt* pMesh){

	int dim, dom, size, i;
	DataCoefficients* pCoeff; pCoeff = 0;

	for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
		for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
			EdgeInfo* edge = iter2->second;

			// initialize pointer to where things will be put
			// ------------------------------------------------------------------------------------------------------------------------
			pCoeff = new DataCoefficients;

			// store all tetrahedral's physical flags to edge. Repeated flags are ignored.
			// ------------------------------------------------------------------------------------------------------------------------
			std::set<int> set_of_flags;
			for (auto iter=edge->set_of_Neighbors.begin(); iter!=edge->set_of_Neighbors.end(); iter++){
				IDX_str* idx = *iter;
				set_of_flags.insert(pMesh->pBase_EDS[idx->base].geom);
			}

			// number of sub-domains sharing edge
			// ------------------------------------------------------------------------------------------------------------------------
			dim = pMesh->dim();
			pCoeff->num_subdomains = (int)set_of_flags.size();
			size = pCoeff->num_subdomains*dim;
			pCoeff->sub_domain_list = new int[pCoeff->num_subdomains];

			//cout << "pCoeff->num_subdomains: " << pCoeff->num_subdomains << endl;

			// fill sub_domain_list with sub_domain flags: C-Style index
			// ------------------------------------------------------------------------------------------------------------------------
			i = 0;
			for(auto it=set_of_flags.begin(); it!=set_of_flags.end(); it++){

				if (i>=pCoeff->num_subdomains){
					throw Exception(__LINE__,__FILE__,"Out of range!");
				}
				pCoeff->sub_domain_list[i] = *it-1;
				i++;
			}
			set_of_flags.clear();

			// initialize Cij to store all coefficients for all sub_domains tetrahedral
			// ------------------------------------------------------------------------------------------------------------------------
//			pCoeff->Cij = new double[size];
//			for(i=0; i<size; i++){
//				pCoeff->Cij[i] = 0;
//			}

			edge->pEdgeStuffes = (DataCoefficients*)pCoeff;
			pCoeff = 0;
		}
	}
}

void initCoefficients_2D(EdgeInfo** edge){
	if ((*edge)->pEdgeStuffes){
		return;
	}

	DataCoefficients* pCoeff = new DataCoefficients;

	// one domain
	if ( (*edge)->numFaces==1 ){
		pCoeff->Cij = new double[2];
		pCoeff->Cij[0] = .0;
		pCoeff->Cij[1] = .0;
		if ( (*edge)->bdry ){
			pCoeff->Dij = new double[2];
			pCoeff->Dij[0] = .0;
			pCoeff->Dij[1] = .0;
		}
	}
	// between two domains
	else{
		pCoeff->Cij = new double[4];
		pCoeff->Cij[0] = .0;
		pCoeff->Cij[1] = .0;
		pCoeff->Cij[2] = .0;
		pCoeff->Cij[3] = .0;
		if ( (*edge)->bdry ){
			pCoeff->Dij = new double[4];
			pCoeff->Dij[0] = .0;
			pCoeff->Dij[1] = .0;
			pCoeff->Dij[2] = .0;
			pCoeff->Dij[3] = .0;
		}
	}

	(*edge)->pEdgeStuffes = (DataCoefficients*)pCoeff;
	pCoeff = 0;
}

void middlePoint(const double* p1, const double* p2, double* mp){
	mp[0] = 0.5*(p1[0]+p2[0]);
	mp[1] = 0.5*(p1[1]+p2[1]);
	mp[2] = 0.5*(p1[2]+p2[2]);
}

void calculate_area(double* I, double* J, double* K, double &area){
	double n[3];
	double a[3] = { J[0]-I[0], J[1]-I[1], J[2]-I[2] };
	double b[3] = { K[0]-I[0], K[1]-I[1], K[2]-I[2] };
	n[0] = a[1]*b[2] - b[1]*a[2];
	n[1] = a[2]*b[0] - a[0]*b[2];
	n[2] = a[0]*b[1] - a[1]*b[0];
	area = .5*sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
}

void init_volumes(ParMeshAdapt* pMesh, int ndom){
	int physical;
	VertexInfo *I, *J, *K, *L;
	for (auto it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);

		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
		physical = pMesh->pBase_EDS[idx->base].physical;

		pMesh->getVertex(elm->id1,&I);
		pMesh->getVertex(elm->id2,&J);
		pMesh->getVertex(elm->id3,&K);
		I->volume[physical] = .0;
		J->volume[physical] = .0;
		K->volume[physical] = .0;

		if (pMesh->dim()==3){
			pMesh->getVertex(elm->id4,&L);
			L->volume[physical] = .0;
		}
	}
}

void calculate_triagles_centroids_and_areas(ParMeshAdapt* pMesh){
#ifdef DEBUG
	//cout << "calculate_triagles_centroids_and_areas START\n";
#endif

	int physical;
	double total = 0;
	VertexInfo *I, *J, *K;
	for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		physical = pMesh->pBase_EDS[idx->base].physical;

		pMesh->getVertex(elm->id1,&I);
		pMesh->getVertex(elm->id2,&J);
		pMesh->getVertex(elm->id3,&K);

		elm->pCentroid = new double[3];
		elm->pCentroid[0] = 0.3333333333333333*(I->coords[0]+J->coords[0]+K->coords[0]);
		elm->pCentroid[1] = 0.3333333333333333*(I->coords[1]+J->coords[1]+K->coords[1]);
		elm->pCentroid[2] = 0;

		calculate_area(I->coords,J->coords,K->coords,elm->volume);

		I->volume[physical] += (double)elm->volume/3.;
		J->volume[physical] += (double)elm->volume/3.;
		K->volume[physical] += (double)elm->volume/3.;

		total += elm->volume;
		//		printElement(elm);
		//		cout << "c: "  << elm->pCentroid[0] << " " << elm->pCentroid[1] << "\tarea: " << elm->volume << endl;
	}
	//	cout << setprecision(10) << "Total: " << total << endl;

#ifdef DEBUG
	//cout << "calculate_triagles_centroids_and_areas END\n";
#endif
}

void calculate_centroids(ParMeshAdapt* pMesh){
	debug_msg("\tcalculate_centroids: START");

	int i,j,physical;
	double total = 0;
	VertexInfo *vertex[4];
	for (auto it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
		pMesh->getVertex(elm->id1,&vertex[0]);
		pMesh->getVertex(elm->id2,&vertex[1]);
		pMesh->getVertex(elm->id3,&vertex[2]);
		pMesh->getVertex(elm->id4,&vertex[3]);

		//cout << "Tetra centroid: ";
		elm->pCentroid = new double[3];
		for(i=0; i<3; i++){
			elm->pCentroid[i] = 0;
			for(j=0; j<4; j++){
				elm->pCentroid[i] += vertex[j]->coords[i];
			}
			elm->pCentroid[i] *= 0.25;
			//cout << elm->pCentroid[i] << " ";
		}
		//cout << endl;
	}

	debug_msg("\tcalculate_centroids: END");
}

void calculate_face_centroids(ParMeshAdapt* pMesh){
	debug_msg("\tcalculate_face_centroids: START");

	int i,j,physical;
	double total = 0, C = 0.33333333333333333333333;
	VertexInfo *vertex[4];
	for (IdxIter it = pMesh->idxlist.begin(); it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		// sort tetra vertices' IDs
		sort_tetra(elm->id1,elm->id2,elm->id3,elm->id4);

		pMesh->getVertex(elm->id1,&vertex[0]);
		pMesh->getVertex(elm->id2,&vertex[1]);
		pMesh->getVertex(elm->id3,&vertex[2]);
		pMesh->getVertex(elm->id4,&vertex[3]);

		// centroids for each one of four tetrahedral faces
		elm->pFaceCentroid = new double*[4];
		for(j=0; j<4; j++){
			elm->pFaceCentroid[j] = new double[3];
		}

		// 0-1-2
		// 0-1-3
		// 1-2-3
		// 0-2-3

		//                              node               node                 node
		// face: 01                      0                   1                    2
		elm->pFaceCentroid[0][0] = C*(vertex[0]->coords[0]+vertex[1]->coords[0]+vertex[2]->coords[0]);
		elm->pFaceCentroid[0][1] = C*(vertex[0]->coords[1]+vertex[1]->coords[1]+vertex[2]->coords[1]);
		elm->pFaceCentroid[0][2] = C*(vertex[0]->coords[2]+vertex[1]->coords[2]+vertex[2]->coords[2]);

		// face: 02                      0                   1                    3
		elm->pFaceCentroid[1][0] = C*(vertex[0]->coords[0]+vertex[1]->coords[0]+vertex[3]->coords[0]);
		elm->pFaceCentroid[1][1] = C*(vertex[0]->coords[1]+vertex[1]->coords[1]+vertex[3]->coords[1]);
		elm->pFaceCentroid[1][2] = C*(vertex[0]->coords[2]+vertex[1]->coords[2]+vertex[3]->coords[2]);

		// face: 03                      1                   2                    3
		elm->pFaceCentroid[2][0] = C*(vertex[1]->coords[0]+vertex[2]->coords[0]+vertex[3]->coords[0]);
		elm->pFaceCentroid[2][1] = C*(vertex[1]->coords[1]+vertex[2]->coords[1]+vertex[3]->coords[1]);
		elm->pFaceCentroid[2][2] = C*(vertex[1]->coords[2]+vertex[2]->coords[2]+vertex[3]->coords[2]);

		// face: 04                      0                   2                    3
		elm->pFaceCentroid[3][0] = C*(vertex[0]->coords[0]+vertex[2]->coords[0]+vertex[3]->coords[0]);
		elm->pFaceCentroid[3][1] = C*(vertex[0]->coords[1]+vertex[2]->coords[1]+vertex[3]->coords[1]);
		elm->pFaceCentroid[3][2] = C*(vertex[0]->coords[2]+vertex[2]->coords[2]+vertex[3]->coords[2]);

//		cout << setprecision(6) << "Face's centroids: ";
//		cout << elm->pFaceCentroid[0][0] << " " << elm->pFaceCentroid[0][1] << " " << elm->pFaceCentroid[0][2] << ", ";
//		cout << elm->pFaceCentroid[1][0] << " " << elm->pFaceCentroid[1][1] << " " << elm->pFaceCentroid[1][2] << ", ";
//		cout << elm->pFaceCentroid[2][0] << " " << elm->pFaceCentroid[2][1] << " " << elm->pFaceCentroid[2][2] << ", ";
//		cout << elm->pFaceCentroid[3][0] << " " << elm->pFaceCentroid[3][1] << " " << elm->pFaceCentroid[3][2] << "\n";
	}

	debug_msg("\tcalculate_face_centroids: END");
}

void cleanup_auxiliary_data(ParMeshAdapt* pMesh){
	for (IdxIter it = pMesh->idxlist.begin(); it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		delete[] elm->pCentroid; elm->pCentroid=0;

		for(int j=0; j<4; j++){
			delete[] elm->pFaceCentroid[j];elm->pFaceCentroid[j]=0;
		}
		delete elm->pFaceCentroid; elm->pFaceCentroid =0;
	}

}

void calculate_volume(ParMeshAdapt* pMesh){

	debug_msg("\tcalculate_volume: START");
	int i,physical;
	VertexInfo *vertex[4];
	for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
		physical = pMesh->pBase_EDS[idx->base].physical;
		pMesh->getVertex(elm->id1,&vertex[0]);
		pMesh->getVertex(elm->id2,&vertex[1]);
		pMesh->getVertex(elm->id3,&vertex[2]);
		pMesh->getVertex(elm->id4,&vertex[3]);
		double x[4] = {vertex[0]->coords[0],vertex[1]->coords[0],vertex[2]->coords[0],vertex[3]->coords[0]};
		double y[4] = {vertex[0]->coords[1],vertex[1]->coords[1],vertex[2]->coords[1],vertex[3]->coords[1]};
		double z[4] = {vertex[0]->coords[2],vertex[1]->coords[2],vertex[2]->coords[2],vertex[3]->coords[2]};
		double a[3] = { x[1]-x[0], y[1]-y[0], z[1]-z[0] } ;
		double b[3] = { x[2]-x[0], y[2]-y[0], z[2]-z[0] } ;
		double c[3] = { x[3]-x[0], y[3]-y[0], z[3]-z[0] } ;
		double normal[3] = {a[1]*b[2] - b[1]*a[2], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
		double volume = (normal[0]*c[0] + normal[1]*c[1] + normal[2]*c[2]);

		//cout << "Tetras's nodes volumes: ";
		for(i=0; i<4; i++){
			vertex[i]->volume[physical] += 0.04166666666666667*volume;
			//cout << vertex[i]->volume[physical] << " ";
		}
		//cout << endl;
	}
	debug_msg("\tcalculate_volume: END");
}

void getMeshDomains(ParMeshAdapt* pMesh, GeomData* pGCData){
	set_int set_domain;
	int i,nelem;

	pMesh->getNumElementsOfComputationalMesh(nelem);
	for (i=0; i<nelem; i++){
		set_domain.insert( pMesh->pBase_EDS[i].physical );
	}

	pGCData->ndom = (int)set_domain.size();
	pGCData->domains = new int[pGCData->ndom];
	pGCData->setNumDomains(pGCData->ndom);

	i = 0;
	for(set_int_Iter it = set_domain.begin(); it!=set_domain.end(); it++){
		pGCData->domains[i] = *it;
		i++;
	}
	set_domain.clear();

	for (i=0; i<pGCData->ndom; i++){
		pGCData->map_domain_flags[ pGCData->domains[i] ] = i;
		//cout << pGCData->domains[i] << "  map to " << i << endl;
	}
}

void print_coefficients(ParMeshAdapt* pMesh){
	cout << setprecision(8) << scientific;
	cout << "printing Cij coefficients ...\n";
	for (EdgeIter iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
		for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){

			int id0 = iter1->first;
			int id1 = iter2->first;
			EdgeInfo* edge = iter2->second;
			DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;

			cout << id0 << " " << id1 << ": Cij at dom "<< pCoeff->dom1 << " = [" << pCoeff->Cij[0] << " " << pCoeff->Cij[1] << "], num faces: " << edge->bdry;
			if (pCoeff->dom1!=pCoeff->dom2 && pCoeff->dom2!=-1){
				cout << "\tCij at dom "<< pCoeff->dom2 << " = [" << pCoeff->Cij[2] << " " << pCoeff->Cij[3] << "]\n";
			}
			else{
				cout << endl;
			}
		}
	}
	cout << "Done\n";

	cout << "printing Dij coefficients ...\n";
	for (EdgeIter iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
		for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){

			int id0 = iter1->first;
			int id1 = iter2->first;
			EdgeInfo* edge = iter2->second;

			if (edge->bdry){
				DataCoefficients*pCoeff = (DataCoefficients*)edge->pEdgeStuffes;

				cout << id0 << " " << id1 << ": Dij at dom "<< pCoeff->dom1 << " = [" << pCoeff->Dij[0] << " " << pCoeff->Dij[1] << "]";
				if (pCoeff->dom1!=pCoeff->dom2 && pCoeff->dom2!=-1){
					cout << "\tDij at dom "<< pCoeff->dom2 << " = [" << pCoeff->Dij[2] << " " << pCoeff->Dij[3] << "]\n";
				}
				else{
					cout << endl;
				}
			}
		}
	}

	//cout << "printing Volumes ...\n";
	VertexInfo *vertex;
	std::map<int,double> volume;
	for(VIter vit = pMesh->VertexDB.begin(); vit!=pMesh->VertexDB.end(); vit++){
		vertex = vit->second;
		volume = vertex->volume;
		//cout << "Volume(s) for vertex "<< vit->first << ": ";
		for(std::map<int,double>::iterator it = volume.begin(); it!=volume.end(); it++){
			cout << "dom["<< it->first <<"] = " << it->second << "\t";
		}
		//cout << endl;
	}
	//cout << "Done\n";
}
