/*
 * PMA_validation.cpp
 *
 *  Created on: 5 de abr de 2016
 *      Author: rogerio
 */

#include "ParMeshAdapt.h"

// structure to attach to mesh edges
struct Coefficients{
	double* Cij;
	double* Dij;
	double length;
};

void initCoefficients(EdgeInfo** edge);
void calcutale_Cij(double* I, double* J, double* pCentroid,double* Cij);
void calcutale_Dij(double* I, double* J, double* pCentroid, bool bdry, double* Dij);
void initVolumes(ParMeshAdapt* pMesh);

void calculateCoeffients();
string check_first(ParMeshAdapt* pMesh);
string check_second(ParMeshAdapt* pMesh);

void validateMeshIntegrity(ParMeshAdapt* pMesh){
	calculateCoeffients();
	cout << "Checking...\n;";
	cout << "\t1: " << check_first(pMesh) << ". \n";
	cout << "\t2: " << check_second(pMesh) << ". \n";
	cout << "Concluded.\n\n";
}

void calculateCoeffients(ParMeshAdapt* pMesh){

	VertexInfo* Vertex[3];
	EdgeInfo* edge[3];
	double pCentroid[3];
	double a,area;

	Coefficients* pCoeff;

	for(auto vit = pMesh->VertexDB.begin(); vit!=pMesh->VertexDB.end(); vit++){
		//vit->second->volume = 0;
	}

	// calculate: Cij, Dij and volumes of control-volumes
	for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
		TEDGEDB edgeDB = pMesh->pEdgeDB[idx->level];

		int ID[3] = {elm->id1, elm->id2, elm->id3};
		sort(ID,ID+3);

		pMesh->getVertex(ID[0],&Vertex[0]);
		pMesh->getVertex(ID[1],&Vertex[1]);
		pMesh->getVertex(ID[2],&Vertex[2]);

//		double* pCoords[3] = {Vertex[0]->coords,Vertex[1]->coords,Vertex[2]->coords};
//		centroid(elm,pCentroid,3);
//
//		pMesh->getEdge(elm->id1,elm->id2,&edge[0],edgeDB);
//		pMesh->getEdge(elm->id2,elm->id3,&edge[1],edgeDB);
//		pMesh->getEdge(elm->id1,elm->id3,&edge[2],edgeDB);

		initCoefficients(&edge[0]);
		initCoefficients(&edge[1]);
		initCoefficients(&edge[2]);

		pCoeff = (Coefficients*)edge[0]->pEdgeStuffes;
		calcutale_Cij(Vertex[0]->coords,Vertex[1]->coords,pCentroid,pCoeff->Cij);
		calcutale_Dij(Vertex[0]->coords,Vertex[1]->coords,pCentroid,edge[0]->bdry,pCoeff->Dij);

		pCoeff = (Coefficients*)edge[1]->pEdgeStuffes;
		calcutale_Cij(Vertex[1]->coords,Vertex[2]->coords,pCentroid,pCoeff->Cij);
		calcutale_Dij(Vertex[1]->coords,Vertex[2]->coords,pCentroid,edge[1]->bdry,pCoeff->Dij);

		pCoeff = (Coefficients*)edge[2]->pEdgeStuffes;
		calcutale_Cij(Vertex[0]->coords,Vertex[2]->coords,pCentroid,pCoeff->Cij);
		calcutale_Dij(Vertex[0]->coords,Vertex[2]->coords,pCentroid,edge[2]->bdry,pCoeff->Dij);

		calculate_area(Vertex[0]->coords,Vertex[1]->coords,Vertex[2]->coords,area);
		a = (double)area/3.0;
//		Vertex[0]->volume += a;
//		Vertex[1]->volume += a;
//		Vertex[2]->volume += a;
	}
}

void initCoefficients(EdgeInfo** edge){
	if ((*edge)->pEdgeStuffes){
		return;
	}

	Coefficients* pCoeff = new Coefficients;
	pCoeff->Cij = new double[3];
	pCoeff->Cij[0] = pCoeff->Cij[1] = .0;
	pCoeff->Dij = new double[3];
	pCoeff->Dij[0] = pCoeff->Dij[1] = .0;

	(*edge)->pEdgeStuffes = (Coefficients*)pCoeff;
	pCoeff = 0;
}

void calcutale_Cij(double* I, double* J, double* pCentroid,double* Cij){
	double edgeCenter[3];
	middlePoint(I,J,edgeCenter);
	double IJ[2] = {J[0]-I[0], J[1]-I[1]};										// edge vector, used as a reference vector
	double v[2] = {edgeCenter[0]-pCentroid[0],edgeCenter[1]-pCentroid[1]};		// vector: from element center to edge middle point,
																				// used as a reference vector
	if ( (v[1]*IJ[0] + (-v[0])*IJ[1]) <= .0 ){									// Cij must point as if I inside the CV and J outside
		v[0] = -v[0];
		v[1] = -v[1];
	}
	Cij[0] += v[0];
	Cij[1] += v[1];
}

void calcutale_Dij(double* I, double* J, double* pCentroid, bool bdry, double* Dij){
	if (bdry){
		Dij[0] = -(J[1]-I[1])/2.0;	// Dij vector is orthogonal to edge (it's unknown Dij orientation)
		Dij[1] =  (J[0]-I[0])/2.0;

//		// make Dij points to outside domain. First, take face that uses edge and its flag
//		int domains[2] = {0,0};
//		for (i=0; i<E_numFaces(edge); i++){
//			face = E_face(edge,i);
//			if (!face){
//				throw Exception(__LINE__,__FILE__,"Null face!\n");
//			}
//			domains[i] = EN_getFlag(face);
//		}

		double edgeCenter[3];
		middlePoint(I,J,edgeCenter);

		// vector: from element center to edge middle point, used as a reference vector
		double v[2] = {edgeCenter[0]-pCentroid[0],edgeCenter[1]-pCentroid[1]};

		// Dij must point to outside element
		if (  (Dij[0]*v[0] + Dij[1]*v[1]) <= .0 ){
			Dij[0] = -Dij[0];
			Dij[1] = -Dij[1];
		}
	}
}

//string check_first(ParMeshAdapt* pMesh){
//
//}
//
//string check_second(ParMeshAdapt* pMesh){
//
//}
