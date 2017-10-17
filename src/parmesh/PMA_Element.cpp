/*
 * PMA_Element.cpp
 *
 *  Created on: 12 de mar de 2016
 *      Author: rogerio
 */


#include "ParMeshAdapt.h"

void ParMeshAdapt::getNumElementsOfComputationalMesh(int &n) const{
	n = (int)idxlist.size();
}

void ParMeshAdapt::getBdryFace_list(list_ElmFromFile **pElmList){
	*pElmList = &tmpTrianglesList;
}

void ParMeshAdapt::getElement(IDX_str &idx, Element **pElm){
	*pElm = &this->pBase_EDS[idx.base].pLeaves[idx.level][idx.leaf];
}

void ParMeshAdapt::createVertex(int ID, VertexInfo* vinfo){
	VertexDB[ID] = vinfo;
	//cerr << "ID: "  << ID << "  maxVertex_ID: " << maxVertex_ID << endl;
	if (ID > maxVertex_ID){
		maxVertex_ID = ID;
	}
}

void ParMeshAdapt::getNumVertices(int &n){
	n = nVertices;
}

// set vertices that belong to the computational mesh
void ParMeshAdapt::updateNumVertices(){
	VertexInfo* vinfo;

	// store momently elements IDS for mapping them.
	std::set<int> ids;
	for(IdxIter it = idxlist.begin(); it!=idxlist.end(); it++){
		IDX_str *idx = &(*it);
		Element *elm = &this->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
		ids.insert(elm->id1);
		ids.insert(elm->id2);
		ids.insert(elm->id3);
	}

	for( TVertexDBIter node_it = VertexDB.begin(); node_it!=VertexDB.end(); node_it++){
		vinfo = node_it->second;
		vinfo->toCompMesh = false;
	}

	int i = 0;
	for (std::set<int>::iterator it = ids.begin(); it!=ids.end(); it++){
		vinfo = VertexDB[*it];
		vinfo->toCompMesh = true;
		vinfo->ID_mapped = i++;
	}

	nVertices = (int)ids.size();
	ids.clear();
}

void ParMeshAdapt::createVertex(int ID, double x, double y, double z, int physical, int geom){
	VertexInfo *vinfo = new VertexInfo;
	vinfo->coords = new Coords[3];
	vinfo->coords[0] = x;
	vinfo->coords[1] = y;
	vinfo->coords[2] = z;
	vinfo->geom = geom;
	vinfo->physical = physical;
	vinfo->toCompMesh = true;
	VertexDB[ID] = vinfo;

	//cerr << "ID: "  << ID << "  maxVertex_ID: " << maxVertex_ID << endl;
	if (ID > maxVertex_ID){
		maxVertex_ID = ID;
	}
}

void ParMeshAdapt::getVertex(int ID, VertexInfo** vinfo){
	std::map<int, VertexInfo*>::iterator iter = VertexDB.find(ID);
	if (iter==VertexDB.end()){
		cout << "WARNING: Vertex does not exist.\n";
		vinfo = 0;
	}
	*vinfo = iter->second;
}

void ParMeshAdapt::setVertex(int ID, int physical, int geom){
	std::map<int, VertexInfo*>::iterator iter = VertexDB.find(ID);
	if (iter==VertexDB.end()){
		cout << "WARNING: vertex not found!\n";
	}
	else{
		VertexInfo *vinfo = iter->second;
		vinfo->geom = geom;
		vinfo->physical = physical;
	}
}

void ParMeshAdapt::getVertices(const int* ID, int n, VertexInfo** vertices){
	for(int i=0; i<n; i++){
		auto iter = VertexDB.find(ID[i]);
		if (iter==VertexDB.end()){
			cout << "Vertex does not exist\n";
			exit(1);
		}
		vertices[i] = iter->second;
	}
}

void ParMeshAdapt::printElm(std::list<IDX_str> &list){
	IDX_str *idx;
	IdxIter it = list.begin();
	for(;it!=list.end(); it++){
		idx = &(*it);
		Element *pElm = &this->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
		cout << "Elm: " << pElm->id1<< " " << pElm->id2 << " " << pElm->id3 << endl;
	}
	cout << "------------------------------\n";
}

void ParMeshAdapt::printEdges(TEDGEDB &edgeDB){
	int i = 0;
	std::map<int, EdgeInfo*>::iterator it2;
	EdgeIter it1 = edgeDB.begin();
	for(;it1!=edgeDB.end();it1++){
		int id1 = it1->first;

		it2 = it1->second.begin();
		for(;it2!=it1->second.end();it2++){
			int id2 = it2->first;
			EdgeInfo* edge = it2->second;
			cout << ++i << "\tedge: [" << id1 << " " << id2;
			cout <<  "]\t" << edge->pNeighbor_1[0] << " " << edge->pNeighbor_1[1] << " " << edge->pNeighbor_1[2];
			cout <<  "\t" << edge->pNeighbor_2[0] << " " << edge->pNeighbor_2[1] << " " << edge->pNeighbor_2[2] << endl;
		}
	}
	cout << "--------------------------------------------------------\n";
}

void printElement(Element* elm){
	cout << "Elm: " << elm->id1 << " " << elm->id2 << " " << elm->id3 << "\n";
}

void printElement(Element* elm, int withVertexID){
	if (elm->id1==withVertexID){
		cout << "Elm: <" << elm->id1 << "> " << elm->id2 << " " << elm->id3 << "\n";
	}
	else if (elm->id2==withVertexID){
		cout << "Elm: " << elm->id1 << " <" << elm->id2 << "> " << elm->id3 << "\n";
	}
	else if (elm->id3==withVertexID){
		cout << "Elm: " << elm->id1 << " " << elm->id2 << " <" << elm->id3 << ">\n";
	}
}

void ParMeshAdapt::printEdge(int id1, int id2){
//	EdgeIter iter1;
//	EdgeInfo* einfo; einfo = 0;
//	bool found0, found1;
//	for (int l=0; l<=3; l++){
//		//		findEdge(id1,id2,found0,found1,iter1,pEdgeDB[0]);
//		if (found0 && found1){
//			//		getEdge(id1,id2,&einfo,pEdgeDB[0]);
//			cout << "edge: " << id1 << " " << id2 << " at level: " << l << "\tNeighbos pointers: " <<
//					einfo->pNeighbor_1 << " " << einfo->pNeighbor_2 << endl;
//			Element *elm = &this->pBase_EDS[einfo->pNeighbor_1[0]].pLeaves[einfo->pNeighbor_1[1]][einfo->pNeighbor_1[2]];
//			cout << "Elements sharing edge:\n";
//			printElement(elm);
//			if (einfo->pNeighbor_2){
//				elm = &this->pBase_EDS[einfo->pNeighbor_2[0]].pLeaves[einfo->pNeighbor_2[1]][einfo->pNeighbor_2[2]];
//				printElement(elm);
//			}
//			cout << endl;
//			break;
//		}
//	}
}


//void middlePoint(const double* p1, const double* p2, double* mp){
//	mp[0] = 0.5*(p1[0]+p2[0]);
//	mp[1] = 0.5*(p1[1]+p2[1]);
//	mp[2] = 0.5*(p1[2]+p2[2]);
//}
//
//void centroid(const double** pointers, double* pCentroid, int number_of_pointer){
//	for (int j=0; j<3; j++){
//		pCentroid[j] = .0;
//		for(int i=0; i<number_of_pointer-1; i++){
//			pCentroid[j] += pointers[i][j];
//		}
//		pCentroid[j] /= (double)(number_of_pointer);
//	}
//}
//
//void calculate_area(double* I, double* J, double* K, double &area){
//	double n[3];
//	double a[3] = { J[0]-I[0], J[1]-I[1], J[2]-I[2] };
//	double b[3] = { K[0]-I[0], K[1]-I[1], K[2]-I[2] };
//	n[0] = a[1]*b[2] - b[1]*a[2];
//	n[1] = a[2]*b[0] - a[0]*b[2];
//	n[2] = a[0]*b[1] - a[1]*b[0];
//	//area = .5*sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
//}
//
//void computeCrossProduct( double* a, double* b, double* n){
//
//}
//
//
//
//
//
//
