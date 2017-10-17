/*
 * EBFV1__pp3D_controlVolumes.cpp
 *
 *  Created on: 2 de dez de 2016
 *      Author: rogerio
 *
 *     Aqui, o objetivo é extrair apenas as coordenadas dos triangulos sobre os quais sao definidos os vetores Cij.
 *     Como isso, pretende-se visualizar as superficies dos volumes de controle.
 */

#include "EBFV1__pre-processors.h"

void pp_3D_controlSurface(ParMeshAdapt* pMesh, GeomData *pGCData){
	debug_msg("Pre-pp_3D_controlSurface: START");

	initCoefficients_3D(pMesh);
	calculate_centroids(pMesh);						// calculate tetrahedral's centroids
	calculate_face_centroids(pMesh);

	VertexInfo *vertex[4], * cs_vertex;
	int id1, id2, id3, id4;
	double *coords_1, *coords_2, *coords_3, *coords_4;
	double *cntr_1, *cntr_2, *cntr_3, *cntr_4;
	coords_1 = 0; cntr_1 = 0;
	coords_2 = 0; cntr_2 = 0;
	coords_3 = 0; cntr_3 = 0;
	coords_4 = 0; cntr_4 = 0;

	// gere um ID para cada aresta
	int counter = -1;
	for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
		for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
			counter++;
			EdgeInfo* edge = iter2->second;
			pMesh->getVertex(iter1->first,&vertex[0]);
			pMesh->getVertex(iter2->first,&vertex[1]);

			cs_vertex = new VertexInfo;
			cs_vertex->coords = new double[3];
			for (int i=0; i<3; i++){
				cs_vertex->coords[i] = .5*(vertex[0]->coords[i] + vertex[1]->coords[i]);
			}
			pMesh->CS_VerticesDB[counter] = cs_vertex;
			edge->MPV_id = counter;

		}
	}

	int face_ID;
	int tetra_ID = counter + 1;

	// loop over all tetra: take each one (of six) and build Cij
	// ---------------------------------------------------------------------------------------------------------------------------------
	for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);															// index to element data structure
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];		// element taken from data structure

		id1 = elm->id1; id2 = elm->id2; id3 = elm->id3; id4 = elm->id4;



		// element's vertices
		pMesh->getVertex(id1,&vertex[0]); pMesh->getVertex(id2,&vertex[1]);
		pMesh->getVertex(id3,&vertex[2]); pMesh->getVertex(id4,&vertex[3]);

		coords_1 = vertex[0]->coords; coords_2 = vertex[1]->coords;
		coords_3 = vertex[2]->coords; coords_4 = vertex[3]->coords;

		cntr_1 = elm->pFaceCentroid[0];
		cntr_2 = elm->pFaceCentroid[1];
		cntr_3 = elm->pFaceCentroid[2];
		cntr_4 = elm->pFaceCentroid[3];

		// ID sobre o tetraedro
		cs_vertex = new VertexInfo;
		cs_vertex->coords = elm->pCentroid;
		pMesh->CS_VerticesDB[tetra_ID] = cs_vertex;

		// compute Cij
		// face 1: 1-2-3 : pFaceCentroid[0]
		// face 2: 1-2-4 : pFaceCentroid[1]
		// face 3: 2-3-4 : pFaceCentroid[2]
		// face 4: 1-3-4 : pFaceCentroid[3]

		face_ID = tetra_ID + 1;
		createface(pMesh,id1,id2,cntr_1,cntr_2,tetra_ID,face_ID);
		createface(pMesh,id1,id3,cntr_1,cntr_4,tetra_ID,face_ID);
		createface(pMesh,id1,id4,cntr_2,cntr_4,tetra_ID,face_ID);
		createface(pMesh,id2,id3,cntr_1,cntr_3,tetra_ID,face_ID);
		createface(pMesh,id2,id4,cntr_2,cntr_3,tetra_ID,face_ID);
		createface(pMesh,id3,id4,cntr_3,cntr_4,tetra_ID,face_ID);
		tetra_ID = face_ID;
	}

	print_CS_VTK(pMesh);

	debug_msg("Pre-pp_3D_controlSurface: END");
	exit(1);
}

void createface(ParMeshAdapt* pMesh, int id1, int id2, double *face_centroid_1, double* face_centroid_2, int tetra_ID, int &face_ID){

	//if (id1!=2 && id2!=2) return;

	EdgeInfo* edge = 0;
	pMesh->getEdge(id1,id2,&edge);
	int ID_1A = edge->MPV_id;
	int ID_2A = tetra_ID;
	int ID_3A = face_ID;

	VertexInfo* cs_vertex = new VertexInfo;
	cs_vertex->coords = face_centroid_1;
	pMesh->CS_VerticesDB[ID_3A] = cs_vertex;

	DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
	pCoeff->face_IDs.push_back( ID_1A );
	pCoeff->face_IDs.push_back( ID_2A );
	pCoeff->face_IDs.push_back( ID_3A );

	ID_3A = ++face_ID;

	cs_vertex = new VertexInfo;
	cs_vertex->coords = face_centroid_2;
	pMesh->CS_VerticesDB[ID_3A] = cs_vertex;

	pCoeff->face_IDs.push_back( ID_1A );
	pCoeff->face_IDs.push_back( ID_2A );
	pCoeff->face_IDs.push_back( ID_3A );
	edge->pEdgeStuffes = (DataCoefficients*)pCoeff;

	face_ID++;
}

void print_CS_VTK(ParMeshAdapt* pMesh){
	ofstream fid;
	fid.open("ControlSurfaces.vtk");

	fid << setprecision(8) << fixed;

	fid << "# vtk DataFile Version 2.0\n";
	fid << "Control Volume Surfaces\n";
	fid << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
	fid << "POINTS " << pMesh->CS_VerticesDB.size() << " float\n";

	// print vertices coordinates
	// ------------------------------------------------------------------------
	VertexInfo* vertex; vertex=0;
	for(auto node_it = pMesh->CS_VerticesDB.begin(); node_it!=pMesh->CS_VerticesDB.end(); node_it++){
		vertex = node_it->second;
		fid << vertex->coords[0] << " " << vertex->coords[1] << " " << vertex->coords[2] << endl;
	}
	fid << endl;

	// count elements
	int counter = 0;
	for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
		for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
			EdgeInfo* edge = iter2->second;

			DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
			counter += (int)((int)pCoeff->face_IDs.size()/3);
		}
	}

	// print connectivities
	fid << "\nCELLS " << counter << " " << 4*counter << endl;
	for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
		for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
			EdgeInfo* edge = iter2->second;
			DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;
			for (auto it = pCoeff->face_IDs.begin(); it!=pCoeff->face_IDs.end();){
				fid << 3 << " ";
				fid << *it << " "; it++;
				fid << *it << " "; it++;
				fid << *it << "\n"; it++;
			}
		}
	}

	// print cell types
	fid << "\nCELL_TYPES " << counter << endl;
	for(int i=0; i<counter; i++){
		fid << 5 << endl;
	}

	fid.close();
}
