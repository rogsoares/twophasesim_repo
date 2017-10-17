/*
 * PMA_refineTests.cpp
 *
 *  Created on: 12 de mar de 2016
 *      Author: rogerio
 */

#include "PMA_refineTests.h"

void refine_ring(ParMeshAdapt* pMesh){

	// REFINEMENT: refine mesh inside ring region
	// =========================================================
	//pMesh->write();
	pMesh->statistics();


	setElementsToRefine(pMesh);
	pMesh->OneLevelSmoothing();					// before refine, neighbor elements must be one level of refinement
												// difference smoothing process base on number of subdivisions
	pMesh->refine();							// now refine!

	for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
		elm->toRefine = false;
		elm->numrefine = 0;
		elm->count = 0;
	}

	pMesh->RefinementSmoothing_part3();			// after refinement, mesh must be smoothing again, but now it's based on the number of refinement levels
	pMesh->buildGreenElements();				// special refinement: eliminate hang-nodes
	pMesh->updateNumVertices();


	pMesh->write();
	for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
		elm->toRefine = false;
		elm->numrefine = 0;
		elm->count = 0;
	}
	//validateMeshIntegrity(pMesh);
	//pMesh->statistics();


	// UNREFINEMENT: Unrefine mesh inside ring region
	// =========================================================

	// set elements of the computational mesh with unrefinement level
	setElementsWithUnRefinementLevel(pMesh);
	pMesh->write();
	pMesh->OneLevelSmoothing();
	pMesh->write();
	pMesh->setElementsTobeRemoved();
	pMesh->updateNumVertices();

	pMesh->write();
}

void setElementsToRefine(ParMeshAdapt* pMesh){

	// Set elements to be refined
	cout << "// Set elements to be refined\n";

	VertexInfo *v1, *v2, *v3;
	int count=0;

	for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
		IDX_str* idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		// get element vertices
		pMesh->getVertex(elm->id1,&v1);
		pMesh->getVertex(elm->id2,&v2);
		pMesh->getVertex(elm->id3,&v3);

		// calculate element baricenter
		double bar[2] = {(v1->coords[0]+v2->coords[0]+v3->coords[0])/3.,(v1->coords[1]+v2->coords[1]+v3->coords[1])/3.};

		double x_0 = 0.5, y_0 = 0.5;
		double x = bar[0];
		double y = bar[1];
		double R_max = 0.3, R_min = 0.2;
		double xy = (x-x_0)*(x-x_0) + (y-y_0)*(y-y_0);
		if (  xy > R_min*R_min && xy < R_max*R_max  ){
			elm->toRefine = true;
			if (x<.5 && y<0.5){
				elm->numrefine = 1;
			}
			else if (x>.5 && y>0.5){
				elm->numrefine = 1;
			}
			else if (x<.5 && y>0.5){
				elm->numrefine = 2;
			}
			else if (x>.5 && y<0.5){
				elm->numrefine = 3;
			}
			elm->toUnRefine = false;
			count++;
		}
		else{
			elm->toRefine = false;
		}
		elm->count = 0;
	}
}

// walk through computational mesh and identify those which level if 1
void setElementsWithUnRefinementLevel(ParMeshAdapt* pMesh){
	int count=0;
	VertexInfo *v1, *v2, *v3;

	for (IdxIter it = pMesh->idxlist.begin(); it!=pMesh->idxlist.end(); it++){
		IDX_str *idx = &(*it);
		Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		// get element vertices
		pMesh->getVertex(elm->id1,&v1);
		pMesh->getVertex(elm->id2,&v2);
		pMesh->getVertex(elm->id3,&v3);

		// calculate element baricenter
		double bar[2] = {(v1->coords[0]+v2->coords[0]+v3->coords[0])/3.,(v1->coords[1]+v2->coords[1]+v3->coords[1])/3.};

		double x_0 = 0.55, y_0 = 0.5;
		double x = bar[0];
		double y = bar[1];
		double R_max = 0.3, R_min = 0.2;
		double xy = (x-x_0)*(x-x_0) + (y-y_0)*(y-y_0);
		if (  xy > R_min*R_min && xy < R_max*R_max && idx->level ){
			elm->toUnRefine = true;
			if (x<.5 && y<0.5){
				elm->numrefine = -1;
			}
			else if (x>.5 && y>0.5){
				elm->numrefine = -1;
			}
			else if (x<.5 && y>0.5){
				elm->numrefine = -2;
			}
			else if (x>.5 && y<0.5){
				elm->numrefine = -3;
			}
			elm->toRefine = false;
			count++;
		}
		else{
			elm->numrefine = 0;
		}
	}
	cout << "count unrefine: " << count << endl;
}
