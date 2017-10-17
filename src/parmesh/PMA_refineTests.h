/*
 * PMA_refineTests.h
 *
 *  Created on: 12 de mar de 2016
 *      Author: rogerio
 */

#ifndef SRC_PMA_REFINETESTS_H_
#define SRC_PMA_REFINETESTS_H_


#include "ParMeshAdapt.h"


// refine (level=1) inside region defined by two circumferences
void refine_ring(ParMeshAdapt* pMesh);

void setElementsToRefine(ParMeshAdapt* pMesh);
void setElementsToUnRefine(ParMeshAdapt* pMesh);

// unrefine (level=1) inside region defined by two circumferences
void setElementsWithUnRefinementLevel(ParMeshAdapt* pMesh);


void validateMeshIntegrity(ParMeshAdapt* pMesh);

#endif /* SRC_PMA_REFINETESTS_H_ */
