/*
 * EBFV1_AssemblyMatVec.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: rogerio
 */

#include "EBFV1_elliptic.h"

namespace PRS{

	int EBFV1_elliptic::wells(GeomData* pGCData, SimulatorParameters* pSimPar, MeshData* pMData){

		if (pSimPar->run_Benchmark() && pSimPar->ss_term){
			const double* coords = NULL;
			double volume, sol;
			int num_free_nodes;
			int *free_nodes_ptr = NULL;
			pMData->get_num_free_nodes(num_free_nodes);
			pMData->get_free_nodes_ptr(&free_nodes_ptr);
			for (int i=0; i<num_free_nodes; i++){
				pGCData->getVolume(free_nodes_ptr[i],volume);
				pGCData->getCoordinates(free_nodes_ptr[i],&coords);
				pSimPar->ss_term(coords,sol);
				VecSetValue(pDStruct->RHS,i,volume*sol,ADD_VALUES);
			}
			coords = NULL;
			VecAssemblyBegin(pDStruct->RHS);
			VecAssemblyEnd(pDStruct->RHS);
		}
		else{
			int node_ID, row, num_well_nodes, well_flag;
			double Vt, Vi, Qi, Qt;

			// for each well flag
			for (map_set_int_Iter mit = pSimPar->mapNodesOnWell.begin(); mit!=pSimPar->mapNodesOnWell.end(); mit++){
				well_flag = mit->first;
				set_int &set_node_IDs = mit->second;

				// only include nodes flagged as Neumann.
				if (pSimPar->isNodeFree(well_flag)){
					Qt = pSimPar->getFlowrateValue(well_flag);		// source/sink term
					Vt = pSimPar->getWellVolume(well_flag);			// total well volume
					num_well_nodes = mit->second.size();			// number of nodes flagged as well

					// get all flagged node IDs for that well
					if (!num_well_nodes){
						throw Exception(__LINE__,__FILE__,"No wells found!");
					}

					// get all nodes'ID for well: well_flag
					for (set_int_Iter sit = set_node_IDs.begin(); sit!=set_node_IDs.end(); sit++){
						node_ID = *sit;						// node ID
						pGCData->getVolume(node_ID-1,Vi);	// node's volume
						Qi = Qt*(Vi/Vt);					// fractional flow rate
						pMData->get_node_condition(node_ID,row);
						VecSetValue(pDStruct->RHS,row,Qi,ADD_VALUES);
					}
				}
			}
			VecAssemblyBegin(pDStruct->RHS);
			VecAssemblyEnd(pDStruct->RHS);
		}
		return 0;
	}
}
