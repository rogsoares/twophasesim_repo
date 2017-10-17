/*
 * GeomData_initialize.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{

	/*
	 * For each domain, calculate number of:
	 *
	 * 	a) edges.
	 * 	b) boundary edges.
	 * 	c) external boundary edges.
	 * 	d) boundary nodes.
	 */
	void GeomData::calculate_1(ParMeshAdapt* pMesh){
		debug_msg("calculate_1: START");

		map_int_int edge_counter;
		map_int_int bdry_edge_counter;
		std::map<int, std::set<int> > bdry_vertices_counter;

		int physical1, physical2;

		// initialize counter
		for (EdgeIter iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;

				physical1 = pMesh->pBase_EDS[edge->pNeighbor_1[0]].physical;
				edge_counter[physical1] = 0;
				if (edge->numFaces==2){
					physical2 = pMesh->pBase_EDS[edge->pNeighbor_2[0]].physical;

					if (physical1!=physical2){
						edge_counter[physical2] = 0;
					}
				}

				if (edge->bdry){
					bdry_edge_counter[physical1] = 0;
					if (edge->numFaces==2){
						if (physical1!=physical2){
							bdry_edge_counter[physical2] = 0;
						}
					}
					else{
						numExtBdryEdges++;
					}
				}
			}
		}

		// count!
		for (EdgeIter iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;
				physical1 = pMesh->pBase_EDS[edge->pNeighbor_1[0]].physical;
				edge_counter[physical1]++;

				if (edge->numFaces==2){
					physical2 = pMesh->pBase_EDS[edge->pNeighbor_2[0]].physical;
					if (physical1!=physical2){
						edge_counter[physical2]++;
					}
				}

				if (edge->bdry){
					int id0 = iter1->first;
					int id1 = iter2->first;

					bdry_edge_counter[physical1]++;

					set_int &set_ids = bdry_vertices_counter[physical1];
					set_ids.insert(id0);
					set_ids.insert(id1);
					bdry_vertices_counter[physical1] = set_ids;

					if (edge->numFaces==2){
						if (physical1!=physical2){
							bdry_edge_counter[physical2]++;
							set_int &set_ids = bdry_vertices_counter[physical2];
							set_ids.insert(id0);
							set_ids.insert(id1);
							bdry_vertices_counter[physical2] = set_ids;
						}
					}
				}
			}
		}

		int i = 0;
		map_int_int_Iter it1 = edge_counter.begin();
		for(;it1!=edge_counter.end(); it1++){
			cout << "# edges for domain " << it1->first << ": " << it1->second << endl;
			numDomEdges[i++] = it1->second;
		}
		edge_counter.clear();

		i = 0;
		for(it1 = bdry_edge_counter.begin(); it1!=bdry_edge_counter.end(); it1++){
			cout << "# bdry edges for domain " << it1->first << ": " << it1->second << endl;
			numDomBDRYEdges[i++] = it1->second;
		}
		edge_counter.clear();

		i = 0;
		map_set_int_Iter it2 = bdry_vertices_counter.begin();
		for(;it2!=bdry_vertices_counter.end(); it2++){
			cout << "# bdry vertices for domain " << it2->first << ": " << it2->second.size() << endl;
			numBdryNodesPerDomain[i++] = (int)it2->second.size();
		}
		debug_msg("calculate_1: END");
	}

	/*
	 * For each domain, calculate number of:
	 *
	 * 	a) elements.
	 * 	b) nodes.
	 */
	void GeomData::calculate_2(ParMeshAdapt* pMesh){
		debug_msg("calculate_2: START");

		Element *elm;
		int physical;
		map_int_int elem_counter;
		std::map<int, std::set<int> > vertices_counter;


		numElem = (int)pMesh->idxlist.size();
//		for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
//			IDX_str* idx = &(*it);
//			physical = pMesh->pBase_EDS[idx->base].physical;
//			elem_counter[physical] = 0;
//		}

		// initialize
		for(int i=0; i<ndom; i++){
			elem_counter[ domains[i] ] = 0;
		}


		for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
			IDX_str* idx = &(*it);

			elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
			physical = pMesh->pBase_EDS[idx->base].physical;
			elem_counter[physical]++;

			vertices_counter[physical].insert(elm->id1);
			vertices_counter[physical].insert(elm->id2);
			vertices_counter[physical].insert(elm->id3);
			if (dim==3){
				vertices_counter[physical].insert(elm->id4);
			}
		}

		int i = 0;
		for(auto it1 = elem_counter.begin(); it1!=elem_counter.end(); it1++){
			numDomElem[i++] = it1->second;
		}

		pMesh->getNumVertices(numNodes);
		i = 0;
		for(auto it2 = vertices_counter.begin(); it2!=vertices_counter.end(); it2++){
			numNodesPerDomain[i] = (int)it2->second.size();
			i++;
		}

		debug_msg("calculate_2: END");
	}

	/*
	 * For each domain, calculate:
	 *
	 * 	a) IJ versor.
	 * 	b) edge length
	 * 	c) smallest edge length.
	 */
	void GeomData::calculate_3(ParMeshAdapt* pMesh){
		debug_msg("calculate_3: START");

		double vec[3];
		double length;
		double delta_x = 1e30;
		int i,dom;

		int physical_1, physical_2;

		VertexInfo *vertex_I; vertex_I = 0;
		VertexInfo *vertex_J; vertex_J = 0;
		int id0, id1;

		int ext_bdry_counter = 0;

		int indices[ndom];
		for (i=0; i<ndom; i++){
			indices[i] = 0;
		}

		for (EdgeIter iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (EdgeIter2 iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;

				// to which domains edge belongs
				physical_1 = pMesh->pBase_EDS[edge->pNeighbor_1[0]].physical;


				id0 = iter1->first;
				id1 = iter2->first;

				vertex_I = pMesh->VertexDB[id0];
				vertex_J = pMesh->VertexDB[id1];

				// create edge vector: IJ;
				for (i=0; i<dim; i++){
					vec[i] = vertex_I->coords[i] - vertex_J->coords[i];
				}

				// calculate IJ vector length
				length = .0;
				for (i=0; i<dim; i++){
					length += vec[i]*vec[i];
				}
				length = sqrt(length);

				// calculate versor
				for (i=0; i<dim; i++){
					vec[i] /= length;
				}

				// take smallest edge length
				if (delta_x > length){
					delta_x = length;
				}

				dom = map_domain_flags[physical_1];
				int &row_1 = indices[dom];
				//cout << "dom: " << physical_1 << ", row_1: "<< row_1 << " ";
				edge_length[dom].setValue(row_1,length);
				edge_versor[dom].setValue(row_1,0,vec[0]);
				edge_versor[dom].setValue(row_1,1,vec[1]);
				edge_versor[dom].setValue(row_1,2,vec[2]);
				row_1++;

				if (edge->numFaces==2){
					physical_2 = pMesh->pBase_EDS[edge->pNeighbor_2[0]].physical;
					if (physical_1 != physical_2){
						dom = map_domain_flags[physical_2];
						int &row_2 = indices[dom];
						//cout << "dom: " << physical_2 << ", row_2: "<< row_2;
						edge_length[dom].setValue(row_2,length);
						edge_versor[dom].setValue(row_2,0,vec[0]);
						edge_versor[dom].setValue(row_2,1,vec[1]);
						edge_versor[dom].setValue(row_2,2,vec[2]);
						row_2++;
					}
				}
				//cout << endl;

				// external boundary edge
				if (edge->numFaces==1){
					versor_ExtBdryElem[0].setValue(ext_bdry_counter,0,vec[0]);
					versor_ExtBdryElem[0].setValue(ext_bdry_counter,1,vec[1]);
					versor_ExtBdryElem[0].setValue(ext_bdry_counter,2,vec[2]);
					ext_bdry_counter++;
				}
			}
		}
		setSmallestEdgeLength(delta_x);

		debug_msg("calculate_3: END");
	}

	// calculate a unit vector normal to each external face
	void GeomData::calculate_extFaceVersor(){
	}
}
