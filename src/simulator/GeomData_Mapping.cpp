#include "GeomData.h"

namespace PRS{

	void GeomData::mapping(ParMeshAdapt* pMesh){
		debug_msg("mapping: START");

		const int ndom = getNumDomains();
		map_int_int mapIDtoIndex, mapBdryIDtoIndex, mapIDtoIndex_global;
		set_int dom_bdry_ids[30], dom_ids[30];

		// set a sequential numbering for each vertex ID: from 0 to n-1, where n is the number of vertices
		set_sequential_numbering(pMesh,mapIDtoIndex_global);

		// dom_ids holds node IDs for each domain
		set_nodeIDs_subdomain(pMesh,dom_ids);

		// dom_bdry_ids holds node IDs for each domain
		if (dim==2){
			set_bdrynodeIDs_subdomain(pMesh,dom_bdry_ids,mapIDtoIndex_global);
		}
		else{
			set_bdrynodeIDs_subdomain_3D(pMesh,dom_bdry_ids,mapIDtoIndex_global);
		}

		for (int k=0; k<ndom; k++){
			set_global_data(pMesh,dom_ids,k,mapIDtoIndex,mapIDtoIndex_global);
			set_bdry_data(pMesh,dom_bdry_ids,k,mapBdryIDtoIndex);
			mapping_edges(pMesh,k,mapIDtoIndex,mapIDtoIndex_global,mapBdryIDtoIndex);
		}


		if (dim==2){
			mapping_elements(pMesh,mapIDtoIndex,mapIDtoIndex_global);
		}

		mapBdryIDtoIndex.clear();
		mapIDtoIndex.clear();
		mapIDtoIndex_global.clear();

		debug_msg("mapping: END");
	}

	void GeomData::set_sequential_numbering(ParMeshAdapt* pMesh, map_int_int& mapIDtoIndex_global){
		debug_msg("\tset_sequential_numbering: START");
		int i = 0;
		int nnodes = (int)pMesh->VertexDB.size();
		alloc_INT_vector(__LINE__,__FILE__,pNodeID,nnodes);
		VertexInfo* vinfo; vinfo = 0;
		for(auto node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
			vinfo = node_it->second;
			pNodeID[i] = node_it->first;
			mapIDtoIndex_global[pNodeID[i]] = i;
			i++;
		}
		debug_msg("\tset_sequential_numbering: END");
	}

	void GeomData::set_nodeIDs_subdomain(ParMeshAdapt* pMesh, set_int *dom_ids){
		debug_msg("\tset_nodeIDs_subdomain: START");
		for (auto it = pMesh->idxlist.begin(); it!=pMesh->idxlist.end(); it++){
			IDX_str* idx = &(*it);
			int dom = map_domain_flags[ pMesh->pBase_EDS[idx->base].physical ];
			Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

			dom_ids[dom].insert(elm->id1);
			dom_ids[dom].insert(elm->id2);
			dom_ids[dom].insert(elm->id3);
			if (dim==3){
				dom_ids[dom].insert(elm->id4);
			}
		}
		debug_msg("\tset_nodeIDs_subdomain: END");
	}

	void GeomData::set_bdrynodeIDs_subdomain(ParMeshAdapt* pMesh, set_int *dom_bdry_ids, map_int_int& mapIDtoIndex_global){
		debug_msg("\tset_bdrynodeIDs_subdomain: START");

		int i = 0;
		VertexInfo* v1; v1 = 0;
		VertexInfo* v2; v2 = 0;

		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;
				if (edge->bdry){

					int id0 = iter1->first;
					int id1 = iter2->first;

					int physical1 = pMesh->pBase_EDS[edge->pNeighbor_1[0]].physical;
					int dom = map_domain_flags[physical1];

					set_int &setids1 = dom_bdry_ids[dom];
					setids1.insert(id0);
					setids1.insert(id1);

					if (edge->numFaces==2){
						int physical2 = pMesh->pBase_EDS[edge->pNeighbor_2[0]].physical;
						if (physical1!=physical2){
							dom = map_domain_flags[physical2];
							set_int &setids2 = dom_bdry_ids[dom];
							setids2.insert(id0);
							setids2.insert(id1);
						}
					}

					if (edge->numFaces==1){
						pMesh->getVertex(id0,&v1);
						pMesh->getVertex(id1,&v2);
						external_bdry_elem[0].setValue(i,0,mapIDtoIndex_global[id0]);
						external_bdry_elem[0].setValue(i,1,mapIDtoIndex_global[id1]);
						external_bdry_elem[0].setValue(i,2,v1->physical);
						external_bdry_elem[0].setValue(i,3,v2->physical);
						i++;
						v1 = 0; v2 = 0;
					}
				}
			}
		}
		debug_msg("\tset_bdrynodeIDs_subdomain: END");
	}

	void GeomData::set_bdrynodeIDs_subdomain_3D(ParMeshAdapt* pMesh, set_int *dom_bdry_ids, map_int_int& mapIDtoIndex_global){
		debug_msg("\tset_bdrynodeIDs_subdomain_3D: START");

		ElementFromFile face;
		list_ElmFromFile *face_list; face_list=0;
		pMesh->getBdryFace_list(&face_list);

		for (int dom=0; dom<_ndom; dom++){
			for(auto it=face_list->begin(); it!=face_list->end(); it++){
				face = *it;
				auto it1 = surface_map.find(face.geom);

				for (auto it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
					int sub_dom = *it2;
					if (dom==sub_dom-1){
						dom_bdry_ids[dom].insert(face.id1);
						dom_bdry_ids[dom].insert(face.id2);
						dom_bdry_ids[dom].insert(face.id3);
					}
				}
			}
		}
		face_list=0;

		debug_msg("\tset_bdrynodeIDs_subdomain_3D: END");
	}

	void GeomData::set_global_data(ParMeshAdapt* pMesh, set_int *dom_ids, int k, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global){
		debug_msg("\tset_global_data: START");

		int i = 0;
		VertexInfo* v1; v1 = 0;
		for (auto iter=dom_ids[k].begin(); iter!= dom_ids[k].end(); iter++){
			int id = *iter;
			mapIDtoIndex[id] = i;		// vertex ID
			//cout << "dom: " << k << " id " << id << " mapped as " << i << " dom_ids[k]: "  << dom_ids[k].size() << endl;
			ID[k].setValue(i,id);		// gives a sequential numbering (0,1,2,...)
			pMesh->getVertex(id,&v1);
			volume[k].setValue(i,v1->volume[domains[k]]);
			nodes[k].setValue(i,mapIDtoIndex_global[id]);
			i++;
		}

		debug_msg("\tset_global_data: END");
	}

	void GeomData::set_bdry_data(ParMeshAdapt* pMesh, set_int *dom_bdry_ids, int k, map_int_int& mapBdryIDtoIndex){
		debug_msg("\tset_bdry_data: START");

		int i = 0;
		VertexInfo* v1; v1 = 0;
		for (auto iter=dom_bdry_ids[k].begin(); iter!= dom_bdry_ids[k].end(); iter++){
			int id = *iter;					// boundary vertex ID
			mapBdryIDtoIndex[id] = i;		// gives a sequential numbering (0,1,2,...) for all boundary vertices ID for domain k
			ID_bdry[k].setValue(i,id);
			pMesh->getVertex(id,&v1);
			volume_bdry[k].setValue(i, v1->volume[domains[k]]);
			i++;
		}
		//dom_bdry_ids[k].clear();
		debug_msg("\tset_bdry_data: END");
	}

	void GeomData::mapping_edges(ParMeshAdapt* pMesh, int k, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global, map_int_int& mapBdryIDtoIndex){
		debug_msg("\tmapping_edges: START");

		if (dim==2){
			mapping_edges_2D(pMesh,k,mapIDtoIndex,mapIDtoIndex_global,mapBdryIDtoIndex);
		}
		else{
			mapping_edges_3D(pMesh,k,mapIDtoIndex,mapIDtoIndex_global);

			ElementFromFile face;
			list_ElmFromFile *face_list; face_list=0;
			pMesh->getBdryFace_list(&face_list);

			int i = 0;
			for(auto it=face_list->begin(); it!=face_list->end(); it++){
				face = *it;
				auto it1 = surface_map.find(face.geom);
				for (auto it2 = it1->second.begin(); it2!=it1->second.end(); it2++){
					int sub_dom = *it2;
					if (k==sub_dom-1){
						faces_bdry[k].setValue(i,0,mapBdryIDtoIndex[face.id1]);		// index number for boundary vertex ID for domain k
						faces_bdry[k].setValue(i,1,mapBdryIDtoIndex[face.id2]);		// index number for boundary vertex ID for domain k
						faces_bdry[k].setValue(i,2,mapBdryIDtoIndex[face.id3]);		// index number for boundary vertex ID for domain k

						faces_bdry[k].setValue(i,3,mapIDtoIndex[face.id1]);			// index number for vertex ID for domain k
						faces_bdry[k].setValue(i,4,mapIDtoIndex[face.id2]);			// index number for vertex ID for domain k
						faces_bdry[k].setValue(i,5,mapIDtoIndex[face.id3]);			// index number for vertex ID for domain k

						faces_bdry[k].setValue(i,6,mapIDtoIndex_global[face.id1]);	// global index number for vertex ID
						faces_bdry[k].setValue(i,7,mapIDtoIndex_global[face.id2]);	// global index number for vertex ID
						faces_bdry[k].setValue(i,8,mapIDtoIndex_global[face.id2]);	// global index number for vertex ID
						i++;
					}
				}
			}
			face_list=0;
		}

		debug_msg("\tmapping_edges: END");

	}
	void GeomData::mapping_edges_2D(ParMeshAdapt* pMesh, int k, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global, map_int_int& mapBdryIDtoIndex){
		int nphysical;
		int counter = 0;
		int bcounter = 0;
		int physical[2];
		VertexInfo* v1; v1 = 0;
		VertexInfo* v2; v2 = 0;

		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;

				int id0 = iter1->first;
				int id1 = iter2->first;
				pMesh->getVertex(id0,&v1);
				pMesh->getVertex(id1,&v2);

				nphysical = 1;
				physical[0] = pMesh->pBase_EDS[edge->pNeighbor_1[0]].physical;
				if (edge->numFaces==2){
					physical[1] = pMesh->pBase_EDS[edge->pNeighbor_2[0]].physical;
					if (physical[0]!=physical[1]){
						nphysical++;
					}
				}

				for (int i=0; i<nphysical; i++){
					int dom = map_domain_flags[physical[i]];
					if (dom==k){
						edges[dom].setValue(counter,0,mapIDtoIndex[id0]);					// index number for vertex ID for domain k
						edges[dom].setValue(counter,1,mapIDtoIndex[id1]); 				// index number for vertex ID for domain k
						edges[dom].setValue(counter,2,mapIDtoIndex_global[id0]);			// global index number for vertex ID
						edges[dom].setValue(counter,3,mapIDtoIndex_global[id1]);			// global index number for vertex ID
						edges[dom].setValue(counter,4,v1->physical);
						edges[dom].setValue(counter,5,v2->physical);
						//cout << "step1: ID:" << id0 << "   "  << id1 << ", local" << mapIDtoIndex[id0] << ","<< mapIDtoIndex[id1] << " global:" << mapIDtoIndex_global[id0]+1 << "  " << mapIDtoIndex_global[id1] + 1 << " in dom " << dom << " counter:" << counter+1 <<  endl;
						counter++;

						if (edge->bdry){
							edges_bdry[dom].setValue(bcounter,0,mapBdryIDtoIndex[id0]);
							edges_bdry[dom].setValue(bcounter,1,mapBdryIDtoIndex[id1]);
							edges_bdry[dom].setValue(bcounter,2,mapIDtoIndex[id0]);
							edges_bdry[dom].setValue(bcounter,3,mapIDtoIndex[id1]);
							edges_bdry[dom].setValue(bcounter,4,mapIDtoIndex_global[id0]);
							edges_bdry[dom].setValue(bcounter,5,mapIDtoIndex_global[id1]);
							bcounter++;
						}
					}
				}
			}
		}
	}

	void GeomData::mapping_edges_3D(ParMeshAdapt* pMesh, int k, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global){
		debug_msg("\tmapping_edges_3D: START");

		int id0,id1, dom, counter = 0;
		VertexInfo* v1; v1 = 0;
		VertexInfo* v2; v2 = 0;

		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;

				id0 = iter1->first;
				id1 = iter2->first;
				pMesh->getVertex(id0,&v1);
				pMesh->getVertex(id1,&v2);

				DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;

//				cout << id0 << "-" << id1 << " k:" << k <<  " pCoeff->num_subdomains: " << pCoeff->num_subdomains << "  ";
//				for (dom=0; dom<pCoeff->num_subdomains; dom++){
//					cout << pCoeff->sub_domain_list[dom] << ", ";
//				}
//				cout << endl;

				for (dom=0; dom<pCoeff->num_subdomains; dom++){
					if (pCoeff->sub_domain_list[dom]==k){
						edges[k].setValue(counter,0,mapIDtoIndex[id0]);					// index number for vertex ID for domain k
						edges[k].setValue(counter,1,mapIDtoIndex[id1]); 				// index number for vertex ID for domain k
						edges[k].setValue(counter,2,mapIDtoIndex_global[id0]);			// global index number for vertex ID
						edges[k].setValue(counter,3,mapIDtoIndex_global[id1]);			// global index number for vertex ID
						edges[k].setValue(counter,4,v1->physical);
						edges[k].setValue(counter,5,v2->physical);
//						cout << "step1: ID:" << id0 << "-"  << id1 << ", local: " <<
//								mapIDtoIndex[id0] << "-"<< mapIDtoIndex[id1] << " global:" <<
//								mapIDtoIndex_global[id0]+1 << "-" << mapIDtoIndex_global[id1] + 1 << " in dom: " <<
//								dom << " counter: " << counter+1 << " k: " << k <<
//								" pCoeff->sub_domain_list[dom]: " << pCoeff->sub_domain_list[dom] <<  endl;
						counter++;
					}
				}
			}
		}

		//cout << "counter: " << counter << endl;

		debug_msg("\tmapping_edges_3D: END");
	}

	void GeomData::mapping_elements(ParMeshAdapt* pMesh, map_int_int& mapIDtoIndex, map_int_int& mapIDtoIndex_global){
		debug_msg("\tmapping_elements: START");

		int i, dom;
		int faces_counter[ndom];
		for (i=0; i<ndom; i++){
			faces_counter[i] = 0;
		}

		for (auto it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
			IDX_str* idx = &(*it);
			dom = map_domain_flags[ pMesh->pBase_EDS[idx->base].physical ];
			Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

			i = faces_counter[dom];
			elem[dom].setValue(i,0,mapIDtoIndex[elm->id1]);
			elem[dom].setValue(i,1,mapIDtoIndex[elm->id2]);
			elem[dom].setValue(i,2,mapIDtoIndex[elm->id3]);
			elem[dom].setValue(i,3,mapIDtoIndex_global[elm->id1]);
			elem[dom].setValue(i,4,mapIDtoIndex_global[elm->id2]);
			elem[dom].setValue(i,5,mapIDtoIndex_global[elm->id3]);
		}

		debug_msg("\tmapping_elements: END");
	}
}
