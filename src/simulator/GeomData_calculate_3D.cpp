/*
 * GeomData_calculate_3D.cpp
 *
 *  Created on: 27 de jun de 2016
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{

	void GeomData::calculate_1__3D(ParMeshAdapt* pMesh){
		debug_msg("calculate_1__3D: START");


		// calculate number of edge per sub-domain
		// ---------------------------------------------------------------------------------------------------------------------

		int i;
		for(i=0; i<ndom; i++){
			numDomEdges[ i ] = 0;
		}

		// count!
		std::set<int> setofflags;
		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;

				// loop over all tetrahedra sharing the same edge
				for (auto it=edge->set_of_Neighbors.begin(); it!=edge->set_of_Neighbors.end(); it++){
					int geom = pMesh->pBase_EDS[(*it)->base].geom;

					if (geom > ndom){
						throw Exception(__LINE__,__FILE__,"You cannot set a geometric flag number higher than the number of sub-domains.");
					}
					setofflags.insert( geom );
				}

				for(auto it=setofflags.begin(); it!=setofflags.end(); it++){
					numDomEdges[ *it-1 ]++;
				}
				setofflags.clear();
			}
		}

		// calculate number of boundary nodes per sub-domain
		// ---------------------------------------------------------------------------------------------------------------------

		std::map<int, std::set<int> > map_faceIDs_subdomains;
		ElementFromFile face;
		list_ElmFromFile *face_list; face_list=0;
		pMesh->getBdryFace_list(&face_list);
		for(auto it=face_list->begin(); it!=face_list->end(); it++){
			face = *it;
			auto it1 = surface_map.find(face.geom);
			for (auto it2 = it1->second.begin(); it2!=it1->second.end(); it2++){

				int sub_domain = *it2;
				map_faceIDs_subdomains[sub_domain].insert(face.id1);
				map_faceIDs_subdomains[sub_domain].insert(face.id2);
				map_faceIDs_subdomains[sub_domain].insert(face.id3);
			}
		}
		face_list=0;

		i = 0;
		auto it = map_faceIDs_subdomains.begin();
		for(; it!=map_faceIDs_subdomains.end(); it++){
			//cout << "Num. bdry nodes [" << it->first << "]:  " << it->second.size() << endl;
			numBdryNodesPerDomain[i] = (int)it->second.size();
			i++;
		}

		//exit(1);
		debug_msg("calculate_1__3D: END");
	}

	void GeomData::calculate_3__3D(ParMeshAdapt* pMesh){
		debug_msg("calculate_3__3D: START");

		double vec[3];
		double length;
		double delta_x = 1e30;
		int i,k,dom;

		VertexInfo *vertex_I; vertex_I = 0;
		VertexInfo *vertex_J; vertex_J = 0;
		int id0, id1;

		int indices[ndom];
		for (i=0; i<ndom; i++){
			indices[i] = 0;
		}

		bool key;
		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;

				key = true;
				DataCoefficients* pCoeff  = (DataCoefficients*)edge->pEdgeStuffes;
				// loop over all tetrahedra sharing the same edge
				for (k=0; k<pCoeff->num_subdomains; k++){

					// if there aremany sub-domains sharing this edge, calculate just once
					//if (key){
						id0 = iter1->first;
						id1 = iter2->first;

						vertex_I = pMesh->VertexDB[id0];
						vertex_J = pMesh->VertexDB[id1];

						// create edge vector: IJ;
						for (i=0; i<3; i++){
							vec[i] = vertex_I->coords[i] - vertex_J->coords[i];
						}

						// calculate IJ vector length
						length = .0;
						for (i=0; i<3; i++){
							length += vec[i]*vec[i];
						}
						length = sqrt(length);

						// calculate versor
						for (i=0; i<3; i++){
							vec[i] /= length;
						}

						// take smallest edge length
						if (delta_x > length){
							delta_x = length;
						}
//						key = false;
//					}

					//cout << "k: " << k << ", pCoeff->sub_domain_list[k]: "<< pCoeff->sub_domain_list[k]<< " " << endl;


					dom = pCoeff->sub_domain_list[k];
					int &row = indices[dom];
					//cout << "dom: " << dom << ", row: "<< row << " ";
					edge_length[dom].setValue(row,length);
					edge_versor[dom].setValue(row,0,vec[0]);
					edge_versor[dom].setValue(row,1,vec[1]);
					edge_versor[dom].setValue(row,2,vec[2]);
					row++;
				}
			}
		}
		setSmallestEdgeLength(delta_x);
		debug_msg("calculate_3__3D: END");
	}
}
