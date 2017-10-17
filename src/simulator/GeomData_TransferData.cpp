#include "GeomData.h"

namespace PRS{

	// transfer only Cij vectors (3D only)
	void GeomData::transfer_Cij(ParMeshAdapt* pMesh){
		debug_msg("transfer_Cij: START");


		int dom;
		double norm;
		double *Cij = NULL;
		int* Cij_indices;

		alloc_INT_vector(__LINE__,__FILE__,Cij_indices,ndom);

		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;

				// get pointer to edge coefficients
				DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;

				for(auto it = pCoeff->C_ij.begin(); it!=pCoeff->C_ij.end(); it++){

					dom = it->first-1; //map_domain_flags[it->first];
					Cij = it->second;
					vector_norm(Cij,3,norm);
					int &row = Cij_indices[dom];

					//cout << dom << "  "  << Cij[0] << "  " << Cij[1] << " - " << row_1 << endl;
					//setCij(dom,row,Cij);
					setCij(dom,row,&Cij);
					setCij_norm(dom,row,norm);
					row++;

					//delete[] Cij; Cij=0;
					Cij = 0;
				}
				pCoeff->C_ij.clear();
			}
		}

		debug_msg("transfer_Cij: END");
	}


	void GeomData::transfer_Cij_Dij(ParMeshAdapt* pMesh){
		debug_msg("transfer_Cij_Dij: START");

		double Cij[3], Dij[3];
		double norm;
		int i;

		int Cij_indices[ndom];
		int Dij_indices[ndom];
		for (i=0; i<ndom; i++){
			Cij_indices[i] = 0;
			Dij_indices[i] = 0;
		}


		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;

				// get pointer to edge coefficients
				DataCoefficients* pCoeff = (DataCoefficients*)edge->pEdgeStuffes;

				/*
				 *  ---------------------
				 *  Cij for domain: dom1
				 *  ---------------------
				 */
				for (i=0; i<dim; i++){
					Cij[i] = pCoeff->Cij[i];
				}

				// calculate Cij norm
				norm = 0;
				for (i=0; i<dim; i++){
					norm += Cij[i]*Cij[i];
				}
				norm = sqrt(norm);


				int &row_1 = Cij_indices[pCoeff->dom1];
				//cout << pCoeff->dom1 << "  "  << Cij[0] << "  " << Cij[1] << " - " << row_1 << endl;
				setCij(pCoeff->dom1,row_1,Cij);
				setCij_norm(pCoeff->dom1,row_1,norm);
				row_1++;

				/*
				 *  ---------------------
				 *  Cij for domain: dom2
				 *  ---------------------
				 */
				if (pCoeff->dom1 != pCoeff->dom2){
					for (i=0; i<dim; i++){
						Cij[i] = pCoeff->Cij[i+2];
					}

					// calculate Cij norm
					norm = 0;
					for (i=0; i<dim; i++){
						norm += Cij[i]*Cij[i];
					}
					norm = sqrt(norm);

					int &row_2 = Cij_indices[pCoeff->dom2];
					setCij(pCoeff->dom2,row_2,Cij);
					setCij_norm(pCoeff->dom2,row_2,norm);
					row_2++;
				}


				if (edge->bdry){
					for (i=0; i<dim; i++){
						Dij[i] = pCoeff->Dij[i];
					}
					int &row_3 = Dij_indices[pCoeff->dom1];
					setDij(pCoeff->dom1,row_3,Dij);
					row_3++;

					if (edge->numFaces==2){
						if (pCoeff->dom1 != pCoeff->dom2){
							for (i=0; i<dim; i++){
								Dij[i] = -Dij[i];
							}
							int &row_4 = Dij_indices[pCoeff->dom2];
							setDij(pCoeff->dom2,row_4,Dij);
							row_4++;
						}
					}
				}
			}
		}

		debug_msg("transfer_Cij_Dij: END");
	}

	void GeomData::transfer_Volume(ParMeshAdapt* pMesh){
		debug_msg("transfer_Volume: START");

		int idx = 0;
		double volume;
		map_int_double_Iter iter;
		VertexInfo* vinfo; vinfo = 0;
		TVertexDBIter node_it = pMesh->VertexDB.begin();
		for(; node_it!=pMesh->VertexDB.end(); node_it++){
			vinfo = node_it->second;
			volume = .0;
			for(iter=vinfo->volume.begin(); iter!=vinfo->volume.end();iter++){
				volume += iter->second;
			}
			setVolume(idx,volume);
			idx++;

			//vinfo->volume.clear();
		}

		debug_msg("transfer_Volume: END");
	}

	void GeomData::transfer_mesh_3D(ParMeshAdapt* pMesh){
		debug_msg("transfer_mesh_3D: START");


		int row = 0;
		double coord[3];
		VertexInfo* vinfo; vinfo = 0;
		for(auto node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
			vinfo = node_it->second;
			for(int i=0; i<3; i++){
				pCoords->setValue(row,i,vinfo->coords[i]);
			}
			row++;

			//delete[] vinfo->coords; vinfo->coords=0;
		}

		row = 0;
		for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
			IDX_str* idx = &(*it);

			Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
			pConnectivities->setValue(row,0,elm->id1-1);	// ID-1: VTK purposes
			pConnectivities->setValue(row,1,elm->id2-1);	// ID-1: VTK purposes
			pConnectivities->setValue(row,2,elm->id3-1);	// ID-1: VTK purposes
			pConnectivities->setValue(row,3,elm->id4-1);	// ID-1: VTK purposes
			row++;

			//delete elm; elm=0;
		}

		debug_msg("transfer_mesh_3D: END");
	}

	void GeomData::transfer_Mesh(ParMeshAdapt* pMesh){
		debug_msg("transfer_Mesh: START");

		int row = 0;
		VertexInfo* vinfo; vinfo = 0;
		for(auto node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
			vinfo = node_it->second;
			for(int i=0; i<3; i++){
				pCoords->setValue(row,i,vinfo->coords[i]);
			}
			row++;
		}

		row = 0;
		for (IdxIter it = pMesh->idxlist.begin() ;it!=pMesh->idxlist.end(); it++){
			IDX_str* idx = &(*it);

			Element *elm = &pMesh->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];
			pConnectivities->setValue(row,0,elm->id1-1);	// ID-1: VTK purposes
			pConnectivities->setValue(row,1,elm->id2-1);	// ID-1: VTK purposes
			pConnectivities->setValue(row,2,elm->id3-1);	// ID-1: VTK purposes
			row++;
		}

		debug_msg("transfer_Mesh: END");
	}

	void GeomData::cleanData(ParMeshAdapt* pMesh){
		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* einfo = iter2->second;
				DataCoefficients* pCoeff = (DataCoefficients*)einfo->pEdgeStuffes;

				if (dim==2){
					delete[] pCoeff->Cij; pCoeff->Cij = 0;
					if (einfo->bdry){
						delete[] pCoeff->Dij; pCoeff->Dij = 0;
					}
				}
				else{
					for (auto it=pCoeff->C_ij.begin(); it!=pCoeff->C_ij.end();it++){
						delete[] it->second; it->second = 0;
					}
					pCoeff->C_ij.clear();
				}
				delete pCoeff; pCoeff = 0;
			}
		}
	}
}

