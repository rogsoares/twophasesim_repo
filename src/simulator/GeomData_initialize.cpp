/*
 * GeomData_initialize.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{

	void GeomData::initialize(ParMeshAdapt* pMesh){
		debug_msg("GeomData::initialize START");

		setMeshDim(pMesh->dim());
		elemtype = (dim==2)?3:4;

		int nnodes_global;
		pMesh->getNumVertices(nnodes_global);

		alloc_INT_vector(__LINE__,__FILE__,numDomEdges,_ndom);
		alloc_INT_vector(__LINE__,__FILE__,numDomBDRYEdges,_ndom);
		alloc_INT_vector(__LINE__,__FILE__,numBdryNodesPerDomain,_ndom);
		alloc_INT_vector(__LINE__,__FILE__,numDomElem,_ndom);
		alloc_INT_vector(__LINE__,__FILE__,numNodesPerDomain,_ndom);

		if (dim==2){
			calculate_1(pMesh);
		}
		else{
			calculate_1__3D(pMesh);
		}

		calculate_2(pMesh);
		allocatePointers(nnodes_global,dim);

		if (dim==2){
			calculate_3(pMesh);
			transfer_Cij_Dij(pMesh);
			transfer_Mesh(pMesh);
		}
		else{
			calculate_3__3D(pMesh);
			transfer_Cij(pMesh);
			transfer_mesh_3D(pMesh);
		}
		transfer_Volume(pMesh);

		cout << "Number of nodes per sub-domain   : ";
		for (int i=0; i<ndom; i++){
			cout << "(" << domains[i] << ")  "  << numNodesPerDomain[i] << "\t";
		}
		cout << endl;

		cout << "Number of edges per sub-domain: ";
		for (int i=0; i<ndom; i++){
			cout << "(" << domains[i] << ") "  << numDomEdges[i] << "\t";
		}
		cout << endl;

		cout << "Number of faces per sub-domain   : ";
		for (int i=0; i<ndom; i++){
			cout << "(" << domains[i] << ")  "  << numDomBDRYFaces[i] << "\t";
		}
		cout << endl;

		cout << "Number of elements per sub-domain: ";
		for (int i=0; i<ndom; i++){
			cout << "(" << domains[i] << ") "  << numDomElem[i] << "\t";
		}
		cout << endl;

		//calculate_extFaceVersor(pMesh);
		mapping(pMesh);

		debug_msg("GeomData::initialize END");
	}

	void GeomData::read_geometry(const char* filename){

		ifstream fid;
		char fname[512];
		sprintf(fname,"%s/geometry.dat",filename);
		open_file(fid,fname,__LINE__,__FILE__);
		seek_position_on_file(fid,"Digite abaixo o que se pede:");

		fid >> num_surfaces;
		int i, surf_number, num_sub_domains, sub_domain_1, sub_domain_2;
		for (i=0; i<num_surfaces; i++){
			fid >> surf_number >> num_sub_domains >> sub_domain_1;
			sub_domain_2 = 0;
			if (num_sub_domains==2){
				fid >> sub_domain_2;
			}

			std::set<int> sub_domains_number_set;
			sub_domains_number_set.insert(sub_domain_1);
			if (num_sub_domains==2){
				sub_domains_number_set.insert(sub_domain_2);
			}
			surface_map.insert( pair<int, std::set<int> >(surf_number,sub_domains_number_set));
		}

//		for (auto it1 = surface_map.begin(); it1!= surface_map.end(); it1++){
//			cout << it1->first << ":  ";
//			for (auto it2 = it1->second.begin(); it2!= it1->second.end(); it2++){
//				cout << *it2 << " ";
//			}
//			cout << endl;
//		}
	}
}
