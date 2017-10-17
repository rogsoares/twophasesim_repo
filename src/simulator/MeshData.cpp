#include "MeshData.h"
#include "GeomData.h"

namespace PRS{

	MeshData::MeshData(){
	}

	MeshData::~MeshData(){
	}

	void MeshData::initialize(ParMeshAdapt* pMesh, SimulatorParameters* pSimPar, GeomData* pGCData){
		debug_msg("MeshData::initialize: START");

		pMesh->getNumVertices(num_global_nodes);

		set_dirichlet_free_nodes_arrays(pMesh,pSimPar);

		int dim = pGCData->getMeshDim();
		int *extnodes_ptr; extnodes_ptr = 0;
		alloc_INT_vector_init(__LINE__,__FILE__,&extnodes_ptr,dim*num_global_nodes);
		ISCreateGeneral(PETSC_COMM_WORLD,num_free_nodes,free_nodes_ptr,PETSC_COPY_VALUES,&IS_free_rows);
		ISCreateGeneral(PETSC_COMM_WORLD,num_dirichlet_nodes,dirichlet_nodes_ptr,PETSC_USE_POINTER,&IS_dirichlet_rows);
		ISCreateGeneral(PETSC_COMM_WORLD,dim*num_global_nodes,extnodes_ptr,PETSC_USE_POINTER,&IS_extendednodes_nodes);
		VecCreate(PETSC_COMM_WORLD,&Vec_Dirichlet);
		VecSetSizes(Vec_Dirichlet,PETSC_DECIDE,num_dirichlet_nodes);
		VecSetFromOptions(Vec_Dirichlet);

		//cout << "num_dirichlet_nodes: " << num_dirichlet_nodes << endl;
		for(int i=0; i<num_dirichlet_nodes; i++){
			//cout << dirichlet_nodes_ptr[i] << "  " << dirichlet_values_ptr[i] << endl;
			VecSetValue(Vec_Dirichlet,i,dirichlet_values_ptr[i],INSERT_VALUES);
		}

		cout << "\nSummary:\n";
		cout << "-------------------------------------\n";

		cout << "Number of global nodes            : " << num_global_nodes << endl;
		cout << "Number of free nodes              : " << num_free_nodes << endl;
		cout << "Number of dirichlet nodes         : " << num_dirichlet_nodes << endl;
		cout << "Dirichlet flags (value prescribed): ";
		MapFlagIter mIter = pSimPar->mapBC.begin();
		for(; mIter!=pSimPar->mapBC.end(); mIter++){
			cout << mIter->first << " (" << mIter->second->val << ")\t";
		}
		cout << "\n\n";

		debug_msg("MeshData::initialize: END");
		//exit(1);
	}

	int MeshData::set_dirichlet_free_nodes_arrays(ParMeshAdapt* pMesh, SimulatorParameters* pSimPar){

		int ID, id0, id1, i, k, j ;
		map_int_double dirichlet;

		// step 1: Geometry nodes
		VertexInfo* vinfo; vinfo = 0;
		for(auto node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
			vinfo = node_it->second;
			if ( !pSimPar->isNodeFree(vinfo->physical) ){
				int ID = node_it->first;
				dirichlet[ID] = pSimPar->getBC_Value(vinfo->physical);
			}
		}

		// step 2: Boundary edges' nodes
		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
				EdgeInfo* edge = iter2->second;
				id0 = iter1->first;
				id1 = iter2->first;
				if (edge->bdry){
					if ( !pSimPar->isNodeFree(edge->physical) ){
						dirichlet[id0] = pSimPar->getBC_Value(edge->physical);
						dirichlet[id1] = pSimPar->getBC_Value(edge->physical);
					}
				}
			}
		}

		if (pMesh->dim()==3){
			ElementFromFile face;
			list_ElmFromFile *face_list;face_list=0;
			pMesh->getBdryFace_list(&face_list);
			for(auto it=face_list->begin(); it!=face_list->end(); it++){
				face = *it;
				int flag = face.physical;
				if ( !pSimPar->isNodeFree(flag) ){
					dirichlet[face.id1] = pSimPar->getBC_Value(flag);
					dirichlet[face.id2] = pSimPar->getBC_Value(flag);
					dirichlet[face.id3] = pSimPar->getBC_Value(flag);
				}
			}
			face_list=0;
		}

		// allocate array for dirichlet nodes
		num_dirichlet_nodes = dirichlet.size();
		num_free_nodes = num_global_nodes - num_dirichlet_nodes;
		alloc_INT_vector(__LINE__,__FILE__,dirichlet_nodes_ptr,num_dirichlet_nodes);
		alloc_DOUBLE_vector(__LINE__,__FILE__,dirichlet_values_ptr,num_dirichlet_nodes);

		i = 0;
		if ( pSimPar->run_Benchmark() ){
			VertexInfo *node; node =0;
			double sol_exact;
			for(auto mit = dirichlet.begin(); mit!=dirichlet.end(); mit++){
				ID = mit->first;
				pMesh->getVertex(ID,&node);
				pSimPar->exact_solution(node->coords,sol_exact);
				dirichlet_nodes_ptr[i] = ID - 1;
				dirichlet_values_ptr[i] = sol_exact;
				i++;
			}
			node =0;
		}
		else{
			for(auto mit = dirichlet.begin(); mit!=dirichlet.end(); mit++){
				dirichlet_nodes_ptr[i] = mit->first - 1;
				dirichlet_values_ptr[i] = mit->second;
				i++;
			}
		}

		k = j = 0;
		alloc_INT_vector(__LINE__,__FILE__,free_dirichlet_global,num_global_nodes);
		for (i = 0; i<num_global_nodes; i++){
			map_int_double_Iter mit = dirichlet.find(i+1);
			if (mit != dirichlet.end()){
				free_dirichlet_global[i] = k++;
			}
			else{
				free_dirichlet_global[i] = j++;
			}
		}
		dirichlet.clear();

		// allocate array for free nodes
		i = 0;
		alloc_INT_vector(__LINE__,__FILE__,free_nodes_ptr,num_free_nodes);
		for(auto node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
			vinfo = node_it->second;
			if ( pSimPar->isNodeFree(vinfo->physical) ){
				free_nodes_ptr[i] = node_it->first-1;
				i++;
			}
		}
		return 0;
	}

	void MeshData::get_node_condition(int ID, int &pos){
		pos = free_dirichlet_global[ID-1];
	}

	void MeshData::get_dirichlet_value(int i, double& val){
		val = dirichlet_values_ptr[i];
	}

	void MeshData::get_num_global_nodes(int &ngn){
		ngn = num_global_nodes;
	}

	void MeshData::get_num_free_nodes(int &nfn){
		nfn = num_free_nodes;
	}

	void MeshData::get_num_dirichlet_nodes(int &ndn){
		ndn = num_dirichlet_nodes;
	}

	void MeshData::get_free_nodes_ptr(int** nfn_ptr){
		*nfn_ptr = free_nodes_ptr;
	}

	void MeshData::get_dirichlet_nodes_ptr(int **ndn_ptr){
		*ndn_ptr = dirichlet_nodes_ptr;
	}

	void MeshData::get_IS_dirichlet_nodes(IS *dir_rows){
		*dir_rows = IS_dirichlet_rows;
	}

	void MeshData::get_IS_free_nodes(IS *frows){
		*frows = IS_free_rows;
	}

	void MeshData::get_IS_extendednodes_nodes(IS *extendednodes_nodes){
		*extendednodes_nodes = IS_extendednodes_nodes;
	}

	void MeshData::get_Vec_dirichlet(Vec *v){
		*v = Vec_Dirichlet;
	}
}
