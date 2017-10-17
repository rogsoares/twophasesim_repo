#include "Physical.h"

Matrix<double> *pGrad_dom;		// each pointer is a matrix nx3
Matrix<double> *SwGrad_dom; 	// DO NOT REMOVE IT FROM THE CODE!!!!!!!!!!
Matrix<double> pressure;
Matrix<double> SwGrad;			// a unique matrix nnodesx3 for whole mesh
Matrix<double> Sw;				// for domain
Matrix<double> Sw_old;			// for domain
Matrix<double> Sw_tmp;
Matrix<double> p_tmp;

namespace PRS{

Physical::Physical(){
	steady_state = false;
	Allocated = false;
}

Physical::~Physical(){
}

void Physical::update_Sw_p(int n){
	allocateTemporaryData(n);
}

void Physical::allocateTemporaryData(int n){
	Sw_tmp.allocateMemory(n);
	Sw_tmp.initialize(0);
	p_tmp.allocateMemory(n);
	p_tmp.initialize(0);
}

void Physical::transferTmpData(){
	// deallocate Sw and pressure matrices first
	Sw.freeMemory();
	pressure.freeMemory();
	Sw_old.freeMemory();

	// allocate with new dimension
	int rows,cols;
	Sw_tmp.getsize(rows,cols);
	Sw.allocateMemory(rows);
	pressure.allocateMemory(rows);
	Sw_old.allocateMemory(rows);
	Sw_old.initialize(0);

	// transfer data
	for(int i=0;i<rows;i++){
		Sw.setValue(i,Sw_tmp.getValue(i));
		pressure.setValue(i,p_tmp.getValue(i));
	}
	// deallocate temporary data
	Sw_tmp.freeMemory();
	p_tmp.freeMemory();
}

void Physical::allocateData(SimulatorParameters *pSimPar, GeomData* pGCData, int numnodes){
	int ndom = pGCData->getNumDomains();
	pGrad_dom = new Matrix<double>[ndom];
	SwGrad_dom = new Matrix<double>[ndom];
	velocity = new Matrix<double>[ndom];
	for (int k=0; k<ndom; k++){
		int nrows = pGCData->getNumNodesPerDomain(k);
		int nedges = pGCData->getNumEdgesPerDomain(k);
		pGrad_dom[k].allocateMemory(nrows,3);
		pGrad_dom[k].initialize(.0);
		SwGrad_dom[k].allocateMemory(nrows,3);
		SwGrad_dom[k].initialize(.0);
		velocity[k].allocateMemory(nedges,6);
		velocity[k].initialize(.0);
	}
	nnodes = numnodes;
	SwGrad.allocateMemory(nnodes,3);
	SwGrad.initialize(.0);
	injectionWell.allocateMemory(nnodes);
	projectedSw_grad.allocateMemory(nnodes);
	nonvisc.allocateMemory(nnodes);

	// if adaptation used, Sw and pressure will be allocated in transferTmpData. This is for the very first time. Beginning of simulation.
	if (!Allocated){
		Sw.allocateMemory(nnodes);
		Sw.initialize(0);
		Sw_old.allocateMemory(nnodes);
		Sw_old.initialize(0);
		pressure.allocateMemory(nnodes);
		pressure.initialize(0);
		Allocated = true;
	}
}

void Physical::deallocateData(SimulatorParameters *pSimPar){
	int ndom = 1; //pSimPar->getNumDomains();
	for (int k=0; k<ndom; k++){
		pGrad_dom[k].freeMemory();
		SwGrad_dom[k].freeMemory();
		velocity[k].freeMemory();
	}
	SwGrad.freeMemory();
	injectionWell.freeMemory();
	projectedSw_grad.freeMemory();
	nonvisc.freeMemory();
	delete[] pGrad_dom; pGrad_dom = 0;
	delete[] SwGrad_dom; SwGrad_dom = 0;
	delete[] velocity; velocity = 0;

	if (!Allocated){
		pressure.freeMemory();
		Sw.freeMemory();
		Sw_old.freeMemory();
	}
}

void Physical::initialize(MeshData *pMData, SimulatorParameters *pSimPar, GeomData* pGCData, ParMeshAdapt* pMesh){

	debug_msg("Physical::initialize() START");

	int flag;
	int nnodes = pGCData->getNumNodes();
	allocateData(pSimPar,pGCData,nnodes);
	Swr = pSimPar->Swr();				// Irreducible water saturation
	Sor = pSimPar->Sor();				// Residual oil saturation
	mi_w = pSimPar->waterViscosity();	// water viscosity
	mi_o = pSimPar->oilViscosity();		// oil viscosity

	// transfer physical flag from edges and faces to nodes
	// ------------------------------------------------------------------------------------------------------
	int ID[3];
	VertexInfo* vertices[3];
	for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
		for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
			EdgeInfo* edge = iter2->second;
			flag = edge->physical;
			if ( pSimPar->isProductionWell(flag) || pSimPar->isInjectionWell(flag)){
				ID[0] = iter1->first;
				ID[1] = iter2->first;
				pMesh->getVertices(ID,2,vertices);
				vertices[0]->physical = edge->physical;
				vertices[1]->physical = edge->physical;
			}

		}
	}

	if (pMesh->dim()==3){
		ElementFromFile face;
		list_ElmFromFile *face_list;face_list=0;
		pMesh->getBdryFace_list(&face_list);

		for(auto it=face_list->begin(); it!=face_list->end(); it++){

			face = *it;
			flag = face.physical;
			if ( pSimPar->isProductionWell(flag) || pSimPar->isInjectionWell(flag)){
				ID[0] = face.id1;
				ID[1] = face.id2;
				ID[2] = face.id3;
				pMesh->getVertices(ID,3,vertices);
				vertices[0]->physical = face.physical;
				vertices[1]->physical = face.physical;
				vertices[2]->physical = face.physical;
			}
		}
		face_list=0;
	}

	// set initial saturation, define number of neumann nodes and number free nodes (not neumann)
	// ------------------------------------------------------------------------------------------------------
	int idx = 0;
	int k = 0;
	nfree = 0, nneumann = 0;
	VertexInfo* vinfo; vinfo = 0;
	for(auto node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
		vinfo = node_it->second;
		int flag = vinfo->physical;

		//cout << "flag: " << flag << endl;
		setSaturation(idx, pSimPar->getInitialSaturation(flag) );

		if ( pSimPar->isProductionWell(flag) ){
			nneumann++;
			injectionWell.setValue(k++,false);
		}
		if ( !pSimPar->isInjectionWell(flag) && !pSimPar->isProductionWell(flag) ){
			nfree++;
			injectionWell.setValue(k++,false);
		}
		idx++;
	}

	if(!nneumann && !pSimPar->run_Benchmark()){
		//throw Exception(__LINE__,__FILE__,"Any well detected!");
		return;
	}

	// define arrays
	// ------------------------------------------------------------------------------------------------------
	pWellsFree_index = new int[nfree];
	pWellsNeumann_index = new int[nneumann];
	idx = 0;
	int free_idx = 0;
	int neumann_idx = 0;
	for(auto node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
		vinfo = node_it->second;
		int flag = vinfo->physical;
		if ( pSimPar->isProductionWell(flag) ){
			pWellsNeumann_index[neumann_idx] = idx;
			neumann_idx++;
		}
		if ( !pSimPar->isInjectionWell(flag) && !pSimPar->isProductionWell(flag) ){
			pWellsFree_index[free_idx] = idx;
			free_idx++;
		}
		idx++;
	}
	ksModel = pSimPar->ksModel();

	debug_msg("Physical::initialize() END");
}

void Physical::setInitialSaturation(SimulatorParameters *simPar, ParMeshAdapt* pMesh){
	//		double Sw;
	//		int idx = 0;
	//		double well = false;
	//
	//		VertexInfo* vinfo; vinfo = 0;
	//		for(auto node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
	//			vinfo = node_it->second;
	//			Sw = simPar->getInitialSaturation(vinfo->physical);
	//			if ( Sw > .0 ){
	//				well = true;
	//			}
	//			setSaturation(idx,Sw);
	//			idx++;
	//		}
	//
	//		// Seeking for flagged edges
	//		// -------------------------------------------------------------------------------------------------------------------------
	//		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
	//			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){
	//				EdgeInfo* edge = iter2->second;
	//				Sw = simPar->getInitialSaturation(edge->physical);
	//				if ( Sw > .0 ){
	//					int id0 = iter1->first;
	//					int id1 = iter2->first;
	//					setSaturation(id0-1,Sw);
	//					setSaturation(id1-1,Sw);
	//					well = true;
	//				}
	//			}
	//		}
	//
	//		// Seeking for flagged faces (3-D only!)
	//		// -------------------------------------------------------------------------------------------------------------------------
	//		if (pMesh->dim()==3){
	//			ElementFromFile face;
	//			list_ElmFromFile face_list;
	//			pMesh->getBdryFace_list(face_list);
	//
	//			for(auto it=face_list.begin(); it!=face_list.end(); it++){
	//				face = *it;
	//				Sw = simPar->getInitialSaturation(face.physical);
	//				if ( Sw > .0 ){
	//					setSaturation(face.id1-1,Sw);
	//					setSaturation(face.id2-1,Sw);
	//					setSaturation(face.id3-1,Sw);
	//					well = true;
	//				}
	//			}
	//		}
	//
	//		if (!well){
	//			throw Exception(__LINE__,__FILE__,"Injection wells are missing!");
	//		}
}

double Physical::getTotalMobility(double Sw){
	double krw = get_ksw(Sw);
	double kro = get_kso(Sw);
	return krw/mi_w + kro/mi_o;
}

double Physical::getFractionalFlux(const double &Sw){
	double krw = get_ksw(Sw);
	double kro = get_kso(Sw);
	return krw/( krw + (mi_w/mi_o)*kro );
}

double Physical::getOilFractionalFlux(const double &Sw){
	double krw = get_ksw(Sw);
	double kro = get_kso(Sw);
	return (kro/mi_o)/( krw/mi_w + kro/mi_o );
}

double Physical::get_ksw(const double &Sw){
	switch(ksModel){
	case 1:
		return Sw;
	case 2:
		return pow(Sw,2);
	case 3:
		return pow((Sw - Swr)/(1. - Swr - Sor),2);
	case 4:
		return pow(Sw,4);
	default:
		throw Exception(__LINE__,__FILE__,"Unknown model for relative permeability.\n");
	}
}

double Physical::get_kso(const double &Sw){
	switch(ksModel){
	case 1:
		return 1. - Sw;
	case 2:
		return pow(1. - Sw,2);
	case 3:
		return pow((1. - Sw - Swr)/(1. - Swr - Sor),2);
	case 4:
		return (1-Sw)*(1-Sw)*(1-Sw*Sw);
	default:
		throw Exception(__LINE__,__FILE__,"Unknown model for relative permeability.\n");
	}
}

void Physical::resetNonvisc(double &alpha_max){
	alpha_max = .0;
	for (int i=0; i<nnodes; i++){
		nonvisc.setValue(i,.0);
	}
}

int Physical::getNumWells() const{
	return numWells;
}

int Physical::getNumWellNodes(int row) const{
	return numWellNodes.getValue(row,0);
}

void Physical::getFlowRate(int well_idx, int ith_node, int& row_idx, double& Qi) const{
	row_idx = idxWellNodeMat[well_idx][ith_node];
	Qi = flowrateWellNodeMat[well_idx][ith_node];
}

void Physical::getAverageMobility(int idx0_global, int idx1_global, double &MobIJ){
	double Sw_I, Sw_J;
	getSaturation(idx0_global,Sw_I);
	getSaturation(idx1_global,Sw_J);
	MobIJ = 0.5*(getTotalMobility(Sw_I) + getTotalMobility(Sw_J));
}

void Physical::calculateNodalFlowRate(GeomData* pGCData){
	//		int flag,n,i,j;
	//		double Qt, Vi, TWV;
	//		std::set<pEntity> nodes_set;
	//		std::map<int,std::set<pEntity> > wellNodes_map;
	//
	//		// search for flagged nodes: we are supposing that flags belonging to range 1 to 1999
	//		pEntity node, edge;
	//		VIter vit = M_vertexIter(theMesh);
	//		while ( (node = VIter_next(vit)) ){
	//			flag = GEN_tag( node->getClassification() );
	//			if (flag>1 && flag<1999){
	//				nodes_set = wellNodes_map[flag];
	//				nodes_set.insert(node);
	//				wellNodes_map[flag] = nodes_set;
	//			}
	//		}
	//		VIter_delete(vit);
	//
	//		EIter eit = M_edgeIter(theMesh);
	//		while ( (edge = EIter_next(eit)) ){
	//			flag = GEN_tag( edge->getClassification() );
	//			if (flag>1 && flag<1999){
	//				nodes_set = wellNodes_map[flag];
	//				nodes_set.insert((pEntity)edge->get(0,0));
	//				nodes_set.insert((pEntity)edge->get(0,1));
	//				wellNodes_map[flag] = nodes_set;
	//			}
	//		}
	//		EIter_delete(eit);
	//
	//		nodes_set.clear();
	//		std::set<pEntity>::iterator set_iter;
	//		std::map<int,std::set<pEntity> >::iterator map_iter = wellNodes_map.begin();
	//		numWells = (int)wellNodes_map.size();		// defines number of wells
	//		flowrateWellNodeMat = new double*[numWells];
	//		idxWellNodeMat = new int*[numWells];
	//		for (i=0; i<numWells; i++){
	//
	//			// first: Allocate memory
	//			flag = map_iter->first;
	//			n = (int)map_iter->second.size();			// size of node set belonging to i_th well
	//			flowrateWellNodeMat[i] = new double[n];		// allocate memory
	//			idxWellNodeMat[i] = new int[n];
	//
	//			// second: compute total well volume (TWV) (TWV = sum(Vi), Vi = volume of nodal control volume which belongs to well p)
	//			// loop over set of well p nodes
	//			TWV = .0;
	//			for (set_iter = map_iter->second.begin(); set_iter!=map_iter->second.end(); set_iter++){
	//				node = *set_iter;
	//				pGCData->getVolume_MEBFV(EN_id(node),Vi);
	//				TWV += Vi;
	//			}
	//
	//			// third: compute weighted flow rate per node for well p.
	//			set_iter = map_iter->second.begin();
	//			getTotalWellFlowRate(flag,Qt);
	//			for (j=0; j<n; j++){
	//				node = *set_iter;
	//				pGCData->getVolume_MEBFV(EN_id(node),Vi);
	//				flowrateWellNodeMat[i][j] = (double)Qt*(Vi/TWV);
	//				idxWellNodeMat[i][j] = EN_id(node)-1;
	//				set_iter++;
	//			}
	//			map_iter++;
	//		}
	//		wellNodes_map.clear();
}

void Physical::getTotalWellFlowRate(int flag, double& Qt){
	Qt = totalFlowRate_map[flag];
}
}
