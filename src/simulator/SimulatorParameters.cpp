#include "SimulatorParameters.h"

namespace PRS{

	SimulatorParameters::SimulatorParameters(){
		TimeStepCounter = 0;
		exportIter = 0;
		cumTS = .0;
		stop_simulation = false;
		useHOApp = false;
		restart = false;
		vtk_step = 0;

		PVI_increment = 0.01;	// every 1% of total simulation time a new VTK file will be printed
		PVI_cumulative = .0;	// summation of all PVI_increments

		allowPrintingVTK = false;
		pctype = PCNONE;
		firstVTKupdate = true;
		EBFV1_pressureSolver_scheme = false;
		vtk_time_frequency = PVI_increment*getSimTime();
		set_adapt_occur(false);
	}

	SimulatorParameters::~SimulatorParameters(){
		MIter_RockProperties mit = mapRockProp.begin();
		for (; mit != mapRockProp.end(); mit++){
			double *K = mit->second->K;
			delete[] K; K = 0;
		}
		mapRockProp.clear();
		mapBC.clear();
	}

	void SimulatorParameters::checkPermeabilityTensor(){
		const double* K = getPermeability(3300);//BACALHO
		// Kxx*Kyy >= Kxy*Kyx
		if (K[0]*K[3] < K[1]*K[2]) throw Exception(__LINE__,__FILE__,"Permeability tensor must obey the following relation: K_xx*K_yy >= K_xy*K_yx\n");
		K_Isotropic = (K[1]*K[2] == .0)?true:false;
	}

	const double* SimulatorParameters::getPermeability(const int &dom){
		MIter_RockProperties mit = mapRockProp.find(dom);
		if (mit != mapRockProp.end()){
			return mit->second->K;
		}
		cout << "Warning: no permeability tensor associated to domain " << dom << ".\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;
		return 0;
	}

	double SimulatorParameters::getBC_Value(const int &flag){
		MapFlagIter mIter = mapBC.find(flag);
		if (mIter != mapBC.end()){
			return mIter->second->val;
		}
		cout << "Warning: getBC_Value() return null value\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;
		return 0;
	}

	bool SimulatorParameters::isNodeFree(const int &flag){
		MapFlagIter mIter = mapBC.find(flag);
		if (mIter == mapBC.end()) return true;
		else{
			BdryConditionType *bct = mIter->second;
			return ( !bct->type.compare("dirichlet") )?false:true;
		}
	}

	bool SimulatorParameters::finishSimulation(){
		std:: cout << setprecision(2) << fixed;
		if (stop_simulation || (getCumulativeSimulationTime() >= getSimTime()) ){
			std::cout << setprecision(1) <<  " " << 100 << "%\n";
			return 1;
		}
		else{
			double advance = (double)100.0*getCumulativeSimulationTime()/getSimTime();
			static double reference = 2.5;
			static bool once = true;

			if (once){
				std::cout << "Simulation progress: 0%";
				once = false;
			}

			if (advance > reference){
				std::cout << setprecision(1) <<  " " << reference << "%";
				reference += 2.5;
			}

			return 0;
		}
	}

	// physical parameters
	double SimulatorParameters::getPorosity(const int &dom){
		MIter_RockProperties mit = mapRockProp.find(dom);
		if (mit != mapRockProp.end()) return mit->second->porosity;

		cout << "Warning: porosity() return null value.";
		cout << "Domain "<< dom << " not found.\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;

		return 0;
	}

	// For a specific well, the flow rate for each node will be weighted by nodal volume.
	// Qt = (V1/Vt)Qt + (V2/Vt)Qt + ... + (Vn/Vt)Qt
	// where:
	// 		Qt 				- total flow rate (read from file)
	// 		Q1 = (V1/Vt)Qt 	- flow rate in node 1
	// 		Q2 = (V2/Vt)Qt 	- flow rate in node 2
	// 		Qn = (Vn/Vt)Qt 	- flow rate in node n
	void SimulatorParameters::WellFlowRate_weighted(ParMeshAdapt* pMesh, GeomData *pGCData){

		VertexInfo *vertex_I; vertex_I = 0;
		int well_flag, ID;
		double Vt;

		cout << "\n\nWeighting well's flow rate by node control volumes:\n";
		cout << "---------------------------------------------------------------------\n";
		cout << setprecision(6) << scientific;

		cout << "well-flag     kind       #nodes  total_volume       Flow-rate\n";
		for (auto miter=mapNodesOnWell.begin(); miter!=mapNodesOnWell.end(); miter++){
			well_flag = miter->first;
			Vt = .0;
			set_int &set_of_IDs = miter->second;
			for (auto siter = set_of_IDs.begin(); siter!=set_of_IDs.end();siter++){
				ID = *siter;
				vertex_I = pMesh->VertexDB[ID];
				for (auto iter = vertex_I->volume.begin(); iter!=vertex_I->volume.end(); iter++ ){
					Vt += iter->second;
				}
			}
			well_info = MWells[well_flag];
			well_info.wellVolume = Vt;
			MWells[well_flag] = well_info;

			if (isInjectionWell(well_flag)){
				cout << well_flag <<  "            injection  "     << set_of_IDs.size() <<   "       "  << Vt << "       "  << well_info.flowRate << endl;
			}
			if (isProductionWell(well_flag)){
				cout << well_flag <<  "            production "     << set_of_IDs.size() <<   "       "  << Vt << "       "  << well_info.flowRate << endl;
			}
		}
		cout << endl;
	}

	// get node ids associated to wells
	void SimulatorParameters::getWells(ParMeshAdapt* pMesh){

		// associate a set container to each well flag
		// EX.: two injection and one production wells provided
		// flag		-	nodes IDS associated to that flag
		//  10		-	12,34,76,285,99,43,19,...
		//  11		-	1,3114,56,2285,399,643,189,...
		//  51		-	133,36,7,2855,929,463,149,...
		// initialize map container with a set container for each flag
		//cout << __LINE__ << "\t" << __FILE__ << endl;

		for (auto miter = MWells.begin(); miter!=MWells.end(); miter++){
			int flag = miter->first;
			setNodesOnWell setNOW;
			mapNodesOnWell[flag] = setNOW;
		}


		// Production or injection wells can be assigned to geometry node, lines or faces (3-D).
		// -------------------------------------------------------------------------------------------------------------------------

		// let's start, seeking for flagged nodes
		// -------------------------------------------------------------------------------------------------------------------------
		bool found_geom_node = false;
		VertexInfo* vinfo; vinfo = 0;
		TVertexDBIter node_it;
		for(node_it = pMesh->VertexDB.begin(); node_it!=pMesh->VertexDB.end(); node_it++){
			vinfo = node_it->second;

			// check for INJECTION/PRODUCTION wells
			if ( isInjectionWell(vinfo->physical) || isProductionWell(vinfo->physical) ){
				MNOWIter = mapNodesOnWell.find(vinfo->physical);

				// insert node ID for flag
				if ( MNOWIter!=mapNodesOnWell.end() ){
					MNOWIter->second.insert(node_it->first);
					found_geom_node = true;
				}
			}
		}

		// Seeking for flagged edges
		// -------------------------------------------------------------------------------------------------------------------------
		bool found_geom_edge = false;
		for (auto iter1 = pMesh->pEdgeDB[0].begin(); iter1 != pMesh->pEdgeDB[0].end(); iter1++){
			for (auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); iter2++){

				EdgeInfo* edge = iter2->second;
				int id0 = iter1->first;
				int id1 = iter2->first;

				if ( isInjectionWell(edge->physical) || isProductionWell(edge->physical) ){
					MNOWIter = mapNodesOnWell.find(edge->physical);
					// insert node ID for flag
					if ( MNOWIter!=mapNodesOnWell.end() ){
						MNOWIter->second.insert(id0);
						MNOWIter->second.insert(id1);
						found_geom_edge = true;
					}
				}

			}
		}

		// Seeking for flagged faces (3-D only!)
		// -------------------------------------------------------------------------------------------------------------------------
		bool found_geom_face = false;
		if (pMesh->dim()==3){
			ElementFromFile face;
			list_ElmFromFile *face_list; face_list=0;
			pMesh->getBdryFace_list(&face_list);

			for(auto it=face_list->begin(); it!=face_list->end(); it++){
				face = *it;
				if ( isInjectionWell(face.physical) || isProductionWell(face.physical) ){
					MNOWIter = mapNodesOnWell.find(face.physical);

					// insert node ID for flag
					if ( MNOWIter!=mapNodesOnWell.end() ){
						MNOWIter->second.insert(face.id1);
						MNOWIter->second.insert(face.id2);
						MNOWIter->second.insert(face.id3);
						found_geom_face = true;
					}
				}
			}
			face_list=0;
		}

		cout << "\nWells found on mesh:\n";
		cout << "-------------------------------------------------------------------------------------------------------------------------\n";

		string str1 = (found_geom_node)?"YES":"NO";
		string str2 = (found_geom_edge)?"YES":"NO";
		string str3 = (found_geom_face)?"YES":"NO";

		cout << "Geometry vertices: " << str1;
		cout << "\nGeometry edges   : " << str2;
		cout << "\nGeometry face    : " << str3;
		if (!found_geom_node && !found_geom_edge && !found_geom_face){
			throw Exception(__LINE__,__FILE__,"Any well flag was passed to mesh entities.\n");
		}
	}

	bool SimulatorParameters::isInjectionWell(int flag) const{
		MWCIter mit = MWells.find( flag );
		return (mit!=MWells.end())?mit->second.isInjection:false;
	}

	bool SimulatorParameters::isProductionWell(int flag) const{
		MWCIter mit = MWells.find( flag );
		return (mit!=MWells.end())?mit->second.isProduction:false;
	}

	double SimulatorParameters::getFlowrateValue(int flag) const{
		MWCIter mit = MWells.find( flag );
		return mit->second.flowRate;
	}

	double SimulatorParameters::getTotalInjectionFlowRate() const{
		double Q = 0.0;
		for (MWCIter mwiter=MWells.begin(); mwiter!=MWells.end(); mwiter++){
			//printf("well flag %d. Flow rate: %f.  Volume: %f\n",mwiter->first,mwiter->second.flowRate,mwiter->second.wellVolume);
			if ( isInjectionWell( mwiter->first ) ){
				Q += mwiter->second.flowRate;
			}
		}
		return Q;
	}

	double SimulatorParameters::getWellVolume(int well_flag) const{
		MWCIter mit = MWells.find( well_flag );
		return mit->second.wellVolume;
	}

	/*
	 * double GeomData::getReservoirVolume() returns the total reservoir volume.
	 * To compute the initial oil volume it's necessary to take account rock po-
	 * rosity that can vary from domain to domain. So, setInitialOilVolume() function
	 * make a loop over all control volumes from all partitions (if there are two
	 * or more)
	 */
	void SimulatorParameters::setInitialOilVolume(GeomData *pGCData){
		int  dom, ndom, nnodes, node;
		double vol_total = .0, vol;
		ndom = pGCData->getNumDomains();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			double vol_domain = .0;
			for (node=0; node<nnodes; node++){
				pGCData->getVolume(dom,node,vol);
				vol_domain += vol;
			}
			vol_total += vol_domain;
		}
		_IOV = 0.2*vol_total*(1.0 - Sw_initial());		// initial oil volume (_IOV)
		//		cout << "initial oil volume (_IOV) = " << vol_total << endl;
		//		STOP();
		if (_IOV < 1e-8){
			throw Exception(__LINE__,__FILE__,"Initial Oil Volume null!");
		}
	}

	// Every N time-steps a new VTK file must be printed. It occurs at every 0.1PVI. If the sum of N time-steps exceeds 0.01PVI,
	// then the last one must be corrected to guarantee the exactness of all PVI
	void SimulatorParameters::correctTimeStep(double &timeStep){
		if (firstVTKupdate){
			updatePrintOutVTKFrequency();
		}

		double timeFrequency = getPrintOutVTKFrequency();
		double cummulative_time = timeStep + getCumulativeSimulationTime();

		if ( cummulative_time > timeFrequency ){
			timeStep = timeFrequency - getCumulativeSimulationTime();
			cummulative_time = timeFrequency;
			updatePrintOutVTKFrequency();
			allowPrintingVTK = true;
		}
	}

	void SimulatorParameters::allowPrintVTK(){
		allowPrintingVTK = true;
	}

	double SimulatorParameters::getPrintOutVTKFrequency(){
		if (firstVTKupdate){
			updatePrintOutVTKFrequency();
		}
		return vtk_time_frequency;
	}

	bool SimulatorParameters::allowPrinting_VTK(){
		return allowPrintingVTK;
	}

	void SimulatorParameters::setNotAllowPrinting_VTK(){
		allowPrintingVTK = false;
	}

	void SimulatorParameters::printOutVTK(void *pData1, void *pData2, void *pData3, pFunc_PrintVTK printVTK){
		if (allowPrintingVTK){
			int theStep = getStepOutputFile();
			incrementeStepOutputFile();
			char fname[256];
			sprintf(fname,"%s-%d.vtk",expofName.c_str(),theStep);
			PetscPrintf(PETSC_COMM_WORLD,"VTK Output: step file #%d\n",theStep);
			printVTK(pData1,pData2,pData3,fname);
			allowPrintingVTK = false;
		}
		//throw 1;
	}

	double SimulatorParameters::getInitialSaturation(int flag){
		MapFlagIter MIter = mapSaturation.find(flag);
		return (MIter == mapSaturation.end())?Sw_initial():(1.0 - Sor());
	}

	void SimulatorParameters::updatePrintOutVTKFrequency(){
		firstVTKupdate = false;

		// a new VTK file will be printed at each 0.05 PVI (5% of total simulation time)
		PVI_cumulative += getPVIincrement();

		// when next vtk file must be print out
		vtk_time_frequency = PVI_cumulative*getSimTime();
	}
}
