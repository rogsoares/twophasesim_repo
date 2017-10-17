#include "EBFV1_elliptic.h"
//#include "calculate_results_Benchmark3D.h"

double cpu_time;

namespace PRS{

EBFV1_elliptic::EBFV1_elliptic(){
	pMAS =0;
	Perform_Assembling = true;
	pDStruct = new Data_struct;
}

EBFV1_elliptic::~EBFV1_elliptic(){
}

// solves system of equation for pressure field
double EBFV1_elliptic::solver(MeshData *pMData, GeomData *pGCData, SimulatorParameters *pSimPar, Physical *pPPData){
	debug_msg("EBFV1_elliptic::solver(): START");

	// if you are interested in viewing the exact solution only skip calculations
	// ---------------------------------------------------------------------------------------------------
	int num_global_nodes, num_free_nodes;
	assembly_EFG_RHS(pMData,pGCData,pPPData,pSimPar);
	//pGCData->deallocatePointers(3);
	pMData->get_num_free_nodes(num_free_nodes);
	pMData->get_num_global_nodes(num_global_nodes);

	double t1 = MPI_Wtime();
	setMatrixFreeOperation(num_free_nodes,num_global_nodes,pGCData->getMeshDim());
	double t2 = MPI_Wtime();
	cpu_time = t2 - t1;

	// update pressure field with calculated or exact solution
	// ---------------------------------------------------------------------------------------------------
	updatePressure(pMData,pGCData,pPPData,pSimPar);
	//calculatePressureGradient(pGCData,pPPData);

	//free_matrices(pGCData);

	debug_msg("EBFV1_elliptic::solver(): END");
	return 0;
}

double EBFV1_elliptic::updatePressure(MeshData *pMData, GeomData *pGCData, Physical *pPPData, SimulatorParameters *pSimPar){
	debug_msg("EBFV1_elliptic::updatePressure(): START");

	const double* coords = NULL;
	double sol, sol_exact, sol_calc, volume;
	int i, num_free_nodes, num_dirichlet_nodes, num_global_nodes;
	int *free_nodes_ptr; free_nodes_ptr = 0;
	int *pDirichletNodes; pDirichletNodes=0;

	// update pressure field with exact solution only and then leave function
	// ---------------------------------------------------------------------------------------------------
	//		if ( pSimPar->exactSolutionExist() ){
	//		{
	//			double x,y,z;
	//			double minimum = 1.0e+8;
	//			pMData->get_num_global_nodes(num_global_nodes);
	//			for(int i=0;i<num_global_nodes;i++){
	//				pGCData->getCoordinates(i,&coords);
	//				pSimPar->exact_solution(coords,sol);
	//				if (sol < minimum){
	//					minimum = sol;
	//					x = coords[0];
	//					y = coords[1];
	//					z = coords[2];
	//				}
	//			}
	//			cout << "Minimum exact: " << minimum << "  in:" << x << " " << y << " " << z << endl;
	//		}
	//		exit(1);
	//			return 0;
	//		}

	// update pressure field with calculated solution
	// ---------------------------------------------------------------------------------------------------
	pMData->get_num_free_nodes(num_free_nodes);
	pMData->get_free_nodes_ptr(&free_nodes_ptr);
	for(i=0; i<num_free_nodes; i++){
		VecGetValues(pDStruct->solution,1,&i,&sol);
		pPPData->setPressure(free_nodes_ptr[i],sol);
	}

	// set prescribed (dirichlet) values to pressure field
	// ---------------------------------------------------------------------------------------------------
	pMData->get_num_dirichlet_nodes(num_dirichlet_nodes);
	pMData->get_dirichlet_nodes_ptr(&pDirichletNodes);
	for (i=0; i<num_dirichlet_nodes; i++){
		if (pSimPar->use_nonhomogneous_bc()){
			pGCData->getCoordinates(pDirichletNodes[i],&coords);
			pSimPar->exact_solution(coords,sol);
		}
		else{
			pMData->get_dirichlet_value(i,sol);
		}
		pPPData->setPressure(pDirichletNodes[i],sol);
	}


	// For NIKITIN analysis
	// -------------------------------------------------------------------
	double max_calc = -1e+8;
	double min_calc = 1e+8;
	//		pMData->get_num_free_nodes(num_free_nodes);
	//		pMData->get_free_nodes_ptr(&free_nodes_ptr);
	//		for(i=0; i<num_free_nodes; i++){
	//			VecGetValues(pDStruct->solution,1,&i,&sol_calc);
	//			// max/min
	//			if ( sol_calc > max_calc ){
	//				max_calc = sol_calc;
	//			}
	//			if ( sol_calc < min_calc ){
	//				min_calc = sol_calc;
	//			}
	//		}
	// NIKITIN OUTPUT
	// -------------------------------------------------------------------
	//			string line;
	//			bool already_edited = false;
	//			ifstream fid_read;
	//			fid_read.open("Nikitin.txt");
	//			fid_read >> line;
	//			if ( strcmp(line.c_str(),"node min max") ){ already_edited = true; }
	//			fid_read.close();
	//			ofstream fid;
	//			fid.open("Nikitin.txt", ios_base::app);
	//			if (!already_edited){fid << "node min max\n";}
	//			fid << num_free_nodes << " " << min_calc << " " << max_calc << endl;
	//			fid.close();
	// -------------------------------------------------------------------
	// -------------------------------------------------------------------

	// calculate L2 norm of error and max/min solutions
	// ---------------------------------------------------------------------------------------------------
	if ( pSimPar->run_Benchmark() ){
		double erms, el2v, el2, el2vr, el2r;
		double error;
		double error_pow;
		double error_max = 0;
		double el2v_sum = 0;
		double el2ve_sum = 0;
		double el2_sum = 0;
		double erms_sum = 0;
		double exact_sum = 0; // it serves for relative errors calculations!!!!!

		double max_exact = -1e+8;
		double min_exact = +1e+8;

		pMData->get_num_free_nodes(num_free_nodes);
		pMData->get_free_nodes_ptr(&free_nodes_ptr);
		for(i=0; i<num_free_nodes; i++){
			pPPData->getPressure(free_nodes_ptr[i],sol_calc);
			pGCData->getCoordinates(free_nodes_ptr[i],&coords);
			pGCData->getVolume(free_nodes_ptr[i],volume);
			pSimPar->exact_solution(coords,sol_exact);


			// error summation weighted by node volume
			error = sol_exact - sol_calc;
			error_pow = std::pow(error,2);
			exact_sum += std::pow(sol_exact,2);

			el2ve_sum += std::pow(sol_exact,2)*volume;
			el2v_sum += error_pow*volume;
			el2_sum += error_pow;
			erms_sum += error_pow;


			if (abs(error) > error_max){ error_max = abs(error); }
			if ( sol_exact > max_exact ){ max_exact = sol_exact; }
			if ( sol_exact < min_exact ){ min_exact = sol_exact; }
			if ( sol_calc > max_calc ){ max_calc = sol_calc; }
			if ( sol_calc < min_calc ){ min_calc = sol_calc; }
		}

		el2 = sqrt(el2_sum);
		el2r = sqrt(el2_sum/exact_sum);
		el2v = sqrt(el2v_sum);		// ponderado pelo volume
		el2vr = sqrt(el2v_sum/el2ve_sum);
		erms = sqrt(erms_sum/num_free_nodes);

		// Other benchmark cases OUTPUT
		// ---------------------------------------------------------------
		ofstream fid;
		string filename;
		switch (pSimPar->getCaseProblem()){
		case CASE_1:
			filename = "output_data_analysis/main/HomogeneoAnisotropico/case_1.dat"; break;
		case CASE_5:
			filename = "output_data_analysis/main/HeterogeneoAnisotropico/case_5.dat"; break;
		case OBLIQUEDRAIN:
			filename = "output_data_analysis/main/LinearPreserving/oblique_drain.dat"; break;

		}
		fid.open(filename, ios_base::app);

		if (!fid.is_open()){
			cout << "Error opening file: " << filename << "\nExiting...\n";
			exit(1);
		}

		pMData->get_num_global_nodes(num_global_nodes);
		fid << "nodes CPU el2 el2r el2v el2vr emax erms minex maxex minca maxca\n";
		fid << setprecision(8) << fixed << scientific;
		fid << num_free_nodes << " " << cpu_time << " "
			<< el2 << " " << el2r << " "
			<< el2v << " " << el2vr << " "
			<< error_max << " " << erms << " "
			<< min_exact << " " << max_exact << " "
			<< min_calc << " " << max_calc << endl;
	}

	debug_msg("EBFV1_elliptic::updatePressure(): END");
	return 0;
}

void EBFV1_elliptic::initialize_MAS(GeomData* pGCData){
	debug_msg("EBFV1_elliptic::initialize_MAS(): START");

	int dim = pGCData->getMeshDim();
	int nedges = 0;
	for(int i=0;i<pGCData->getNumDomains(); i++){
		nedges += pGCData->getNumEdgesPerDomain(i);
	}
	//cout << "nedges: " << nedges << endl;

	pMAS = new MAS;
	pMAS->Eij = new double*[nedges];
	pMAS->Gij = new double*[nedges];
	pMAS->indices = new int*[nedges];
	pMAS->edge_lambda = new double[nedges];

	for(int i=0;i<nedges;i++){
		pMAS->Eij[i] = new double[4*dim];
		pMAS->Gij[i] = new double[4];
		pMAS->indices[i] = new int[2];
	}

	debug_msg("EBFV1_elliptic::initialize_MAS(): END");
}

void EBFV1_elliptic::finalize_MAS(GeomData* pGCData){
	int nedges = 0;
	for(int i=0;i<pGCData->getNumDomains(); i++){
		nedges += pGCData->getNumEdgesPerDomain(i);
	}
	for(int i=0;i<nedges;i++){
		delete[] pMAS->Eij[i]; pMAS->Eij[i] = 0;
		delete[] pMAS->Gij[i]; pMAS->Gij[i] = 0;
		delete[] pMAS->indices[i]; pMAS->indices[i] = 0;
	}
	delete[] pMAS->Eij; pMAS->Eij = 0;
	delete[] pMAS->Gij; pMAS->Gij = 0;
	delete[] pMAS->indices; pMAS->indices = 0;
	delete pMAS;
}

int EBFV1_elliptic::free_matrices(GeomData* pGCData){
	debug_msg("EBFV1_elliptic::free_matrices: START");
	int i, j, k, counter;
	int ndom = pGCData->getNumDomains();
	int dim = pGCData->getMeshDim();
	int idxm[2], idxn[2*dim], pos1, pos2;
	int nrows = 2;
	int ncols = 2*dim;
	double Eij[4*dim];
	double Gij[4];

	counter = 0;
	MatDestroy(&pDStruct->G_freenodes);
	//VecZeroEntries(pDStruct->solution);
	VecZeroEntries(pDStruct->RHS);
	for(i=0; i<ndom; i++){
		MatDestroy(&pDStruct->E_freenodes[i]);
		MatDestroy(&pDStruct->F_freenodes[i]);
		int nedges = pGCData->getNumEdgesPerDomain(i);
		for(j=0; j<nedges; j++){

			idxm[0] = pMAS->indices[counter][0] - 1;
			idxm[1] = pMAS->indices[counter][1] - 1;
			pos1 = dim*idxm[0];
			pos2 = dim*idxm[1];

			for (k=0; k<dim; k++){
				idxn[k] = pos1+k;
				idxn[dim+k] = pos2+k;
			}

			for(k=0; k<4; k++){
				Gij[k] = 0;
			}

			for(k=0; k<4*dim; k++){
				Eij[k] = 0;
			}
			counter++;

			MatSetValues(pDStruct->G,2,idxm,2,idxm,Gij,INSERT_VALUES);
			MatSetValues(pDStruct->E[i],nrows,idxm,ncols,idxn,Eij,INSERT_VALUES);
		}
	}
	MatAssemblyBegin(pDStruct->G,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(pDStruct->G,MAT_FINAL_ASSEMBLY);
	for (int dom=0; dom<ndom; dom++){
		MatAssemblyBegin(pDStruct->E[dom],MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(pDStruct->E[dom],MAT_FINAL_ASSEMBLY);
	}
	debug_msg("EBFV1_elliptic::free_matrices: END");
	return 0;
}
}
