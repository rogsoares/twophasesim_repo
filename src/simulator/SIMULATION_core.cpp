#include "SIMULATION_core.h"

namespace PRS {

	SIMULATION_core::SIMULATION_core(){
		pElliptic_eq = 0;
		pHyperbolic_eq = 0;
	}

	SIMULATION_core::~SIMULATION_core(){
	}

	int SIMULATION_core::initialize(int argc, char **argv){
		debug_msg("SIMULATION_core::initialize(): START");

		PetscInitialize(&argc,&argv,(char *)0,(char *)0);

		if ( argc<2 ){
			char msg[256];
			sprintf(msg,"You MUST type: ./PADAMEC_AMR.exe a b, where:\na = 0 or 1 (Steady State or Transient\n");
			throw Exception(__LINE__,__FILE__, msg );
		}

		simFlag = atoi(argv[1]);
		// printSimulationHeader();

		// Initialize simulation pointers
		pPPData = new Physical;
		pGCData = new GeomData;
		pSimPar = new SimulatorParameters;
		pMData = new MeshData;

		ParMeshAdapt* pMesh = new ParMeshAdapt;

		// Load data from files
		pSimPar->setCommandLineArguments(argv,argc);
		pSimPar->defineExactSolution();
		pSimPar->getParametersDirectoryPath();
		pSimPar->loadNumericParameters();
		pGCData->read_geometry(pSimPar->parametersPath.c_str());

		// read mesh from file and create edge data structure
		pMesh->read(pSimPar->prepfName.c_str());
		pMesh->initialize();
		pMesh->statistics();

		// Load data from files
		pSimPar->loadPhysicalParameters( pMesh->dim() );
		pSimPar->loadSolverParameters();
		pSimPar->loadHighOrderParameter();

		// start pre-processor
		//pp_3D_controlSurface(pMesh,pGCData);
		pre_processor(pMesh,pGCData);

		// Initialize Geometric parameters: finite volume coefficients
		pGCData->initialize(pMesh);

		pSimPar->getWells(pMesh);
		pSimPar->WellFlowRate_weighted(pMesh,pGCData);
		if (pSimPar->getEllipticSolver()==1){
			pSimPar->setInitialOilVolume(pGCData);
		}

		// initialize physical properties
		pPPData->initialize(pMData,pSimPar,pGCData,pMesh);

		// initialize auxiliary vectors to matrix assembling
		pMData->initialize(pMesh,pSimPar,pGCData);

		pGCData->cleanData(pMesh);

		// pMesh is not necessary anymore
		pMesh->cleanupMemory();

		// Oil production output
		//pOilProduction = new OilProductionManagement(pSimPar->getOutputPathName(),pSimPar->getInitialOilVolume(),pSimPar->getTotalInjectionFlowRate(),pSimPar->useRestart());

		// Initialize elliptic and hyperbolic solver pointers
		pElliptic_eq = new EBFV1_elliptic;
		pHyperbolic_eq = new EBFV1_hyperbolic;

		debug_msg("SIMULATION_core::initialize(): END");

	//	exit(1);
		return 0;
	}

	int SIMULATION_core::finalize(){
		// Write to file oil production output. Only rank 0 is in charge of it.
		string path = pSimPar->getOutputPathName();

		// free memory
		delete pElliptic_eq;
		delete pHyperbolic_eq;
		delete pPPData;
		delete pSimPar;
		delete pGCData;
		delete pMData;
		return 0;
	}
}
