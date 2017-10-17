#ifndef _IMPES_FORMULATION_H_
#define _IMPES_FORMULATION_H_

#include "EBFV1__pre-processors.h"
#include "solverEquation_methods.h"
#include "write_solution_VTK.h"


//!PRS: Petroleum Reservoir Simulator
namespace PRS{

void initializeParameters();

	class SIMULATION_core{
	public:

		SIMULATION_core();
		~SIMULATION_core();

		int initialize(int argc, char **argv);
		int start_simulation();
		int finalize();

	private:

		/// pointer to solver pressure and saturation fields
		EBFV1_elliptic* pElliptic_eq;
		EBFV1_hyperbolic* pHyperbolic_eq;
		Physical *pPPData;
		SimulatorParameters *pSimPar;
		GeomData *pGCData;
		MeshData *pMData;
		OilProductionManagement *pOilProduction;

		int simFlag;
		enum SIMULATION_States{STEADY_STATE, TRANSIENT, MIMPES_ADAPT};
	};
}
#endif
