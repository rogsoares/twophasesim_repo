/*
 * EBFV1-hyperbolic.cpp
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{

	EBFV1_hyperbolic::EBFV1_hyperbolic(){
		_cumulativeOil = .0;
				setCumulativeOilProd();
	}

	EBFV1_hyperbolic::~EBFV1_hyperbolic(){
	}

	double EBFV1_hyperbolic::solver(MeshData *pMData, GeomData *pGCData, SimulatorParameters *pSimPar, Physical *pPPData, double &timeStep){
		debug_msg("EBFV1_hyperbolic::solver(): START");

		static int timestep_counter = 0;			// counts number of time steps every new VTK
		int dim = pGCData->getMeshDim();

		// initialize time step with a very high number
		timeStep = 1.0e+10;
		int ndom = pGCData->getNumDomains();
		for (int dom=0; dom<ndom; dom++){
			calculateVelocityField(pGCData,pPPData,pSimPar,dom,dim);
		}
		//exit(1);

		// calculate saturation gradient if adaptation or high order approximation were required
		if ( pSimPar->useHOApproximation()){
			calculateSaturationGradient(pGCData,pSimPar,pPPData,dim);
		}

		pPPData->resetNonvisc(alpha_max);
		for (int dom=0; dom<ndom; dom++){
			calculateIntegralAdvectiveTerm(pGCData,pSimPar,pPPData,dom,dim,timeStep);
		}

		pSimPar->correctTimeStep(timeStep);											// correct time-step value to print out the desired simulation moment
		pSimPar->saveCurrentSimulationTimes();										// it must be called before cumulative simulation time
		pSimPar->setCumulativeSimulationTime(timeStep); 							// AccSimTime = AccSimTime + timeStep
		calculateExplicitAdvanceInTime(pGCData,pSimPar,pPPData,timeStep);			// Calculate saturation field: Sw(n+1)

		timestep_counter++;
		// oil production output
//		if (pSimPar->timeToPrintVTK()){
//			pOPManager->printOilProduction(timeStep,pSimPar->getCumulativeSimulationTime(),pSimPar->getSimTime(),getRecoveredOil(),getCumulativeOil(),timestep_counter);
//			timestep_counter = 0;
//		}

		//pSimPar->printOutVTK(pPPData,pSimPar,pGCData,exportSolutionToVTK);

		debug_msg("EBFV1_hyperbolic::solver(): END");
		return 0;
	}

	void EBFV1_hyperbolic::setCumulativeOilProd(){
//		// set cumulative oil production
//		if (pSimPar->useRestart()){
//			string expofn;
//			pSimPar->getExportFileName(expofn);
//			char str[512]; sprintf(str,"%s_oil-production.csv",expofn.c_str());
//			ifstream fid;
//			fid.open(str);
//
//			string strline;
//			for (int i=0; i<=pSimPar->getLastPVI(); i++){
//				getline(fid,strline);
//				//cout << "i = " << i << "  " << strline << endl;
//			}
//			string data[5];
//			fid >> data[0] >> data[1] >> data[2] >> data[3] >> _cumulativeOil >> data[4];
//			cout << "\t" << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << _cumulativeOil << "\t" << data[4];
//			_cumulativeOil *= pOPManager->getInitialOilVolume();
//			cout << "\n_cumulativeOil = " << _cumulativeOil << endl;
//			cout << "pSimPar->getLastPVI() = " << pSimPar->getLastPVI() << endl;
//			fid.close();
//		}
	}
}
