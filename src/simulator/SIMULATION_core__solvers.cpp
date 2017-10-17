/*
 * SIMULATION_core__solvers.cpp
 *
 *  Created on: 27/08/2012
 *      Author: rogsoares
 */

#include "SIMULATION_core.h"

namespace PRS{

	int SIMULATION_core::start_simulation(){

		double t1 = MPI_Wtime();
		char filename[512];
		double timeStep;

		if (simFlag==STEADY_STATE){
			PetscPrintf(PETSC_COMM_WORLD,"\nStart simulation: Steady State\n\n");
			pElliptic_eq->solver(pMData,pGCData,pSimPar,pPPData);
			sprintf(filename,"steady-state.vtk");
			write_solution_VTK(filename,pPPData,pGCData);
		}
		else if (simFlag==TRANSIENT){
			PetscPrintf(PETSC_COMM_WORLD,"\nStart simulation: Transient\n\n");
			cout << setprecision(8);
			while ( !pSimPar->finishSimulation() ){
				pElliptic_eq->solver(pMData,pGCData,pSimPar,pPPData);
				pHyperbolic_eq->solver(pMData,pGCData,pSimPar,pPPData,timeStep);

//				if (pSimPar->allowPrinting_VTK()){
//					sprintf(filename,"%s-%d.vtk",pSimPar->expofName.c_str(),pSimPar->getStepOutputFile());
//					write_solution_VTK(filename,pPPData,pGCData);
//					pSimPar->incrementeStepOutputFile();
//					pSimPar->setNotAllowPrinting_VTK();
//				}
			}
		}
		PetscPrintf(PETSC_COMM_WORLD,"\nFinished!\n");

		double t2 = MPI_Wtime();
		double h,m,s;
		convertSecToTime(t2-t1,&h,&m,&s);
		cout << setprecision(0) << fixed << "\n\nCPU time elapsed: " << h << "h " << m << "m " << s << "s\n\n";
		return 0;
	}

}
