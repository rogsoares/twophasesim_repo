/*
 * EBFV1-hyperbolic.h
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#ifndef EBFV1HYPERBOLI_H_
#define EBFV1HYPERBOLI_H_

#include "OilProductionManagement.h"
#include "Physical.h"
#include "MeshData.h"

namespace PRS{

	class EBFV1_hyperbolic{
	public:

		EBFV1_hyperbolic();
		~EBFV1_hyperbolic();
		double solver(MeshData *pMData, GeomData *pGDCata, SimulatorParameters *pSimPar, Physical *pPPData, double &timestep);

	protected:

		// Main functions related to the advective equation presented in order that they must be called.
		double calculateVelocityField(GeomData *pGDCata, Physical *pPPData, SimulatorParameters *pSimPar, int dom, int dim);
		double calculateIntegralAdvectiveTerm(GeomData *pGDCata, SimulatorParameters *pSimPar, Physical *pPPData, int dom, int dim, double &timeStep);
		void calculateExplicitAdvanceInTime(GeomData *pGCData, SimulatorParameters *pSimPar, Physical *pPPData, double delta_T);

		void nodeWithOut_Wells(GeomData *pGCData, SimulatorParameters *pSimPar, Physical *pPPData, double delta_T);
		void nodeWith_Wells(GeomData *pGCData, SimulatorParameters *pSimPar, Physical *pPPData, double delta_T);
		
		/// For mesh adaptation, saturation field interpolation between old and new mesh over Sw_t, and not Sw_t+1
		void saveSwField();

		void setRecoveredOil(double val){
			oilRecovered = val;
		}

		double getCumulativeOil() const{
			return _cumulativeOil;
		}

		void setCumulativeOil(double co){
			_cumulativeOil = co;
		}

		double getRecoveredOil() const{
			return oilRecovered;
		}

		double oilRecovered;

		// Numeric parameter to compute nonvisc term
		double alpha_max;

		// Saturation gradient is calculated for all domains at once. Nodes on boundary domains contains one gradient vector for each domain.
		void calculateSaturationGradient(GeomData *pGDCata, SimulatorParameters *pSimPar, Physical *pPPData, int dim);

	private:
		// every new time-step, nodal gradient must be set to zero and start a new calculation
		void resetSaturationGradient(GeomData *pGDCata, Physical *pPPData);
		
		// loop over all edges (omega domain
		void calc_Sw_grad_1(GeomData *pGDCata, Physical *pPPData, int dom, int dim);
		
		// loop over all boundary edges (2-D) or all external faces (3-D)
		void calc_Sw_grad_2(GeomData *pGDCata, Physical *pPPData, int dom, int dim);
		
		// Averaging by Total Volume in 3-D (area in 2-D problems)
		void calc_Sw_grad_3(GeomData *pGCData, Physical *pPPData, int dim);

		// Imposition of Homogeneus Neumman Boundary Conditions
		void calc_Sw_grad_4(GeomData *pGDCata, Physical *pPPData, SimulatorParameters *pSimPar, int dim);
		
		// If restart is required, then oil production file data is used to update the cumulated oil value
		void setCumulativeOilProd();

		double _cumulativeOil;

		PetscErrorCode ierr;

		bool PRINT_DEBUG;

		ofstream fid_PRINT_DEBUG;
	};
}


#endif /* EBFV1HYPERBOLI_H_ */
