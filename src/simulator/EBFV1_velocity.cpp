/*
 * velocityField.cpp
 *
 *  Created on: 14/01/2009
 *      Author: rogerio
 *
 *
 *
 *      vel = -lambda*K*grad_p
 *
 *      vel = vel_N + vel_P
 *
 *      vel_N =
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{

	double EBFV1_hyperbolic::calculateVelocityField(GeomData *pGCData, Physical *pPPData, SimulatorParameters *pSimPar, int dom, int dim){
		debug_msg("EBFV1_hyperbolic::calculateVelocityField(): START");

		const double* Cij = NULL;
		const double* pGrad_I = NULL;
		const double* pGrad_J = NULL;
		double* versor = NULL;

		double p_I, p_J, MobIJ, length, dp_dt;
		double dot, vel[3], pw_grad_IJ[3], unit_vector[3], VIJN[3], VIJFD[3], aux[3];
		int i, j, idx0_global, idx1_global, idx0, idx1, id0, id1;

		int nedges = pGCData->getNumEdgesPerDomain(dom);
		for(int edge=0; edge<nedges; edge++){
			pGCData->getCij(dom,edge,Cij);
			pGCData->getEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
			pGCData->getID(dom,idx0,idx1,id0,id1);
			pPPData->getPressure(idx0_global,p_I);
			pPPData->getPressure(idx1_global,p_J);
			pPPData->get_pw_Grad_const(dom,idx0,pGrad_I);
			pPPData->get_pw_Grad_const(dom,idx1,pGrad_J);
			pGCData->getLength(dom,edge,length);
			pGCData->getVersor(dom,edge,&versor);

			for (i=0; i<dim; i++){
				unit_vector[i] = -versor[i];
			}

			// average grad_p
			// ------------------------------------------------------------------
			for (i=0; i<dim; i++){
				pw_grad_IJ[i] = .5*(pGrad_I[i] + pGrad_J[i]);
			}

			// orthogonal velocity component to edge
			// ------------------------------------------------------------------
			dot = .0;
			for (i=0; i<dim; i++){
				dot += unit_vector[i]*pw_grad_IJ[i];
			}
			for (i=0; i<dim; i++){
				VIJN[i] = pw_grad_IJ[i] - dot*unit_vector[i];
			}

			// parallel velocity component to edge
			// ------------------------------------------------------------------
			dp_dt = (p_J - p_I)/length; // (Finite difference approximation)
			for (i=0; i<dim; i++){
				VIJFD[i] = unit_vector[i]*dp_dt;
			}

			// total velocity
			// ------------------------------------------------------------------
			for (i=0; i<dim; i++){
				aux[i] = VIJFD[i]+VIJN[i];
			}

			// MIMPES
//			pPPData->getVelocity_new(dom,edge,vel);
//			pPPData->setVelocity_old(dom,edge,vel);

			const double *K = pSimPar->getPermeability(pGCData->domains[dom]);
			for (i=0; i<dim; i++){

				// initialize velocity
				vel[i] = 0;
				for (j=0; j<dim; j++){
					vel[i] += K[dim*i+j]*aux[j];
				}
			}

			pPPData->getAverageMobility(idx0_global,idx1_global,MobIJ);
			for (i=0; i<dim; i++){
				vel[i] *= -MobIJ;
			}
			pPPData->setVelocity_new(dom,edge,vel);

			cout << setprecision(8);
//			cout << "dom: " << dom << " edge: " << id0 << "," << id1 << "  vel:" << vel[0] << "  " << vel[1] << "  pGrad:";
//			cout << pw_grad_IJ[0] << "  " << pw_grad_IJ[1] << "  versor:";
//			cout << unit_vector[0] << "  " << unit_vector[1] << "  length:";
//			cout << length << "  MobIJ:" << MobIJ << "  " << endl;
		}

		debug_msg("EBFV1_hyperbolic::calculateVelocityField(): END");
//		exit(1);
		return 0;
	}
}

