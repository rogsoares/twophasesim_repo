/*
 * gradientSaturation.cpp
 *
 *  Created on: 14/05/2009
 *      Author: rogerio
 */

#include "EBFV1_elliptic.h"

namespace PRS{
	void EBFV1_elliptic::calculatePressureGradient(GeomData *pGCData, Physical *pPPData){
		//CPU_Profile::Start();
		debug_msg("EBFV1_elliptic::calculatePressureGradient: START");

		resetPressureGradient(pGCData,pPPData);
		int ndom = pGCData->getNumDomains();
		for (int dom=0; dom<ndom; dom++){
			calc_p_grad_1(dom,pGCData,pPPData);
			calc_p_grad_2(dom,pGCData,pPPData);
		}
		calc_p_grad_3(pGCData,pPPData);

		//CPU_Profile::End("calculatePressureGradient");
		debug_msg("EBFV1_elliptic::calculatePressureGradient: END");
	}

	void EBFV1_elliptic::calc_p_grad_1(int dom, GeomData *pGCData, Physical *pPPData){
		double p_I, p_J, val;
		const double* Cij = NULL;
		double* p_grad_I = NULL;
		double* p_grad_J = NULL;
		int i, nedges, edge, idx0, idx1, idx0_global, idx1_global, id0, id1, dim;

		dim = pGCData->getMeshDim();
		nedges = pGCData->getNumEdgesPerDomain(dom);
		for (edge=0; edge<nedges; edge++){
			pGCData->getCij(dom,edge,Cij);
			pGCData->getEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
			pGCData->getID(dom,idx0,idx1,id0,id1);
			pPPData->getPressure(idx0_global,p_I);
			pPPData->getPressure(idx1_global,p_J);
			pPPData->get_pw_Grad(dom,idx0,p_grad_I);
			pPPData->get_pw_Grad(dom,idx1,p_grad_J);

			//cout << "\nstep1: edge " << idx0 << "  " << idx1 << " in dom " << dom;
			val = 0.5*(p_I + p_J);
			for (i=0; i<dim; i++){
				p_grad_I[i] += val*Cij[i];
				p_grad_J[i] += -val*Cij[i];
			}
		}
	}

	void EBFV1_elliptic::calc_p_grad_2(int dom, GeomData *pGCData, Physical *pPPData){
		const double* Dij = NULL;
		double* p_grad_I = NULL;
		double* p_grad_J = NULL;
		double* p_grad_K = NULL;

		double p_I, p_J;
		int nedges, edge, idx0, idx1, idx2, idx0_global, idx1_global, idx2_global, dim;
		cout << setprecision(8) << fixed;
		dim = pGCData->getMeshDim();
		if (dim==2){
			nedges = pGCData->getNumBDRYEdgesPerDomain(dom);
			for (edge=0; edge<nedges; edge++){
				pGCData->getDij(dom,edge,Dij);
				pGCData->getBdryEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
				pPPData->getPressure(idx0_global,p_I);
				pPPData->getPressure(idx1_global,p_J);
				pPPData->get_pw_Grad(dom,idx0,p_grad_I);
				pPPData->get_pw_Grad(dom,idx1,p_grad_J);

				//cout << "dom: " << dom << "\tDij: " << Dij[0] << "  " << Dij[1] << endl;

				// 0.16666667 = 1/6 (avoid round-off errors)
				for (int i=0; i<dim; i++){
					p_grad_I[i] += ((5.*p_I + p_J)/6.0)*Dij[i];
					p_grad_J[i] += ((p_I + 5.*p_J)/6.0)*Dij[i];
				}
			}
		}
		else{
			int nfaces = pGCData->getNumBdryFacesPerDomain(dom);
			double line1[3] = {6., 1., 1.};
			double line2[3] = {1., 6., 1.};
			double line3[3] = {1., 1., 6.};
			double dot1, dot2, dot3, p_vec[3], aux[3];
			for (int face=0; face<nfaces; face++){
				pGCData->getDij(dom,face,Dij);
				pGCData->getBdryFace(dom,face,idx0,idx1,idx2,idx0_global,idx1_global,idx2_global);
				pPPData->getPressure(idx0_global,p_vec[0]);
				pPPData->getPressure(idx1_global,p_vec[1]);
				pPPData->getPressure(idx2_global,p_vec[2]);
				pPPData->get_pw_Grad(dom,idx0,p_grad_I);
				pPPData->get_pw_Grad(dom,idx1,p_grad_J);
				pPPData->get_pw_Grad(dom,idx2,p_grad_K);

				dot1 = dot2 = dot3 = .0;
				for (int i=0; i<3; i++){
					dot1 += line1[i]*p_vec[i];
					dot2 += line2[i]*p_vec[i];
					dot3 += line3[i]*p_vec[i];
					aux[i] = Dij[i]/8;		// aux = Dij / 8 (avoid round-off error)
				}

				for (int i=0; i<3; i++){
					p_grad_I[i] += aux[i]*dot1;
					p_grad_J[i] += aux[i]*dot2;
					p_grad_K[i] += aux[i]*dot3;
				}
			}
		}
	}

	void EBFV1_elliptic::calc_p_grad_3(GeomData *pGCData, Physical *pPPData){
		double* p_grad = NULL;
		double vol;//, p;
		int i, dom, ndom, nnodes, node, dim;

		dim = pGCData->getMeshDim();
		ndom = pGCData->getNumDomains();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pPPData->get_pw_Grad(dom,node,p_grad);
				pGCData->getVolume(dom,node,vol);

//				cout << setprecision(8) << fixed;
//				cout << "dom[" << dom <<  "], node: " << node << "\tvolume: " << vol << " pgrad: ";
				for (i=0; i<dim; i++){
					p_grad[i] /= vol;
					//cout << p_grad[i] << " ";
				}
//				pPPData->getPressure(node,p);
//				cout << " pressao: " << p << endl;
			}
		}
//		exit(1);
	}

	void EBFV1_elliptic::resetPressureGradient(GeomData *pGCData, Physical *pPPData){
		double* p_grad = NULL;
		int dom, ndom, nnodes, node, i, dim;

		dim = pGCData->getMeshDim();
		ndom = pGCData->getNumDomains();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pPPData->get_pw_Grad(dom,node,p_grad);
				for (i=0; i<dim; i++){
					p_grad[i] = .0;
				}
			}
		}
	}
}
