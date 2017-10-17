#include "EBFV1_hyperbolic.h"

namespace PRS{

	void EBFV1_hyperbolic::calculateExplicitAdvanceInTime(GeomData *pGCData, SimulatorParameters *pSimPar, Physical *pPPData, double delta_T){
		debug_msg("EBFV1_hyperbolic::calculateExplicitAdvanceInTime(): START");

		//saveSwField();					// save Sw field (Sw_t) before calculate Sw_t+1
		nodeWithOut_Wells(pGCData,pSimPar,pPPData,delta_T);		// Advance time for all free node not located in wells
		nodeWith_Wells(pGCData,pSimPar,pPPData,delta_T);		// Advance time for nodes in production well (fractional-steps)
	//	exit(1);
		debug_msg("EBFV1_hyperbolic::calculateExplicitAdvanceInTime(): END");
	}

	void EBFV1_hyperbolic::nodeWithOut_Wells(GeomData *pGCData, SimulatorParameters *pSimPar, Physical *pPPData, double delta_T){
		double Sw,Sw_old,nonvisc,volume;
		int idx;
		int nnode = pPPData->getNumNodeFreeWells();
		double phi = 0.2;//pSimPar->getPorosity(1);

		for(int i=0; i<nnode; i++){
			pPPData->getFreeIndex(i,idx);
			pPPData->getSaturation(idx,Sw_old);
			pPPData->getNonvisc(idx,nonvisc);
			pGCData->getVolume(idx,volume);

	//		cout << setprecision(8) << "idx: " << idx << " nonvisc: " << nonvisc << " volume: " << volume << endl;

			// alterar para usar valor de porosidade lido de arquivo
			Sw = Sw_old - (delta_T/(phi*volume))*nonvisc;

			if ( fabs(Sw)<1e-8 ){
				 Sw = .0;
			}
			pPPData->setSaturation(idx,Sw);

			if (Sw > 1 || Sw < .0){
				char msg[256]; sprintf(msg,"Sw[%d/%d] = %.8e. 0.0 <= Sw <= 1.0 \n",i,nnode,Sw);
				throw Exception(__LINE__,__FILE__,msg);
			}
		}
	}

	void EBFV1_hyperbolic::nodeWith_Wells(GeomData *pGCData, SimulatorParameters *pSimPar, Physical *pPPData, double delta_T){
		const int N = pSimPar->getWellTimeDiscretion();		// get number of delta subdivisions
		double dt_well = (double)delta_T/1000;					// time step for well nodes saturation

		int i,j, well_idx;
		double Sw,Sw0,Sw_old,Vt,Qi,Vi,Qwi,Qt,fw,fo,nonvisc,wp,cml_oil,Qo,Qw;

		// TODO: nao usar constantes para identificar pocos!
		Qt = pSimPar->getFlowrateValue(51);					// source/sink term
		Vt = pSimPar->getWellVolume(51);					// for node i, Qi is a fraction of total well flow rate (Qt)

		cout << setprecision(8);

		cml_oil = .0;
		Qo = .0;
		Qw = .0;
		double phi = 0.2;//pSimPar->getPorosity(1);
		int nnodes = pPPData->getNumNodesWells();

		for (i=0; i<nnodes; i++){
			pPPData->getNeumannIndex(i,well_idx);
			pPPData->getNonvisc(well_idx,nonvisc);
			pPPData->getSaturation(well_idx,Sw_old);
			pGCData->getVolume(well_idx,Vi);
			Qi = Qt*(Vi/Vt);								// Fluid (water+oil) flow rate through node i
			wp = phi*Vi;

			for (j=0; j<N; j++){
				Sw0 = Sw_old - (dt_well/wp)*(nonvisc);
				fw = pPPData->getFractionalFlux(Sw0);
				Qwi = fabs(fw*Qi);
				Sw = Sw0 - dt_well*(Qwi/wp);
				Sw_old = Sw;
			}

			//cout << "Qt: " << Qt << " Vt: " << Vt << " nnodes: " << nnodes << "  Qi: " << Qi << "  wp: " << wp << " Sw: " << Sw << endl;

			if (Sw > 1.01 || Sw < .0){
				char msg[256]; sprintf(msg,"Water saturation is out-of-bound [0 1]. Sw = %.4f\n",Sw);
			}

			pPPData->setSaturation(well_idx,Sw);
			fw = pPPData->getFractionalFlux(Sw);	    // oil fractional flux
			fo = pPPData->getOilFractionalFlux(Sw);		// oil fractional flux
			Qo += fabs(Qi*fo);
			Qw += fabs(Qi*fw);
			cml_oil += Qo;
		}

		setRecoveredOil(Qo/(Qo+Qw));
		cml_oil = cml_oil*delta_T + getCumulativeOil();
		setCumulativeOil(cml_oil);
	}

	void EBFV1_hyperbolic::saveSwField(){
//		double Sw;
//		int i, nnodes;
//		pGCData->getMeshNodes(nnodes);
//		for(i=0; i<nnodes; i++){
//			pPPData->getSaturation(i,Sw);
//			pPPData->setSw_old(i,Sw);
//		}
	}
}
