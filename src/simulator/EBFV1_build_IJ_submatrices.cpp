/*
 * EBFV1_build_IJ_submatrices.cpp
 *
 *  Created on: 7 de jun de 2016
 *      Author: rogerio
 *
 *  Build geometric sub matrices: E_ij, G_ij, F_ij
 */

#include "EBFV1_elliptic.h"

namespace PRS{

int EBFV1_elliptic::build_global_matrices(GeomData* pGCData, SimulatorParameters* pSimPar){
		debug_msg("EBFV1_elliptic::build_global_matrices(): START");

		const double* Cij = NULL;
		const int* indices = NULL;
		double *versor = NULL;
		double length;
		int i, nedges, dom, id0, id1, dom_flag, counter, dim, ndom;

		counter = 0;
		dim = pGCData->getMeshDim();
		ndom = pGCData->getNumDomains();

		for (dom=0; dom<ndom; dom++){
			dom_flag = pGCData->domains[dom];
			const double *K = pSimPar->getPermeability(dom+1);
			bool is_K_Isotropic = pSimPar->is_K_Isotropic();
			nedges = pGCData->getNumEdgesPerDomain(dom);
			//cout << "dom: " << dom << " nedges: " << nedges << endl;
			for (i=0; i<nedges; i++){
				pGCData->getVersor(dom,i,&versor);
				pGCData->getCij(dom,i,Cij);
				pGCData->getEdge(dom,i,indices);
				pGCData->getLength(dom,i,length);
				pGCData->getID(dom,indices[0],indices[1],id0,id1);
				divergence_E(Cij,i,dom,dom_flag,indices[2],indices[3],id0,id1,dim,counter,versor,K,is_K_Isotropic);
				divergence_G(Cij,i,dom,dom_flag,indices[2],indices[3],id0,id1,dim,counter,length,versor,K,is_K_Isotropic);
				gradient_F_edges(pDStruct->F[dom],Cij,dom,indices[0],indices[1],indices[2]+1,indices[3]+1,dim,pGCData);
				pMAS->indices[counter][0] = indices[2]+1;
				pMAS->indices[counter][1] = indices[3]+1;
				counter++;
				//cout << dom << "\tCij: " << Cij[0] << " " << Cij[1] << " " << Cij[2] << "\n";
			}
			gradient_F_bdry(pDStruct->F[dom],dom,pGCData,pSimPar);
		}
		//cout << "counter: " << counter << endl;

		Cij = NULL;
		indices = NULL;
		versor = NULL;

		//exit(1);
		debug_msg("EBFV1_elliptic::build_global_matrices(): END");
		return 0;
	}

	int EBFV1_elliptic::divergence_E(const double *Cij, int edge, int dom, int dom_flag, int idx0_global, int idx1_global, int id0, int id1, int dim, int counter, double *versor, const double *K, bool is_K_Isotropic){
		const double I2D[4] = {1.0,.0,.0,1.0};
		const double I3D[9] = {1.0,.0,.0,.0,1.0,.0,.0,.0,1.0};
		int i, j, k;

		double matLij[dim*dim];
		k = 0;
		for (i=0; i<dim; i++){
			for (j=0; j<dim; j++){
				matLij[k++] = versor[i]*versor[j];
			}
		}

		//cout << "Cij: " << Cij[0] << " "  << Cij[1] << " "  << Cij[2] << endl;

		const double *Identity = (dim==2)?I2D:I3D;
		double matSubtrac_ILij[dim*dim];
		for (k=0; k<dim*dim; k++){
			matSubtrac_ILij[k] = -0.5*(Identity[k] - matLij[k]);
		}

		double EA[dim*dim];
		for (i=0; i<dim*dim; i++){
			EA[i] = .0;
		}

		k = 0;
		int pos = 0;
		if ( is_K_Isotropic ){
			int pos1 = 0;
			int pos2 = 0;
			for (i=0; i<dim; i++){
				for (j=0; j<dim; j++){
					EA[k++] = K[pos1]*matSubtrac_ILij[pos2+j];
				}
				pos1 += dim+1;
				pos2 += dim;
			}
		}else{
			for (i=0; i<dim; i++){
				for (j=0; j<dim; j++){
					for (k=0; k<dim; k++)
						EA[pos] += K[dim*i+k]*matSubtrac_ILij[dim*k+j];
					pos++;
				}
			}
		}

		double Eij_part1[3] = {.0,.0,.0};
		for (i=0; i<dim; i++){
			for (j=0; j<dim; j++){
				Eij_part1[i] += Cij[j]*EA[dim*j+i];
			}
		}

		double Eij_row1[2*dim], Eij_row2[2*dim];
		for (i=0; i<dim; i++){
			Eij_row1[i] = Eij_part1[i];
			Eij_row1[i+dim] = Eij_row1[i];
			Eij_row2[i] = -Eij_row1[i];
			Eij_row2[i+dim] = -Eij_row1[i];
		}

		K = 0;
		for (i=0;i<2*dim; i++){
			pMAS->Eij[counter][i] = Eij_row1[i];
			pMAS->Eij[counter][i+2*dim] = Eij_row2[i];
		}

		return 0;
	}

	int EBFV1_elliptic::divergence_G(const double *Cij, int edge, int dom, int dom_flag, int idx0_global, int idx1_global, int id0, int id1, int dim, int counter, double length, double *versor, const double *K, bool is_K_Isotropic){
		int i, j;

		double KL[3] = {.0, .0, .0};
		int pos = 0;
		if ( is_K_Isotropic ){
			for (i=0; i<dim; i++){
				KL[i] = K[pos]*versor[i];
				pos += dim+1;
			}
		}
		else{
			for (i=0; i<dim; i++){
				for (j=0; j<dim; j++){
					KL[i] += K[dim*i+j]*versor[j];
				}
			}
		}

		double aux = .0;
		for (i=0; i<dim; i++){
			aux += Cij[i]*KL[i];
		}
		aux /= -length;

		pMAS->Gij[counter][0] = aux;
		pMAS->Gij[counter][1] = -aux;
		pMAS->Gij[counter][2] = -aux;
		pMAS->Gij[counter][3] = aux;
		return 0;
	}

	int EBFV1_elliptic::gradient_F_bdry(Mat F, int dom, GeomData* pGCData, SimulatorParameters* pSimPar){

		double* Dij; Dij=0;
		double Fij_column1[4], Fij_column2[4], volumeI, volumeJ, volumeK;
		int nedges, nfaces, id0, id1, id2, idx_0, idx_1, idx_2, idx0_global, idx1_global, idx2_global, i, j, pos1, pos2;
		const double C1 = 0.83333333333333333; // 5 divided 6
		const double C2 = 0.16666666666666667; // 1 divided 6


		if (pGCData->getMeshDim()==2){
			nedges = pGCData->getNumBDRYEdgesPerDomain(dom);
			for (j = 0; j<nedges; j++){
				//pGCData->getBdryEdge(dom,j,idx_0,idx_1);
				pGCData->getBdryEdge(dom,j,idx_0,idx_1,idx0_global,idx1_global);
				pGCData->getBdryID(dom,idx_0,idx_1,id0,id1);
				pGCData->getBdryVolume(dom,idx_0,idx_1,volumeI,volumeJ);
				pGCData->getDij(dom,j,&Dij);

				//cout << " dom: " << dom << "\t" << idx0_global+1 << " " << idx1_global+1 << "\tDij: " << Dij[0] << " " << Dij[1] << endl;

				Fij_column1 [0] = C1*Dij[0]/volumeI;
				Fij_column1 [1] = C1*Dij[1]/volumeI;
				Fij_column1 [2] = C2*Dij[0]/volumeJ;
				Fij_column1 [3] = C2*Dij[1]/volumeJ;
				Fij_column2 [0] = C2*Dij[0]/volumeI;
				Fij_column2 [1] = C2*Dij[1]/volumeI;
				Fij_column2 [2] = C1*Dij[0]/volumeJ;
				Fij_column2 [3] = C1*Dij[1]/volumeJ;
				pos1 = 2*idx0_global;
				pos2 = 2*idx1_global;
				int idxm[4] = {pos1,pos1+1,pos2,pos2+1};
				int idxn[2] = {idx0_global,idx1_global};
				MatSetValues(F,4,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
				MatSetValues(F,4,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
			}
		}
		else{
			nfaces = pGCData->getNumBdryFacesPerDomain(dom);
			//cout << "nfaces: " << nfaces << endl;
			// APROXIMACAO DE PRIMEIRA ORDEM PARA O CONTORNO
			if (pSimPar->use_first_order_approximation_on_contour()){

				//cout << "// APROXIMACAO DE PRIMEIRA ORDEM PARA O CONTORNO\n";
				for (j = 0; j<nfaces; j++){
					pGCData->getBdryFace(dom,j,idx_0,idx_1,idx_2,idx0_global,idx1_global,idx2_global);
					pGCData->getBdryID(dom,idx_0,idx_1,idx_2,id0,id1,id2);
					pGCData->getBdryVolume(dom,idx_0,idx_1,idx_2,volumeI,volumeJ,volumeK);
					pGCData->getDij(dom,j,&Dij);
					double Fij_1[3] = {Dij[0]/volumeI, Dij[1]/volumeI, Dij[2]/volumeI};
					double Fij_2[3] = {Dij[0]/volumeJ, Dij[1]/volumeJ, Dij[2]/volumeJ};
					double Fij_3[3] = {Dij[0]/volumeK, Dij[1]/volumeK, Dij[2]/volumeK};
					int pos1 = 3*(id0-1);
					int pos2 = 3*(id1-1);
					int pos3 = 3*(id2-1);
					int idxm_1[3] = {pos1,pos1+1,pos1+2};
					int idxm_2[3] = {pos2,pos2+1,pos2+2};
					int idxm_3[3] = {pos3,pos3+1,pos3+2};
					int idxn[3] = {id0-1, id1-1, id2-1};
					MatSetValues(F,3,idxm_1,1,&idxn[0],Fij_1,ADD_VALUES);
					MatSetValues(F,3,idxm_2,1,&idxn[1],Fij_2,ADD_VALUES);
					MatSetValues(F,3,idxm_3,1,&idxn[2],Fij_3,ADD_VALUES);
				}
			}
			else{// FORMULACAO PARA CONTORNO ORIGINAL
				//cout << "// FORMULACAO PARA CONTORNO ORIGINAL\n";
				for (j = 0; j<nfaces; j++){
					pGCData->getBdryFace(dom,j,idx_0,idx_1,idx_2,idx0_global,idx1_global,idx2_global);
					pGCData->getBdryID(dom,idx_0,idx_1,idx_2,id0,id1,id2);
					pGCData->getBdryVolume(dom,idx_0,idx_1,idx_2,volumeI,volumeJ,volumeK);
					pGCData->getDij(dom,j,&Dij);
					double tmp[3] = {1./(8.*volumeI), 1./(8.*volumeJ), 1./(8.*volumeK)};
					double aux[3][3] = {{6.*tmp[0],tmp[0],tmp[0]},{tmp[1],6.*tmp[1],tmp[1]},{tmp[2],tmp[2],6.*tmp[2]}};
					double Fij_column1[9], Fij_column2[9], Fij_column3[9];
					for (i=0; i<3; i++){
						Fij_column1[3*i] =   aux[i][0]*Dij[0];
						Fij_column1[3*i+1] = aux[i][0]*Dij[1];
						Fij_column1[3*i+2] = aux[i][0]*Dij[2];
						Fij_column2[3*i] = aux[i][1]*Dij[0];
						Fij_column2[3*i+1] = aux[i][1]*Dij[1];
						Fij_column2[3*i+2] = aux[i][1]*Dij[2];
						Fij_column3[3*i] = aux[i][2]*Dij[0];
						Fij_column3[3*i+1] = aux[i][2]*Dij[1];
						Fij_column3[3*i+2] = aux[i][2]*Dij[2];
					}
					int pos1 = 3*(id0-1);
					int pos2 = 3*(id1-1);
					int pos3 = 3*(id2-1);
					int idxm[9] = {pos1,pos1+1,pos1+2, pos2,pos2+1,pos2+2, pos3,pos3+1,pos3+2};
					int idxn[3] = {id0-1, id1-1, id2-1};
					MatSetValues(F,9,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
					MatSetValues(F,9,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
					MatSetValues(F,9,idxm,1,&idxn[2],Fij_column3,ADD_VALUES);
				}
			}
		}

		return 0;
	}

	int EBFV1_elliptic::gradient_F_edges(Mat F, const double *Cij, int dom, int idx0, int idx1, int id0, int id1, int dim, GeomData* pGCData){
		int i;
		double volumeI, volumeJ;

		//cout << "Cij " << id0 << " " << id1 << " " << Cij[0] << " " << Cij[1] <<  " " << Cij[2] << endl;

		pGCData->getVolume(dom,idx0,idx1,volumeI,volumeJ);
		double aux1 = 0.5/volumeI;
		double aux2 = -0.5/volumeJ;
		double Fij_column1[2*dim], Fij_column2[2*dim];
		for (i=0; i<dim; i++){
			Fij_column1[i] = Cij[i]*aux1;
			Fij_column2[i] = Fij_column1[i];
			Fij_column1[i+dim] = Cij[i]*aux2;
			Fij_column2[i+dim] = Fij_column1[i+dim];
		}

		int pos1 = dim*(id0-1);
		int pos2 = dim*(id1-1);
		int idxm[2*dim];				// says on which rows Fij must be assembled into Fg
		for (i=0; i<dim; i++){
			idxm[i] = pos1+i;
			idxm[dim+i] = pos2+i;
		}

		int idxn[2] = {id0-1,id1-1};	// says on which columns Fij must be assembled into Fg
//		if (dom==1){
//			cout << "idxn: " << idxn[0] << "  " << idxn[1] << endl;
//			cout << "volumeI: " << volumeI << "  volumeJ " << volumeJ << endl;
//			cout << "Cij: " << Cij[0] << "  " << Cij[1] << endl;
//			cout << "aux1: " << aux1 << "  aux2:" << aux2 << endl;
//		}

		MatSetValues(F,2*dim,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
		MatSetValues(F,2*dim,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
		return 0;
	}
}
