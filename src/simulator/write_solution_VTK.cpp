
#include "write_solution_VTK.h"

void write_solution_VTK(const char* filename, Physical* pPPData, GeomData* pGCData){
	debug_msg("EBFV1_elliptic::write_solution_VTK: START");
	ofstream fid;
	fid.open(filename);
	//checkFileOpened(fid);

	print_headers(fid,pGCData->getNumNodes());
	print_coordinates(fid,pGCData);
	print_connectivities(fid,pGCData);
	print_celltype(fid,pGCData->getMeshDim(),pGCData->getNumElements());

	// nodal properties:
	print_saturation(fid,pPPData,pGCData->getNumNodes());
	print_pressure(fid,pPPData);
	fid.close();
	debug_msg("EBFV1_elliptic::write_solution_VTK: END");
}

void print_headers(ofstream &fid, int numNodes){
	fid << "# vtk DataFile Version 2.0\n";
	fid << "Two phases flow simulation\n";
	fid << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
	fid << "POINTS " << numNodes << " float\n";
}

void print_coordinates(ofstream &fid, GeomData* pGCData){
	const double *coords = NULL;
	fid << setprecision(8) << scientific;
	for(int i=0; i<pGCData->getNumNodes(); i++){
		pGCData->getCoordinates(i,&coords);
		fid << coords[0] << " " << coords[1] << " " << coords[2] << endl;
	}
}

void print_connectivities(ofstream &fid, GeomData* pGCData){
	const int* connectivities = NULL;
	int dim = pGCData->getMeshDim();
	int numElem = pGCData->getNumElements();
	int size = (dim==2)?3:4;

	fid << "\nCELLS " << numElem << " " << (dim+2)*numElem << endl;
	for (int i=0; i<numElem; i++){
		pGCData->getConnectivities(i,&connectivities);
		fid << size << " ";
		for (int j=0; j<size; j++){
			fid << connectivities[j] << " ";
		}
		connectivities = 0;
		fid << endl;
	}
}

void print_celltype(ofstream& fid, int dim, int numElem){
	fid << "\nCELL_TYPES " << numElem << endl;
	int type = (dim==2)?5:10;				// 5 - triangles;	10 - tetrahedra
	for(int i=0; i<numElem; i++){
		fid << type << endl;
	}
}

void print_saturation(ofstream &fid, Physical* pPPData, int num_nodes){
	double Sw;
	fid << "\nPOINT_DATA "<< num_nodes << endl;
	fid << "SCALARS Saturation float 1\n";
	fid << "LOOKUP_TABLE default\n";
	for (int idx=0; idx<pPPData->getNumNodes(); idx++){
		pPPData->getSaturation(idx,Sw);
		fid << Sw << endl;
	}
}


void print_pressure(ofstream &fid, Physical* pPPData){
	double p;
	fid << "SCALARS Pressure float 1\n";
	fid << "LOOKUP_TABLE default\n";
	fid << setprecision(8) << scientific;
	for (int idx=0; idx<pPPData->getNumNodes(); idx++){
		pPPData->getPressure(idx,p);
		fid << p << endl;
	}
}

//void checkFileOpen(ofstream& fid){
//	if ( !fid.is_open() ){
//		char msg[512];
//		sprintf(msg,"File '%s' could not be opened or it does not exist.\nCheck if directory was typed correctly.\n",filename.c_str());
//		throw Exception(__LINE__,__FILE__,msg);
//	}
//}
