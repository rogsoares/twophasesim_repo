/*
 * PMA_read.cpp
 *
 *  Created on: 29 de fev de 2016
 *      Author: rogerio
 */


#include "ParMeshAdapt.h"

void ParMeshAdapt::read(const char* filename){

	double t1 = MPI_Wtime();

	cout << "\nReading mesh from file:\n";
	cout << "-------------------------------------------------\n";
	// read data from file
	ifstream fid;
	fid.open(filename);
	char msg[512]; sprintf(msg,"ERROR:\nFile could not be opened or it does not exist.\n"
			                   "Check if file name and/or path was typed correctly.\n"
			                   "File: <%s>\n",filename);

	if ( !fid.is_open() ){
		cerr << msg;
		exit(1);
	}

	char line[256];
	fid.getline (line,256);

	if ( !strcmp("$NOD",line) ){
		cerr << "\nERROR: Cannot read Gmsh v.1 file format! Exiting....\n\n";
		exit(1);
	}

	fid.getline (line,256);
	fid.getline (line,256);
	fid.getline (line,256);

	cout << "Reading nodes...    ";
	int ID=0, numElements;
	double x, y, z;
	fid >> nVertices;
	for(int i=0;i<nVertices;i++){
		fid >> ID;
		fid >> x >> y >> z;
		createVertex(ID,(double)x,y,z,-1,-1);
	}
	cout << " Done!\n";

	fid.getline (line,256);
	fid.getline (line,256);
	fid.getline (line,256);

	fid >> numElements;
	int flags[10];
	int elemType,physical,iElm,geom, id1,id2,numFlags;
//	int eType;
	nElem = 0;
	nBdryEdges = 0;
	nBdryFaces = 0;

	cout << "Reading elements... ";
	for (int i=0; i<numElements; i++){
		fid >> iElm >> elemType >> numFlags;

		for(int i=0;i<numFlags;i++){
			fid >> flags[i];
		}

		physical = flags[0];
		geom = flags[1];

		switch(elemType){
		// precribed nodes
		case 15:
			fid >> id1;
			setVertex(id1,physical,geom);
			break;
			// prescribed edges
		case 1:
			fid >> id1 >> id2;

			if (id1>id2){
				swap(id1,id2);
			}

			EdgeFromFile edge;
			edge.id1 = id1;
			edge.id2 = id2;
			edge.physical = physical;
			edge.geom = geom;
			tmpEdgeList.push_back(edge);
			nBdryEdges++;
			break;
			// prescribed triangles
		case 2:

			// Create a triangle element:
			ElementFromFile tri;
			fid >> tri.id1 >> tri.id2 >> tri.id3;
			tri.physical = physical;
			tri.geom = geom;
			tri.etype = TRIANGLE;
			etype = TRIANGLE;
			tmpTrianglesList.push_back(tri);
			nBdryFaces++;
			break;
		case 3:
			break;
		case 4:
			ElementFromFile elm;
			fid >> elm.id1 >> elm.id2 >> elm.id3 >> elm.id4;
			elm.physical = physical;
			elm.geom = geom;
			elm.etype = TETRA;
			etype = TETRA;
			tmpTetraList.push_back(elm);
			break;
		}
	}
	cout << " Done!\n";
	double t2 = MPI_Wtime();

	double frac,h,m,s;
	double t = t2-t1;
	frac = modf(t/3600.,&h);
	frac = modf(frac*60.,&m);
	frac = modf(frac*60.,&s);
	cout << setprecision(0) << fixed << "Elapsed Time: " << h << "h " << m << "min " << s << "s\n";

	_dim = 2;
	if (tmpTetraList.size()){
		_dim = 3;
	}

//	cout << "\nMesh statistics\n";
//	cout << " 1) Number of nodes: " << nVertices << endl;
//	cout << " 2) Boundary edges : " << tmpEdgeList.size() << endl;
//	cout << " 3) Triangles      : " << tmpTrianglesList.size() << endl;
//	cout << " 4) Tetrahedral    : " << tmpTetraList.size() << endl;
//	cout << "-------------------------------------------------\n\n\n";
}

void ParMeshAdapt::write(const char* filename){

	ofstream fid;
	fid.open(filename);

	fid << "# vtk DataFile Version 2.0\n";
	fid << "PMA: Parallel Mesh Adaptation\n";
	fid << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
	fid << "POINTS " << nVertices << " float\n";

	// print vertices coordinates
	// ------------------------------------------------------------------------
	VertexInfo* vertex; vertex=0;
	TVertexDBIter node_it = VertexDB.begin();
	for( ; node_it!=VertexDB.end(); node_it++){
		vertex = node_it->second;
		if (vertex->toCompMesh){
			fid << vertex->coords[0] << " " << vertex->coords[1] << " " << vertex->coords[2] << endl;
		}
	}

	// print element connectivities
	// ------------------------------------------------------------------------
	Element *elm;
	int numElem, dim = 2;

	getNumElementsOfComputationalMesh(numElem);
	//cout << "NumElementsOfComputationalMesh: " << numElem << endl;

	if (!numElem){
		cerr << "Error: no elements found to print! Exiting....\n";
		exit(1);
	}

	fid << "\nCELLS " << numElem << " " << (dim+2)*numElem << endl;
	for ( IdxIter it = idxlist.begin(); it!=idxlist.end(); it++){
		getElement(*it,&elm);
		fid << 3 << " " << elm->id1-1 << " " << elm->id2-1 << " " << elm->id3-1 << "\n";
	}

	fid << "\nCELL_TYPES " << numElem << endl;
	int type = (dim==2)?5:10;				// 5 - triangles;	10 - tetrahedra
	for(int i=0; i<numElem; i++){
		fid << type << endl;
	}

	fid << "\nCELL_DATA " << numElem << endl;
	fid << "SCALARS RefUnrefLevel int 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	for ( IdxIter it = idxlist.begin(); it!=idxlist.end(); it++){
		getElement(*it,&elm);
		fid << elm->numrefine << "\n";
	}
	fid.close();
}

void ParMeshAdapt::write(){

	static int counter = 1;
	char fname[256];
	sprintf(fname,"%s_%d.vtk",output_filename,counter);
	counter++;

	ofstream fid;
	fid.open(fname);

	fid << "# vtk DataFile Version 2.0\n";
	fid << "PMA: Parallel Mesh Adaptation\n";
	fid << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
	fid << "POINTS " << nVertices << " float\n";

	// print vertices coordinates
	// ------------------------------------------------------------------------
	VertexInfo* vertex; vertex=0;
	for( TVertexDBIter node_it = VertexDB.begin(); node_it!=VertexDB.end(); node_it++){
		vertex = node_it->second;
		if (vertex->toCompMesh){
			fid << vertex->coords[0] << " " << vertex->coords[1] << " " << vertex->coords[2] << endl;
		}
	}


	Element *elm;
	int numElem, dim = 2;
	getNumElementsOfComputationalMesh(numElem);

	if (!numElem){
		cerr << "Error: no elements found to print! Exiting....\n";
		exit(1);
	}

	// print element connectivities
	fid << "\nCELLS " << numElem << " " << (dim+2)*numElem << endl;
	for ( IdxIter it = idxlist.begin(); it!=idxlist.end(); it++){
		IDX_str *idx = &(*it);
		Element *elm = &this->pBase_EDS[idx->base].pLeaves[idx->level][idx->leaf];

		fid << 3 << " ";

		vertex = VertexDB[elm->id1];
		fid << vertex->ID_mapped << " ";
		vertex = VertexDB[elm->id2];
		fid << vertex->ID_mapped << " ";
		vertex = VertexDB[elm->id3];
		fid << vertex->ID_mapped << " ";
		fid << endl;
	}

	fid << "\nCELL_TYPES " << numElem << endl;
	int type = (dim==2)?5:10;				// 5 - triangles;	10 - tetrahedra
	for(int i=0; i<numElem; i++){
		fid << type << endl;
	}

	fid << "\nCELL_DATA " << numElem << endl;
	fid << "SCALARS RefUnrefLevel int 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	for ( IdxIter it = idxlist.begin(); it!=idxlist.end(); it++){
		getElement(*it,&elm);
		fid << elm->numrefine << "\n";
	}

	fid.close();
}

void ParMeshAdapt::statistics(){

	int nedges=0;
	getNumEdges(nedges);

	cout << "\n";
	cout << "Mesh statistics\n";
	cout << "-------------------------------------------------\n";
	cout << "1) nodes      : " << nVertices << endl;
	cout << "2) edges      : " << nedges << endl;
	if (_dim==2){
		cout << "3) bdry edges : " << nBdryEdges << endl;
		cout << "4) triangles  : " << idxlist.size() << endl;
	}
	else{
		cout << "3) bdry edges : " << nBdryEdges << endl;
		cout << "3) bdry faces : " << nBdryFaces << endl;
		cout << "4) tetrahedra : " << idxlist.size() << endl;
	}

	cout << "\n";

}
