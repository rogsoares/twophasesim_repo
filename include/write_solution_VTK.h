#ifndef WRITE_SOLUTION_VTK_H_
#define WRITE_SOLUTION_VTK_H_

#include "Physical.h"
#include "GeomData.h"
using namespace PRS;


void write_solution_VTK(const char* filename, Physical* pPPData, GeomData* pGCData);
void print_headers(ofstream &fid, int);
void print_coordinates(ofstream&, GeomData*);
void print_connectivities(ofstream&, GeomData*);
void print_celltype(ofstream&,int,int);
void printCellTypeList(ofstream &, int, int);

void print_saturation(ofstream &fid, Physical* pPPData, int num_nodes);
void print_pressure(ofstream&, Physical*);

#endif
