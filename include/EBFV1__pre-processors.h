#ifndef EBFV1PREPROCESSORS_H
#define EBFV1PREPROCESSORS_H

#include "ParMeshAdapt.h"
#include "GeomData.h"

using namespace PRS;

void pre_processor(ParMeshAdapt* pMesh, GeomData *pGCData);
void pp_2D(ParMeshAdapt* pMesh, GeomData *pGCData);
void pp_3D(ParMeshAdapt* pMesh, GeomData *pGCData);
void pp_3D_controlSurface(ParMeshAdapt* pMesh, GeomData *pGCData);


void initCoefficients_2D(EdgeInfo** edge);
void calcutale_Cij(double* I, double* J, double* pCentroid, double* Cij, int i);
void init_volumes(ParMeshAdapt* pMesh, int n);
void calculate_triagles_centroids_and_areas(ParMeshAdapt* pMesh);
void getMeshDomains(ParMeshAdapt* pMesh, GeomData* pGCData);
void validate_coefficients(ParMeshAdapt* pMesh, GeomData* pGCData);
void validate_3D(ParMeshAdapt* pMesh, GeomData* pGCData);
void print_coefficients(ParMeshAdapt* pMesh);

void calculate_face_centroids(ParMeshAdapt* pMesh);
void calculate_centroids(ParMeshAdapt* pMesh);
void calculate_volume(ParMeshAdapt* pMesh);
void initCoefficients_3D(ParMeshAdapt* pMesh);
void sort_tetra(int& id1,int& id2,int& id3,int& id4);
void compute_Cij(ParMeshAdapt* pMesh, int id0, int id1, int physical, const double* I_coords, const double* J_coords, const double* face1_centroid, const double* face2_centroid, const double* tetra_centroid);
void calculate_Cij(ParMeshAdapt* pMesh);
void calculate_Dij(ParMeshAdapt* pMesh, GeomData *pGCData);
void calculate_NumFacesPerSubDomain(ParMeshAdapt* pMesh, GeomData *pGCData);
void cleanup_auxiliary_data(ParMeshAdapt* pMesh);


// construct control surfaces
void createface(ParMeshAdapt* pMesh, int id1, int id2, double *face_centroid_1, double* face_controid_2, int tetra_ID, int &face_ID);
void print_CS_VTK(ParMeshAdapt* pMesh);

#endif
