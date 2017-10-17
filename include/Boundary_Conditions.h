/*
 * Boundary_Conditions.h
 *
 *  Created on: Jan 29, 2015
 *      Author: rogerio
 */

#ifndef BOUNDARY_CONDITIONS_H_
#define BOUNDARY_CONDITIONS_H_

#include "auxiliar.h"

/*
 * When prescribed value (Dirichlet) and prescribed flux (Neumann) were necessary, you MUST define them here.
 */

enum BENCHMARK {CASE_1, CASE_5, OBLIQUEDRAIN, NIKITIN};

/*
 * Artigo: 3D BENCHMARK ON DISCRETIZATION SCHEMES FOR ANISOTROPIC DIFFUSION PROBLEMS ON GENERAL GRIDS
 * Autores: Eymard, Henry, Herbin, Hubert, Klofkorn, Manzini
 * ---------------------------------------------------------------------------------------------------------
 */
/* Case 1: */
/* Exact Solution  : */ void Benchmark3D_case1__ES(const double *coords, double &sol);
/* Source/sink term: */ void Benchmark3D_case1__SST(const double *coords, double &sol);

/* Case 5: */
/* Exact Solution  : */ void Benchmark3D_case5__ES(const double *coords, double &sol);
/* Source/sink term: */ void Benchmark3D_case5__SST(const double *coords, double &sol);

/* Oblique-drain: */
/* Exact Solution  : */ void Benchmark3D_obliquedrain__ES(const double *coords, double &sol);


#endif /* BOUNDARY_CONDITIONS_H_ */
