/*
 * Boundary_Conditions.cpp
 *
 *  Created on: Jan 29, 2015
 *      Author: rogerio
 */

#include "Boundary_Conditions.h"

void Benchmark3D_case1__ES(const double *coords, double &sol){
	double x = coords[0];
	double y = coords[1];
	double z = coords[2];
	sol = 1.0 + sin(pi*x) * sin(pi*(y+.5)) * sin(pi*(z+.33333333333));
}

void Benchmark3D_case1__SST(const double *coords, double &sol){
	double x = coords[0];
	double y = coords[1];
	double z = coords[2];
	double A = double(sin(pi*x));
	double B = double(sin(pi*(y+.5)));
	double C = double(sin(pi*(z+.3333333333)));
	double D = double(cos(pi*x));
	double E = double(cos(pi*(y+.5)));
	double F = double(cos(pi*(z+.3333333333)));
	sol = pi*pi*( 3.0*A*B*C - D*E*C - A*E*F );
}

void Benchmark3D_case5__ES(const double *coords, double &sol){
	double alpha = 0;
	double x = coords[0];
	double y = coords[1];
	double z = coords[2];
	if (y<=.5 && z<=.5){
		alpha = 0.1;
	}
	else if (y>.5 && z<=.5){
		alpha = 10;
	}
	else if (y>.5 && z>.5){
		alpha = 100;
	}
	else if (y<=.5 && z>.5){
		alpha = 0.01;
	}

	double excsol = (double)(alpha*sin(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*z));
//	if (fabs(excsol)<1e-15){
//		excsol = .0;
//	}
	sol = excsol;
}

void Benchmark3D_case5__SST(const double* coords, double &sol){
	double alpha=0, ax=0, ay=0, az=0;
	double x = coords[0];
	double y = coords[1];
	double z = coords[2];
	if (y<=.5 && z<=.5){
		alpha = 0.1;
		ax = 1.;
		ay = 10.;
		az = 0.01;
	}
	else if (y>.5 && z<=.5){
		alpha = 10.;
		ax = 1.;
		ay = 0.1;
		az = 100.;
	}
	else if (y>.5 && z>.5){
		alpha = 100.;
		ax = 1.;
		ay = 0.01;
		az = 10.;
	}
	else if (y<=.5 && z>.5){
		alpha = 0.01;
		ax = 1.;
		ay = 100.;
		az = 0.1;
	}

	double c = (4.*pi*pi*alpha*sin(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*z));
//	if (fabs(c)<1e-15){
//		c = 0;
//	}
	sol = (ax + ay+ az)*c;
}

void Benchmark3D_obliquedrain__ES(const double *coords, double &sol){
	double x = coords[0];
	double y = coords[1];
	sol = (double)(-x - 0.2*y);
}
