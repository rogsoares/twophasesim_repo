/*
 * SimulatorParameters2.cpp
 *
 *  Created on: Jan 29, 2015
 *      Author: rogerio
 */

#include "SimulatorParameters.h"

namespace PRS{

	bool SimulatorParameters::run_Benchmark() const{
			return run_benchmark;
		}
	bool SimulatorParameters::use_first_order_approximation_on_contour(){
		return ( command_line.find("use_first_order_approximation_on_contour") != std::string::npos )?true:false;
	}

	// says if boundary conditions values are defined by a function or not
	bool SimulatorParameters::use_nonhomogneous_bc() const{
		return bc_external_definition;
	}

	void SimulatorParameters::defineExactSolution(){
		bc_external_definition = true;

		// put all command line strings into one string object
		for(int i=0; i<__argc; i++){
			command_line.append(__argv[i]);
		}

		/*
		 * Run benchmark case using the exact solution as output.
		 */
		run_benchmark = ( command_line.find("run_benchmark") != std::string::npos )?true:false;

		//cout << "Print exact solution is " << exact_sol_exist << endl;

		if ( command_line.find("benchmark3d_case1") != std::string::npos ){
			exact_solution = Benchmark3D_case1__ES;
			ss_term = Benchmark3D_case1__SST;
			case_problem = CASE_1;
			cout << "Benchmark: case 1\n";
		}
		else if ( command_line.find("benchmark3d_case5") != std::string::npos ){
			exact_solution = Benchmark3D_case5__ES;
			ss_term = Benchmark3D_case5__SST;
			case_problem = CASE_5;
			cout << "Benchmark: case 5\n";
		}
		else if ( command_line.find("benchmark3d_obliquedrain") != std::string::npos ){
			exact_solution = Benchmark3D_obliquedrain__ES;
			ss_term = 0;
			case_problem = OBLIQUEDRAIN;
			cout << "Benchmark: Oblique drain\n";
		}
		else{
			bc_external_definition = false;
		}

		if (use_first_order_approximation_on_contour()){
			cout << "Using first order approximation on contour\n";
		}
	}
}
