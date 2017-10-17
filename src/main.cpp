#include "SIMULATION_core.h"

int main(int argc, char **argv){

	PRS::SIMULATION_core sim;
	try{
		sim.initialize(argc,argv);
		sim.start_simulation();
	}
	catch (Exception excp) {
		excp.showExceptionMessage();
	}

	// If an exception is thrown, finalize all necessary things before terminate running.
	sim.finalize();
	return 0;
}
