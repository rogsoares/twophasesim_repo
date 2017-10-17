#ifndef SIMULATORPARAMETERS_H_
#define SIMULATORPARAMETERS_H_

#include "Boundary_Conditions.h"
#include "typedef.h"
#include "GeomData.h"


namespace PRS{


/**
 * SimulatorParameters class provides all IO data stream during simulation. Numeric, physical, pre-processor data input are made through this
 * class. When a specific class object is not applied to do so within SimulatorParameters class, it performs the data reading itself. VTK prin-
 * ting files are also made using SimulatorParameters.
 */
class SimulatorParameters{
public:

	SimulatorParameters();
	~SimulatorParameters();

	// read data from files.
	void deallocateData();
	void loadNumericParameters();
	void loadPhysicalParameters(int);
	void loadSolverParameters();
	void loadHighOrderParameter();

	// returns true if fimulation time was reached
	bool finishSimulation();

	// finishes simulation for some specific reason
	void stopSimulation() { stop_simulation = true; }
	int getEllipticSolver() const { return ellipticSolver; }
	int getHyperbolicSolver() const { return hyperbolicSolver; }
	string exportFilename() const { return expofName; }
	string getOutputPathName() const { return expofName; }
	string prepFilename() const { return prepfName; }

	// numeic parameters
	double CFL() const { return _CFL; }
	double getBC_Value(const int &flag);	// return value associated to flag
	bool isNodeFree(const int &flag);		// returns false if dirichlet
	double getWellTimeDiscretion() const { return _timeStepWell; } // Returns a number that is used to divide
	double getInitialSaturation(int flag);
	// the current time step (fractional time-step)


	/**
	 *  solver settings parameters
	 *  -------------------------------------------------------------------------
	 */
	double abstol() const { return _abstol; }
	double rtol() const { return _rtol; }
	double rtol2() const { return _rtol2; }
	double dtol() const { return _dtol; }
	int maxit() const { return _maxit;	}
	bool useDefectCorrection() const { return EBFV1_pressureSolver_scheme; }
	void setUseDefectCorrection() { EBFV1_pressureSolver_scheme = true; }
	PCType getPCType() const { return pctype; }
	double getPVIincrement() const { return PVI_increment; }
	int getKrylov_restart() const { return _Krylov_restart; }


	// physical parameters
	// -------------------------------------------------------------------------
	double oilDensity() const { return _oil_density; }
	double waterDensity() const { return _water_density; }
	double oilViscosity() const { return _oil_viscosity; }
	double waterViscosity() const { return _water_viscosity; }
	double getPorosity(const int &dom);
	const double* getPermeability(const int &dom);
	bool is_K_Isotropic() const { return K_Isotropic; }
	int ksModel() const { return _ksModel; }	// return relative permeability model

	// initial conditions
	// -------------------------------------------------------------------------
	double Sw_initial() const { return _Sw_initial; }
	double Sor() const { return _Sor; }
	double Swr() const { return _Swr; }
	double getInitialOilVolume() const { return _IOV; }


	/**
	 * Well management
	 * --------------------------------------------------------------------
	 */
	MapWells MWells;
	// for each well flag will be several nodes flagged
	// as edge loop is performed a node can be counted twice. So a set
	// container is used to store the ID from node
	typedef set<int> setNodesOnWell;
	map<int,setNodesOnWell> mapNodesOnWell;
	map<int,setNodesOnWell>::iterator MNOWIter;
//	bool isInjectionWell(pEntity node){ return isInjectionWell(GEN_tag(node->getClassification())); }
//	bool isProductionWell(pEntity node){ return isProductionWell(GEN_tag(node->getClassification())); }
	bool isInjectionWell(int flag) const;
	bool isProductionWell(int flag) const;
	double getTotalInjectionFlowRate() const;
	double getWellVolume(int well_flag) const;	// returns the sum of nodal volumes of the nodes well
	void getWells(ParMeshAdapt* pMesh);
	void WellFlowRate_weighted(ParMeshAdapt* pMesh, GeomData *);
	double getFlowrateValue(int flag) const;

	// reservoir geometric dimension used to make some physical properties dimensionless
	void getReservoirGeometricDimensions(double &L, double &H) const {
		H = reservoir_Height; L = reservoir_Length;
	}
	double dimensionlessFactorForDeltaTImplicit(double dom);
	double getDVTOL() const { return _dvtol; }
	void setInitialOilVolume(GeomData*);	// this work for multi-domains

	/*
	 * High order managemet
	 * --------------------------------------------------------------------
	 */
	bool useHOApproximation() const { return useHOApp; }
	NSLF getNodeSlopeLimitFunc() const { return slf_method_Node; }
	ESLF getEdgeSlopeLimitFunc() const { return slf_method_Edge; }
	double get_koef() const{ return koef; }
	double get_WoodfieldDelta() const { return _WFdelta; }

	// skip comments or any other unnecessary string until a specific point be reached
	void setPositionToRead(ifstream &fid, string str);

	/*
	 * Restart stuff
	 * --------------------------------------------------------------------
	 */
	bool useRestart() const { return restart; }
	string getRestartFilename() const { return restartFilename; }
	double getCumulativeSimulationTime() const { return cumTS; }
	void setCumulativeSimulationTime(double time_increment) { cumTS += time_increment; }

	int getStepOutputFile() const { return vtk_step; }
	void setStepOutputFile(int s) { vtk_step = s; }
	void incrementeStepOutputFile() { vtk_step++; }

	void setCPU_time(double cput) { cpu_time=cput; }
	double getCPU_time() const { return cpu_time; }
	void getSimParFiles();
	//int getStepCounter() const { return stepCounter; }
	int getLastPVI() const { return lastpvi; }
	void getExportFileName(string &expofn) const { expofn = expofName; }

	// Simulation must start from where it was before mesh adaptation.
	// Simulation time advance must be done over new adapted mesh
	// The time step calculated before mesh adaptation must be ignored
	void saveCurrentSimulationTimes(){
		currentST = cumTS;
	}

	void retrieveSimulationTime(){
		cumTS = currentST;
	}


	/*
	 * Simulation time monitoring
	 * --------------------------------------------------------------------
	 */
	double getPVIaccumulated() const { return PVI_cumulative; }
	void setPVIaccumulated(double pviacc) { PVI_cumulative = pviacc; }
	double getSimTime() const { return _ST; }
	double getPVI() const { return _PVI; }


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * If it's desired to print out a VTK file every fraction of PVI or every
	 * n days, timeStep must be corrected to reproduce the exact time evolution.
	 * For example, if every 0.1PVI a new VTK file must be generated it has:
	 *
	 * sum(dt) = n*PVI
	 * sum(dt) = dt_1 + dt_2 + ... + dt_n
	 *
	 * if ( sum(dt)>n*PVI ) dt_n = dt_n - (ST - n*PVI)
	 *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	void allowPrintVTK();
	void correctTimeStep(double&);
	void printOutVTK(void *pData1, void *pData2, void *pData3, pFunc_PrintVTK printVTK);
	double getPrintOutVTKFrequency();
	void setPrintOutVTKFrequency(double vtf) { vtk_time_frequency = vtf; }
	void updatePrintOutVTKFrequency();
	string getFilename(string f);
	void generateFilename(const string &meshFilename, char *filename);
	bool allowPrinting_VTK();
	void setNotAllowPrinting_VTK();

	void setTimetoPrintVTK(){
		allowPrintingVTK = true;
	}

	bool timeToPrintVTK()const{
		return allowPrintingVTK;
	}

//	bool SimulationHas_BC_ExternalDefinition() const{
//		return bc_external_definition;
//	}
	string command_line;
	bool use_first_order_approximation_on_contour();

	// says if boundary conditions values are defined by a function or not
	bool use_nonhomogneous_bc() const;

	// function pointer for exact solution
	void (*exact_solution)(const double *coords, double &sol);

	// function pointer for source/sink term
	void (*ss_term)(const double *coords, double &sol);

	void defineExactSolution();

	BENCHMARK getCaseProblem() const{
		return case_problem;
	}

	void setCommandLineArguments(char** argv, int argc){
		__argv = argv;
		__argc = argc;
	}

	// do not compute numeric solution. go to print exact solution in VTK
	bool run_Benchmark() const;

	bool adaptation_ocurred() const{
		return adapoccur;
	}

	void set_adapt_occur(bool k){
		adapoccur = k;
	}

	bool adapoccur;
	bool run_benchmark;
	bool bc_external_definition;
	char** __argv;
	int    __argc;

	BENCHMARK case_problem;


	// stores how many time steps are performed within every fraction of PVI
	std::list<int> TSCountingList;

	/*
	 * The following three variables below store user's decision to use, from
	 * a set of options, a specific implementation. The chosen option will be
	 * executed through a switch/case coding.
	 */
	NSLF slf_method_Node;		// related to high order method
	ESLF slf_method_Edge;		// related to high order method
	FRACTIONALFLUX fw_method; // related to fractional flux

	std::map<string,NSLF> mapSL_node;
	std::map<string,ESLF> mapSL_edge;
	std::map<string,double> map_koef;

	/*
	 *  koef: coefficient for Sw extrapolation used in high order implementation
	 *  Darlan's thesis (Carvalho, 2005), Eqs. (4.78) and (4.79), p92.
	 */
	double koef;

	// delta for woodfield high order method
	double _WFdelta;

	// store flags and their values into a specific map container
	void mapFlag(ifstream &fid, string whatmap);
	//void checkIfRankHasProductionWell();
	MapRockProperties mapRockProp;
	MapFlag mapBC;
	MapFlag mapSaturation;
	int numNodesOnInjectWell;

	// numeric parameters
	// -------------------------------------------------------------------------
	double _CFL;
	unsigned int _timeStepWell;
	unsigned int TimeStepCounter;

	PCType pctype;
	double _abstol;	// the absolute convergence tolerance
	double _rtol;	// the relative convergence tolerance
	double _rtol2;	// the relative convergence tolerance for external iteration
	double _dtol;	// the divergence tolerance
	int _maxit;		// maximum number of iterations
	bool EBFV1_pressureSolver_scheme;

	// physical parameters
	// -------------------------------------------------------------------------
	double _phi;
	double _oil_density;
	double _water_density;
	double _oil_viscosity;
	double _water_viscosity;
	double _Sw_initial;	// water saturation at t=0
	double _Sor;		// residual oil saturation
	double _Swr;		// irreducible water saturation
	int _ksModel;		// relative permeability model
	double _IOV;			// Initial Oil Volume

	bool _rankHasProductionWell;


	double _PVI;	// pore volume injected
	double _ST;		// simulation time
	bool well;		// inside isNodeFree() well is set to inform if node has a well or not
	int exportIter;

	void printParameters();

	/*
	 * Get from file and mesh wells (injection and production)
	 */
	void getWells(ifstream &fid);

	bool stop_simulation;
	int _Krylov_restart;
	double _dvtol;
	double reservoir_Height;
	double reservoir_Length;
	double initialOilVolume;
	int ellipticSolver;
	int hyperbolicSolver;

	// pointer for mesh data structure

	// strings for input/output files
	string expofName;
	string prepfName;

	/*
	 * If high order approximation is required by user, set useHOApp=true
	 */
	void setUseHOApproximation() { useHOApp = true; }
	bool useHOApp;

	/*
	 * Restart stuff
	 */
	bool restart;
	string restartFilename;
	double cumTS;			// cumulative time step: summation of all time steps
	double currentST;
	int vtk_step;
	double cpu_time;
	//int tsnumber;
	int stepCounter;

	/*
	 * string variable for file name
	 */
	string numeric_Filename;
	string physical_Filename;
	string slopeLimiter_Filename;
	string solver_Filename;

	/// Says when the next vtk file must be printed out
	double vtk_time_frequency;
	double PVI_increment;
	double PVI_cumulative;
	bool allowPrintingVTK;
	bool firstVTKupdate;
	int lastpvi;

	// map preconditioners
	void setPreconditionerDataBase();
	std::map< std::string, PCType> mapPC;

	// says where numeric.dat physical.dat slope_limiters.dat solve.dat are located
	void getParametersDirectoryPath();
	string parametersPath;

	void checkPermeabilityTensor();
	bool K_Isotropic;

	WellInfo well_info;

	// domains
	int *numNodesDom;
	int *numEdgesDom;
	map_set_int _mapNodesDomain;	      	// for each domain -> there is a set o of global node IDs
};
}
#endif /*SimulatorParameters_H_*/

