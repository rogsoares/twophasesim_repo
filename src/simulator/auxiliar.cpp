#include "auxiliar.h"


void debug_msg(string str){
//#ifdef DEBUG
	cout << str << endl;
//#endif
}

void vector_norm(const double* v, int size, double &norm){
	norm = 0;
	for(int i=0; i<size; i++){
		norm += v[i]*v[i];
	}
	norm = sqrt(norm);
}

void create_vector(const double* a, const double* b, double* v){
	v[0] = b[0]-a[0];
	v[1] = b[1]-a[1];
	v[2] = b[2]-a[2];
}

void dot_product(const double* a, const double* b, double &dot){
	dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void cross_product(const double* a, const double* b, double* cross){
	cross[0] = a[1]*b[2] - b[1]*a[2];
	cross[1] = a[2]*b[0] - a[0]*b[2];
	cross[2] = a[0]*b[1] - a[1]*b[0];
}

void throw_exception(bool condition, string msg, int line, const char* sourcefile){
	if (condition){
		throw Exception(line,sourcefile,msg.c_str());
	}
}

void open_file(ofstream& fid, string filename, int line, const char* sourcefile){

	fid.open(filename.c_str());
	char msg[512]; sprintf(msg,"File '%s' could not be opened or it does not exist.\n"
			"\tCheck if directory was typed correctly.\n",filename.c_str());

	if ( !fid.is_open() ){
		throw Exception(line,sourcefile,msg);
	}
}

void open_file(ifstream& fid, string filename, int line, const char* sourcefile){

	fid.open(filename.c_str());
	char msg[512]; sprintf(msg,"File '%s' could not be opened or it does not exist.\n"
			"\tCheck if directory was typed correctly.\n",filename.c_str());

	if ( !fid.is_open() ){
		throw Exception(line,sourcefile,msg);
	}
}

void alloc_BOOL_vector(int LINE, const char* FILE, bool* &p, int size){
	try{
		p = new bool[size];
	}
	catch  (std::bad_alloc& ba){
		throw Exception(LINE,FILE,"Memory allocation failed. Bad allocation.");
	}
}

void dealloc_BOOL_vector(bool* &p){
	delete[] p; p = 0;
}

void alloc_INT_vector_init(int LINE, const char* FILE, int** p, int size){
	try{
		*p = new int[size];
		for (int i=0; i<size; i++){
			(*p)[i] = i;
		}
	}
	catch  (std::bad_alloc& ba){
		throw Exception(LINE,FILE,"Memory allocation failed. Bad allocation.");
	}
}

void alloc_INT_vector(int LINE, const char* FILE, int* &p, int size){
	try{
		p = new int[size];
		for (int i=0; i<size; i++){
			p[i] = 0;
		}
	}
	catch  (std::bad_alloc& ba){
		throw Exception(LINE,FILE,"Memory allocation failed. Bad allocation.");
	}
}

void dealloc_INT_vector(int* &p){
	delete[] p; p = 0;
}

void alloc_DOUBLE_vector(int LINE, const char* FILE, double* &p, int size){
	try{
		p = new double[size];
		for (int i=0; i<size; i++){
			p[i] = 0;
		}
	}
	catch  (std::bad_alloc& ba){
		throw Exception(LINE,FILE,"Memory allocation failed. Bad allocation.");
	}
}

void dealloc_DOUBLE_vector(double* &p){
	delete[] p; p = 0;
}


void failOpeningFile(string filename, int line, const char* cppfile){
	char msg[256]; sprintf(msg,"File '%s' could not be opened or do not exist.\n ",filename.c_str());
	throw Exception(line,cppfile,msg);
}

// converts a string read from file to a double
double strToDouble(string &str){
	return strtod(getSubString(str), 0);
}

// converts a string read from file to an interger
int strToInteger(string &str){
	return atoi(getSubString(str));
}

// filters string the numeric part - internal use
const char* getSubString(string &str){
	string::size_type loc = str.find( "=", 0 );
	string numberstr = str.substr(loc+1, str.size()-loc);
	return numberstr.c_str();
}

void seek_position_on_file(ifstream &fid, string str2){
	int count = 0;
	string str1;
	do{
		getline(fid,str1);
		if (++count > 100){

			cout << "\n\nWARNING: Could not find <" << str2
					<<">.\nProbably you are using an out of date input file."
					" Default value will be used instead.\n\n";
			break;
		}
	}while( str1.compare(str2) );
}

void replaceAllOccurencesOnString(string &theString, string::size_type size,string seekFor, string replaceFor){
	string::size_type pos = 0;
	while (true){
		pos = theString.find(seekFor,pos);
		if (pos != string::npos)
			theString.replace(pos,size,replaceFor);
		else
			break;
	}
}

void printSimulationHeader(){
	std::cout<< "\n\n\t\t\t=====================================================================================\n"
			"\t\t\t\t\tUNIVERSIDADE FEDERAL DE PERNAMBUCO\n"
			"\t\t\t\t\tCENTRO ACADEMICO DO AGRESTE - CAA\n"
			"\t\t\t\t\tNUCLEO DE TECNOLOGIA - NT\n"
			"\t\t\t\t\tPADMEC - DEMEC - CNPq\n"
			"\t\t\t\t\t2008-2014\n"
			"\n\t\t\t\t\tPRS: Parallel Reservoir Simulator\n"
			"\t\t\t\t\tA Finite Volume Scheme for the Two-Phase Oil-Water Flow in 2-D and\n\t\t\t\t\t3-D heterogeneous/anisotropic porous media\n"
			"\n\t\t\t\t\tAuthor: Rogerio Soares da Silva.\n"
			"\t\t\t\t\tNO WARRANTY. IT MEANS: USE IT AT YOUR OWN RISK\n"
			"\t\t\t=====================================================================================\n\n";
}

void convertSecToTime(double t, double *h, double *m, double *s){
	double frac;
	frac = modf(t/3600.,h);
	frac = modf(frac*60.,m);
	frac = modf(frac*60.,s);
}

void convertSecToTime(double t){
	double frac,h,m,s;
	frac = modf(t/3600.,&h);
	frac = modf(frac*60.,&m);
	frac = modf(frac*60.,&s);
	cout << setprecision(0) << "Time: " << h << "h " << m << "min " << s << "s\n";
}


void LogFiles(double timeStep, double assemblyT, double solverT, double gradT, int KSPiter, double hyperbolicCPU,LOG_FILES LG, string path,
		bool restart, int last_step, double cumulativeSTime_Restart, double CPUTime_Restart){

	static int step_counter = (restart)?last_step:0;
	static double cumulativeSTime = (restart)?cumulativeSTime_Restart:0;
	static double cumulativeCPU = (restart)?CPUTime_Restart:0; // cumulative CPU time.
	static ofstream fid;

	switch (LG){
	case OPENLG:{
		// print results on file
		char filename[256]; sprintf(filename,"%s_simulation-monitor-%d.csv",path.c_str(),1);
		if (restart){
			fid.open(filename,ios_base::app);
		}
		else{
			fid.open(filename);
			fid << "#step time-step PVI cumulative-TS assembly-CPU solver-CPU grad-CPU KSP-iter hyperbolic-CPU cumulative-CPU\n";
		}
		std::cout << "\n--------------------------------------------------\nStart Simulation\n--------------------------------------------------\n\n";
	}
	break;
	case UPDATELG:{
		// rank 0 is in charge to print the output CPU-time

		cumulativeSTime += timeStep;
		cumulativeCPU += assemblyT + solverT + gradT + hyperbolicCPU;
		fid << scientific << setprecision(8);
		fid << ++step_counter << " " << timeStep << " " << (int)(100*cumulativeSTime/0.07) << " " << cumulativeSTime << " " << assemblyT << " " <<  solverT
				<< " " << gradT << " " << KSPiter << " " << hyperbolicCPU << " " << cumulativeCPU << endl;

	}
	break;
	case CLOSELG:{

		fid.close();

	}
	}
}

void STOP(){
	MPI_Barrier(MPI_COMM_WORLD); throw 1;
}

int printMatrixToFile(Mat& m,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	MatView(m,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

int printVectorToFile(Vec& v,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	VecView(v,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

