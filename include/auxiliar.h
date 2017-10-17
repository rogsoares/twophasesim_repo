#ifndef AUXILIAR_H_
#define AUXILIAR_H_

#include "includes.h"

typedef std::map<int,double> map_int_double;
typedef std::map<int,double>::iterator map_int_double_Iter;

typedef std::map<int,double*> map_int_pDouble;
typedef std::map<int,double*>::iterator map_int_pDouble_Iter;

typedef std::map<int,int> map_int_int;
typedef std::map<int,int>::iterator map_int_int_Iter;
typedef std::set<int> set_int;
typedef std::set<int>::iterator set_int_Iter;

typedef std::map<int,set<int> > map_set_int;
typedef std::map<int,set<int> >::iterator map_set_int_Iter;

void open_file(ofstream& fid, string filename, int line, const char* sourcefile);
void open_file(ifstream& fid, string filename, int line, const char* sourcefile);
void throw_exception(bool,string, int line, const char* sourcefile);

void debug_msg(string str);

void alloc_BOOL_vector(int LINE, const char* FILE, bool* &p, int size);
void dealloc_BOOL_vector(bool* &p);
void alloc_INT_vector(int LINE, const char* FILE, int* &p, int size);
void dealloc_INT_vector(int* &p);
void alloc_DOUBLE_vector(int LINE, const char* FILE, double* &p, int size);
void dealloc_DOUBLE_vector(double* &p);

void alloc_INT_vector_init(int LINE, const char* FILE, int** p, int size);


enum FIELD {PRESSURE, SATURATION};
const double pi = 3.14159265358979323846;
double strToDouble(string &str);
void convertSecToTime(double t);
int strToInteger(string &str);
const char* getSubString(string &str);
void replaceAllOccurencesOnString(string &, string::size_type, string, string);
void seek_position_on_file(ifstream &fid, string str2);


//void F_getEdges(pMesh, pEntity, std::vector<pEntity> &);

const double qsi = 1e-10; // qsi is used to avoid division by zero
void printSimulationHeader();
void checklinepassing(int, const char*, std::string);
enum LOG_FILES {OPENLG, UPDATELG, CLOSELG};
void LogFiles(double timeStep, double assemblyT, double solverT, double gradT, int KSPiter,
		      double hyperbolicCPU,LOG_FILES LG, string path, bool restart, int last_step,
			  double cumulativeSTime_Restart, double CPUTime_Restart);
void failOpeningFile(string, int, const char *);
void convertSecToTime(double t, double *h, double *m, double *s);
void STOP();
int printMatrixToFile(Mat&,const char*);
int printVectorToFile(Vec&,const char*);

void cross_product(const double* a, const double* b, double* cross);
void dot_product(const double* a, const double* b, double &dot);
void create_vector(const double* a, const double* b, double* v);
void vector_norm(const double* v, int size, double &norm);

///*
// * Define types for arrays of pointer functions (Scalars)
// */
//typedef double(*GetPFunction)(pEntity);
//typedef void(*SetPFunction)(pEntity,double);
//
//typedef GetPFunction* GetFunctionArray;
//typedef SetPFunction* SetFunctionArray;
//
///*
// * Define types for arrays of pointer functions (scalars/gradients)
// */
//typedef double (*GetPFuncScalar)(pEntity);
//typedef void (*SetPFuncScalar)(pEntity,double);
//typedef void (*GetPFuncGrad)(int, int, double*);
//typedef void (*SetPFuncGrad)(int, int, double*);
//
//typedef GetPFuncScalar* GetFuncScalarArray;
//typedef SetPFuncScalar* SetFuncScalarArray;
//typedef GetPFuncGrad* GetFuncGradArray;
//typedef SetPFuncGrad* SetFuncGradArray;

//typedef void(*FuncPointer_GetGradient)(FIELD,int,int,int,double*);


#endif /*AUXILIAR_H_*/
