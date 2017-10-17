/*
 * Matrix.h
 *
 *  Created on: May 26, 2011
 *      Author: rogsoares
 */

#include "auxiliar.h"

#ifndef MATRIX_H_
#define MATRIX_H_

/*
 * This template class provides a simple way to create a matrix of any type and to access its values efficiently.
 * Values are stored in a array of pointers to arrays and not like objects, which will lead to a big overhead for
 * problems where data are continuously read like those found in numeric simulations. Mesh information like nodes
 * id and coordinates MUST NOT be stored like objects.
 */

template <class T>
class Matrix {

public:

	Matrix (){
		_rows = _cols = 0;
	}

	Matrix (int rows, int cols):_rows(rows),_cols(cols){
		mat = new T*[_rows];
		for (int row=0; row<_rows; row++){
			mat[row] = new T[_cols];
		}
	}

	Matrix (int rows, int cols,T val):_rows(rows),_cols(cols){
		mat = new T*[_rows];
		for (int row=0; row<_rows; row++){
			mat[row] = new T[_cols];
		}
		initialize(val);
	}

	~Matrix (){
		if (mat){
			for (int row=0; row<_rows; row++){
				delete mat[row];
			}
			delete [] mat; mat = 0;
		}
	}

	void allocateMemory(int rows, int cols){
		_rows = rows;
		_cols = cols;
		allocateMemory();
	}

	void allocateMemory(int rows){ // we are supposing matrix has one column: matrix is  a vector
		_rows = rows;
		_cols = 1;
		allocateMemory();
	}

	// allocate an array of pointers to reutilize memory already allocated.
	// it avoids allocate new memory space to store the same data.
	void allocate_array_pointers(int rows, int cols){
		_rows = rows;
		_cols = cols;
		try{
			mat = new T*[_rows];
			for (int row=0; row<_rows; row++){
				mat[row] = NULL;
			}
		}
		catch  (std::bad_alloc& ba){
			throw Exception(__LINE__,__FILE__,"Memory allocation failed. Bad allocation.");
		}
	}

	void allocateMemory(){
		try{
			mat = new T*[_rows];
			for (int row=0; row<_rows; row++){
				mat[row] = new T[_cols];
			}
		}
		catch  (std::bad_alloc& ba){
			throw Exception(__LINE__,__FILE__,"Memory allocation failed. Bad allocation.");
		}
	}

	void freeMemory(){
		for (int row=0; row<_rows; row++){
			delete[] mat[row];
		}
		delete[] mat;
		mat = 0;
	}

	void initialize(T val){
		if (!_rows || !_cols){
			throw Exception(__LINE__,__FILE__,"WARNING: Cannot initialize matrix with null size.\n");
		}
		for (int row=0; row<_rows; row++){
			for (int col=0; col<_cols; col++){
				mat[row][col] = val;
			}
		}
	}

	// set functions
	// ---------------------------------------------------------------------------------------------------

	void setValue(int row, int col, const T &val){
#ifdef MATRIX_DEBUG
		if (row<0 || row>=_rows || col<0 || col>=_cols){
			throw Exception(__LINE__,__FILE__,"Attempt of getting value out of bound\n");
		}
#endif //MATRIX_DEBUG
		mat[row][col] = val;
	}

	// reutilize memory already allocated
	void set_row(int row, T **vec){
#ifdef MATRIX_DEBUG
		if (row<0 || row>=_rows){
			throw Exception(__LINE__,__FILE__,"Attempt of getting value out of bound\n");
		}
#endif //MATRIX_DEBUG
		mat[row] = (*vec);
	}

//	void setRow(int row,const T *vec){
//		for (int row=0; row<_cols; row++){
//			mat[row][row] = vec[row];
//		}
//	}

	void setValue(int row, const T &val){
#ifdef MATRIX_DEBUG
		if (row<0 || row>=_rows){
			char msg[512]; sprintf(msg,"Attempt of setting value out of bound. row: %d. #rows: %d\n",row,_rows);
			throw Exception(__LINE__,__FILE__,msg);
		}
#endif //MATRIX_DEBUG
		mat[row][0] = val;
	}


	// get functions
	// ---------------------------------------------------------------------------------------------------
	void getsize(int &m, int&n){
		m = _rows;
		n = _cols;
	}

//	T& operator()( int row, int col ) {
//		if ((row<0 && row>=_rows) || (col<0 && col>=_cols))
//			throw Exception(__LINE__,__FILE__,"Attempt of getting value out of bound\n");
//		else
//			return mat[row][col];
//	}
//
//	T& operator()( int row) {
//		if ((row<0 && row>=_rows)){
//			throw Exception(__LINE__,__FILE__,"Attempt of getting value out of bound\n");
//		}
//		else{
//			return mat[row][0];
//		}
//	}

	void getRow(int row,const T **vec){
#ifdef MATRIX_DEBUG
		if (row<0 || row>=_rows){
			char msg[512]; sprintf(msg,"Attempt of getting value out of bound. row: %d. #rows: %d\n",row,_rows);
			throw Exception(__LINE__,__FILE__,msg);
		}
#endif //MATRIX_DEBUG
		*vec = mat[row]; 
	}

	const T* getrowconst(int row){
#ifdef MATRIX_DEBUG
		if (row<0 || row>=_rows){
			char msg[512]; sprintf(msg,"Attempt of getting value out of bound. row: %d. #rows: %d\n",row,_rows);
			throw Exception(__LINE__,__FILE__,msg);
		}
#endif //MATRIX_DEBUG
		return mat[row];
	}

	T* getrow(int row){
#ifdef MATRIX_DEBUG
		if (row<0 || row>=_rows){
			char msg[512]; sprintf(msg,"Attempt of getting value out of bound. row: %d. #rows: %d\n",row,_rows);
			throw Exception(__LINE__,__FILE__,msg);
		}
#endif //MATRIX_DEBUG
		return mat[row];
	}

	const T getValue(int row, int col) const{
#ifdef MATRIX_DEBUG
		if (row<0 || row>=_rows || col<0 || col>=_cols){
			char msg[512]; sprintf(msg,"Attempt of getting value out of bound. row: %d. #rows: %d. col: %d. #cols: %d\n",row,_rows,col,_cols);
			throw Exception(__LINE__,__FILE__,msg);
		}
#endif //MATRIX_DEBUG
		return mat[row][col];
	}

	T getValue(int row){
#ifdef MATRIX_DEBUG
		if (row<0 || row>=_rows){
			char msg[512]; sprintf(msg,"Attempt of getting value out of bound. row: %d. #rows: %d\n",row,_rows);
			throw Exception(__LINE__,__FILE__,msg);
		}
#endif //MATRIX_DEBUG
		return mat[row][0];
	}

	// print functions
	// ---------------------------------------------------------------------------------------------------

	void printNumRowsCols() const{
		cout << "rows: " << _rows <<"\tcols: " << _cols << endl;
	}

	void printfMatrix() const{
		for (int row=0; row<_rows; row++){
			for (int col=0; col<_cols; col++){
				cout << mat[row][col] << "  ";
			}
			cout << endl;
		}

	}

	void print(const char* filename) const{
		ofstream fid;
		fid.open(filename);
		fid << setprecision(6) << fixed;
		for (int row=0; row<_rows; row++){
			for (int col=0; col<_cols; col++){
				fid << mat[row][col] << "  ";
			}
			fid << endl;
		}
		fid.close();
	}


private:

	int _rows;
	int _cols;
	T **mat;

};

#endif /* MATRIX_H_ */
