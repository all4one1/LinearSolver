#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>

#include "FromOuterSparse/SparseMatrix.h"


using std::cout;
using std::endl;
using std::ofstream;

struct IterativeSolver
{
	int k = 0, write_i = 0, limit = 1000;
	double eps_iter = 1e-5;
	ofstream w;
	IterativeSolver();

	void solveGS(double* f, double* f0, double* bb, int NN, SparseMatrix& M);
	void solveJacobi(double* f, double* f0, double* bb, int NN, SparseMatrix& M);

	void write();
	void auto_test();
};
