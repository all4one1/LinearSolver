#include "LinearSolver.h"

IterativeSolver::IterativeSolver()
{

}

void IterativeSolver::solveGS(double* f, double* f0, double* bb, int NN, SparseMatrix& M)
{
	for (k = 1; k < 100000; k++)
	{
		k++;
		double s = 0;
		for (int j = 0; j < NN; j++)
		{
			s = M.line(j, f);
			f[j] = f[j] + (bb[j] - s) / M[j][j];
		}

		double max = 0;
		double dif;
		for (int i = 0; i < NN; i++)
		{
			dif = abs(f0[i] - f[i]);
			if (dif > max)
				max = dif;
		}
		for (int j = 0; j < NN; j++)
			f0[j] = f[j];

		if (max < eps_iter)	break;
		if (k % 1000 == 0) cout << "host k = " << k << ", eps = " << max << endl;
	}
}

void IterativeSolver::solveJacobi(double* f, double* f0, double* bb, int NN, SparseMatrix& M)
{
	for (k = 1; k < 100000; k++)
	{
		double s = 0;
		for (int j = 0; j < NN; j++)
		{
			s = M.line(j, f0);
			f[j] = f0[j] + (bb[j] - s) / M[j][j];
		}

		double max = 0;
		double dif;
		for (int i = 0; i < NN; i++)
		{
			dif = abs(f0[i] - f[i]);
			if (dif > max)
				max = dif;
		}
		for (int j = 0; j < NN; j++)
			f0[j] = f[j];

		if (max < eps_iter)	break;
		if (k % 1000 == 0) cout << "host k = " << k << ", eps = " << max << endl;
	}

}

void IterativeSolver::write()
{
	if (write_i == 0) w.open("iter_solver.dat");
	w << " " << k << endl;
	write_i++;
}

void IterativeSolver::auto_test()
{
	using Matrix = std::vector<std::vector<double>>;
	Matrix A =
	{ { 30,3,4,0,0,0 },
	 { 4,22,1,3,0,0 },
	 { 5,7,33,6,7,0 },
	 { 0,1,2,42,3,3 },
	 { 0,0,2,11,52,2 },
	 { 0,0,0,3,9,26 } };

	SparseMatrix SM(6);
	SM = A;


	double b[6] = { 1, 2, 3, 3, 2, 1 };
	int n = 6;
	double* x = new double[6];
	double* x0 = new double[6];

	for (int i = 0; i < n; i++)
	{
		x[i] = x0[i] = 0;
	}

	solveGS(x, x0, b, n, SM);

	double cg[6] =	{ 0.1826929218e-1,	0.7636750835e-1,	0.5570467736e-1,	0.6371099009e-1,	0.2193724104e-1,	0.2351661001e-1 };
	
	cout << "x should be: ";	for (int i = 0; i < n; i++)		cout << cg[i] << " ";   cout << endl;
	cout << "auto test:   ";	for (int i = 0; i < n; i++)		cout << x[i] << " ";	cout << endl;

	cout << "b should be: 1 2 3 3 2 1" << endl;
	cout << "auto test:   ";
	for (int i = 0; i < n; i++)
	{
		double s = 0;
		for (int j = 0; j < n; j++)
		{
			s += A[i][j] * x[j];
		}
		cout << s << " ";
	} cout << endl;
}