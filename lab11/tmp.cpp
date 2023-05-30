#include "calerf.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;
const int D = 1;
const double t_min = 0.0; //czas poczatek
const double t_max = 2.0; //czas max
const double r_init = 1.0; // r podany parametr
const double a = 10.0; // podany parametr a
const double x_max = r_init + a; // gorna granica x
const double x_min = r_init; // dolna granica x
const double D_LAMBDA = 0.4;
const double I_LAMBDA = 1.0; //podane lambda
double h = 0.1;

// inicjuje wektor
double* initialize_vector(int n) {
	return new double[n];
}

// inicjuje macierz
double** initialize_matrix(int n, int m) {
	double** matrix = new double* [n];
	for (int i = 0; i < n; i++) {
		matrix[i] = new double[m];
	}
	return matrix;
}
void delete_vector(double* vector) {
	delete[] vector;
}
void delete_matrix(double** matrix, int n) {
	for (int i = 0; i < n; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
}
void vector_to_file(double* vector, int n, string filename) {
	fstream file(filename.c_str(), ios::out);
	if (file.is_open()) {
		for (int i = 0; i < n; i++) {
			file << vector[i] << endl;
		}
	}
}
void matrix_to_file(double** matrix, int n, int m, string filename) {
	fstream file(filename.c_str(), ios::out);
	if (file.is_open()) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				file << matrix[i][j] << ";";
			}
			file << "\n";
		}
	}
	file.close();
}

void plot_to_file(double** matrix, double* _x, int n, int m, double dt, string filename) {
	fstream file(filename.c_str(), ios::out);
	int t_05 = 0.5 / dt;
	int t_12 = 1.2 / dt;
	int t_18 = 1.8 / dt;

	// t = 0.5, t = 1.2, t = 1.8
	if (file.is_open()) {
		for (int x = 0; x < m; ++x) {
			file << _x[x] << " ";
			file << matrix[t_05][x];
			file << " " << matrix[t_12][x];
			file << " " << matrix[t_18][x] << endl;
		}
	}
	file.close();
}

// zwraca maksymalna wartosc z wektora
double vector_max(double* vector, int n) {
	double max = fabs(vector[0]);

	for (int i = 1; i < n; i++) {
		if (max < fabs(vector[i])) {
			max = vector[i];
		}
	}
	return max;
}

// rozwiazanie ukladu rownan dla thomasa
void thomas_algorithm(double* l, double* d, double* u, double* b, double* x, int m) {
	for (int i = 1; i < m-1; ++i) {
		d[i] = d[i] - (l[i - 1] / d[i - 1]) * u[i - 1];
		b[i] = b[i] - (l[i - 1] / d[i - 1]) * b[i - 1];
	}
	x[m - 1] = b[m - 1] / d[m - 1];
	for (int i = m - 2; i >= 0; --i) {
		x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
	}
}

// Warunek poczatkowy
void initial_condition(double** matrix, int n, int m) {
	//double x = x_min;
	for (int i = 0; i < m; i++) {
		matrix[0][i] = 1.0;
		//x += h;
	}
}

// warunki brzegowe
void boundary_conditions(double** matrix, int n, int m) {
	for (int i = 1; i < n; i++) {
		matrix[i][0] = 0.0;
		matrix[i][m - 1] = 1.0;
	}
}

// wartosc dokladna
double real_value(double x, double t) {
	return 1.0 - (r_init / x) * calerfpack::erfc_l( (x - r_init) / (2 * sqrt(D * t)));
}

// funkcja ustawiajaca siatke dokladnych rozwiazan
double** get_analytic_solution(int n, int m, double h, double dt) {
	double** analytic = initialize_matrix(n, m);
	double x = x_min;
	double t = t_min;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			analytic[i][j] = real_value(x, t);
			x += h;
		}
		x = x_min;
		t += dt;
	}
	analytic[0][0] = 1.0;
	return analytic;
}

//==============================================================
void CN_thomas_method(double** A, int n, int m) {
	double C_LAMBDA = 1.0 + I_LAMBDA;
	double* l = initialize_vector(m);
	double* d = initialize_vector(m);
	double* u = initialize_vector(m);
	double* b = initialize_vector(m);
	double* x = initialize_vector(m);


	for (int k = 1; k < n; k++) {

		l[0] = 0.0;
		d[0] = 1.0;
		u[0] = 0.0;
		//b[0] = A[k - 1][0];
		b[0] = 0.0;

		for (int i = 1; i < m - 1; i++) {
			l[i] = (I_LAMBDA / 2.0) * (1.0 - (h / (i * h)));
			d[i] = -C_LAMBDA;
			u[i] = (I_LAMBDA / 2.0) * (1.0 + (h / (i * h)));
			b[i] = -( (I_LAMBDA / 2.0 - (I_LAMBDA / 2.0) * (h / (i * h)) ) * A[k - 1][i - 1] + (1.0 - I_LAMBDA) * A[k - 1][i] + ((I_LAMBDA / 2.0) + (I_LAMBDA / 2.0) * (h / (i * h))) * A[k - 1][i + 1]);
		}

		l[m - 1] = 0.0;
		d[m - 1] = 1.0;
		u[m - 1] = 0.0;
		b[m - 1] = 1.0;

		thomas_algorithm(l, d, u, b, x, m);

		for (int i = 1; i < m - 1; i++) {
			A[k][i] = x[i];
		}

	}
	delete_vector(l);
	delete_vector(d);
	delete_vector(u);
	delete_vector(b);
	delete_vector(x);
}
//==============================================================
void laasonen_thomas_method(double** A, int n, int m) {
	double L_LAMBDA = 1.0 + 2.0 * I_LAMBDA;
	double* l = initialize_vector(m);
	double* d = initialize_vector(m);
	double* u = initialize_vector(m);
	double* b = initialize_vector(m);
	double* x = initialize_vector(m);


	for (int k = 1; k < n; k++) {

		l[0] = 0.0;
		d[0] = 1.0;
		u[0] = 0.0;
		//b[0] = A[k - 1][0];
		b[0] = 0.0;

		for (int i = 1; i < m - 1; i++) {
			l[i] = I_LAMBDA * (1.0 - (h / (i * h)));
			d[i] = -L_LAMBDA;
			u[i] = I_LAMBDA * (1.0 + (h / (i * h)));
			b[i] = -A[k - 1][i];
		}

		l[m - 1] = 0.0;
		d[m - 1] = 1.0;
		u[m - 1] = 0.0;
		b[m - 1] = 1.0;

		thomas_algorithm(l, d, u, b, x, m);

		for (int i = 1; i < m - 1; i++) {
			A[k][i] = x[i];
		}

	}
	delete_vector(l);
	delete_vector(d);
	delete_vector(u);
	delete_vector(b);
	delete_vector(x);
}
//========================================
double* residuum(double** a, double* b, double* x, int m) {
	double sum = 0.0;
	double* temp = new double[m];

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			sum += a[i][j] * x[j];
		}
		temp[i] = sum - b[i];
		sum = 0.0;
	}

	return temp;
}
//========================================
double est(double* x_new, double* x_old, int m) {
	double max = 0.0;
	double* p = new double[m];

	for (int i = 0; i < m; i++) {
		p[i] = x_new[i] - x_old[i];
	}

	if (fabs(p[0]) > fabs(p[1])) {
		max = fabs(p[0]);
	}
	else {
		max = fabs(p[1]);
	}

	for (int i = 0; i < m; i++) {
		if (fabs(p[i]) > max) {
			max = fabs(p[i]);
		}
	}

	delete[] p;
	return max;
}
void jacobi(double** A, double* b, double* x, int n, int m) {
	double TOL = 1e-8; 
	double* x_new = initialize_vector(m);
	double sum = 0.0;

	for (int i = 0; i < m; i++) {
		x_new[i] = x[i];
	}
	//x_new[0] = x[0];

	for (int iter = 0; iter < 70; iter++) {
		// (L+U)x
		for (int i = 0; i < m; i++) {
			sum = 0.0;
			for (int j = 0; j < m; j++) {
				if (j != i) {
					sum += (A[i][j] * x_new[j]);
				}
			}
			x_new[i] = (-sum + b[i]) / A[i][i];
		}

		// x_n => x_{n+1}
		for (int i = 0; i < m; i++) {
			x[i] = x_new[i];
		}

		if ((fabs((vector_max(residuum(A, b, x, m), m))) < TOL) && (fabs((est(x_new, x, m))) < TOL)) break;
	}
}
//==============================================================
void laasonen_jacobi_method(double** A, int n, int m) {
	double L_LAMBDA = 1.0 + 2.0 * I_LAMBDA;
	double* b = new double[m];
	double* wyn = new double[m];

	for (int i = 0; i < m; i++) {
		b[i] = 0.0; wyn[i] = 0.0;
	}

	double** new_A = initialize_matrix(m, m);

	// inicjalizacja macierzy zerami
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			new_A[i][j] = 0.0;
		}
	}

	for (int k = 1; k < n; k++) {

		new_A[0][0] = 1.0;
		//b[0] = A[k - 1][0];
		b[0] = 0.0; 

		// + | - -> 1
		for (int i = 1; i < m - 1; i++) {
			new_A[i][i + 1] = I_LAMBDA * (1.0 - (h / (i * h)));
			new_A[i][i] = -L_LAMBDA; //d
			new_A[i][i - 1] = I_LAMBDA * (1.0 + (h / (i * h)));
			b[i] = -A[k - 1][i];
		}

		b[m - 1] = 1.0;
		new_A[m - 1][m - 1] = 1.0;

		jacobi(new_A, b, wyn, n, m);

		// wypelnienie rzedu k-tego poza brzegami
		for (int i = 1; i < m - 1; i++) {
			A[k][i] = wyn[i];
		}
	}
}
//==============================================================
void CN_gauss_method(double** A, int n, int m) {
	double C_LAMBDA = 1.0 + I_LAMBDA;
	double* b = new double[m];
	double* wyn = new double[m];

	for (int i = 0; i < m; i++) {
		b[i] = 0.0; wyn[i] = 0.0;
	}

	double** new_A = initialize_matrix(m, m);

	// inicjalizacja macierzy zerami
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			new_A[i][j] = 0.0;
		}
	}

	for (int k = 1; k < n; k++) {

		new_A[0][0] = 1.0;
		//b[0] = A[k - 1][0];
		b[0] = 0.0;

		for (int i = 1; i < m - 1; i++) {
			new_A[i][i + 1] = (I_LAMBDA / 2.0) * (1.0 - (h / (i * h)));
			new_A[i][i] = -C_LAMBDA; //d
			new_A[i][i - 1] = (I_LAMBDA / 2.0) * (1.0 + (h / (i * h)));
			b[i] =  -( (I_LAMBDA/2.0 - (I_LAMBDA / 2.0) * (h / (i * h))) * A[k - 1][i-1] + (1.0 - I_LAMBDA) * A[k - 1][i]  + ( (I_LAMBDA / 2.0) + (I_LAMBDA / 2.0)* (h / (i * h))) * A[k - 1][i + 1]);
		}

		b[m - 1] = 1.0;
		new_A[m - 1][m - 1] = 1.0;

		jacobi(new_A, b, wyn, n, m);

		// wypelnienie rzedu k-tego poza brzegami
		for (int i = 1; i < m - 1; i++) {
			A[k][i] = wyn[i];
		}
	}
}
// siatka rozwiązania laasonen thomas
double** get_laasonen_thomas(int n, int m) {
	double** laasonen_thomas = initialize_matrix(n, m);

	initial_condition(laasonen_thomas, n, m);
	boundary_conditions(laasonen_thomas, n, m);

	laasonen_thomas_method(laasonen_thomas, n, m);

	return laasonen_thomas;
}

// siatka rozwiązania laasonen Gauss-Seidel
double** get_laasonen_gauss(int n, int m) {
	double** laasonen_gauss = initialize_matrix(n, m);

	initial_condition(laasonen_gauss, n, m);
	boundary_conditions(laasonen_gauss, n, m);

	laasonen_jacobi_method(laasonen_gauss, n, m);
	return laasonen_gauss;
}

// siatka rozwiązania laasonen thomas
double** get_CN_thomas(int n, int m) {
	double** CN_thomas = initialize_matrix(n, m);

	initial_condition(CN_thomas, n, m);
	boundary_conditions(CN_thomas, n, m);

	CN_thomas_method(CN_thomas, n, m);

	return CN_thomas;
}
// siatka rozwiązania laasonen thomas
double** get_CN_gauss(int n, int m) {
	double** CN_thomas = initialize_matrix(n, m);

	initial_condition(CN_thomas, n, m);
	boundary_conditions(CN_thomas, n, m);

	CN_gauss_method(CN_thomas, n, m);

	return CN_thomas;
}

// zwraca macierz bledow na siatce dyskretnej
double** get_error(double** realSolution, double** numericSolution, int n, int m) {
	double** errors = initialize_matrix(n, m);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			errors[i][j] = fabs(numericSolution[i][j] - realSolution[i][j]);
		}
	}

	return errors;
}
// zwraca wektor najwiekszych bledow w kazdym kroku czasowym
double* max_error(double** errors, int n, int m) {
	double* max_error = initialize_vector(n);

	for (int i = 0; i < n; i++) {
		max_error[i] = vector_max(errors[i], m);
	}

	return max_error;
}

// zwraca dt ze wzroru LAMBDA = D * dt/h^2
double get_dt(double lambda, double h, double D) {
	return (lambda * h * h) / D;
}

// wypelnienie wektora x
double* get_space_steps(double dt, int n, int m) {
	double* v = initialize_vector(m);
	double x = x_min;
	for (int i = 0; i < m; i++) {
		v[i] = x;
		x += h;
	}
	return v;
}

// wypelnienie wektora t
double* get_time_steps(double dt, int n, int m) {
	double* v = initialize_vector(n);
	double t = t_min;
	for (int i = 0; i < n; i++) {
		v[i] = t;
		t += dt;
	}
	return v;
}
double get_time_step(double dt, int n) {
	double t = t_min;
	for (int i = 0; i < n; i++) {
		t += dt;
	}
	return t;
}

int main() {
	cout << "GCC" << endl;

	double dt = get_dt(I_LAMBDA, h, D);

	// ilosc krokow t
	int n = ((t_max - t_min) / dt) + 1;
	// ilposc krokox x
	int m = ((x_max - x_min) / h) + 1;

	// zapisanie do pliku kroków 
	double* _x = get_space_steps(dt, n, m);
	double* _t = get_time_steps(dt, n, m);


	//-----------Rozwiazanie analityczne------------------------------------

	double** analytic = get_analytic_solution(n, m, h, dt);
	matrix_to_file(analytic, n, m, "LAB11_analytical.csv");
	plot_to_file(analytic, _x, n, m, dt, "LAB11_t_analytic.csv");

	//----------------------------------------------------------------------

#ifdef dobre
	//-----------laasonen_thomas-------------------------------------------------
	double** laasonen_thomas = get_laasonen_thomas(n, m);
	matrix_to_file(laasonen_thomas, n, m, "LAB11_laasonen_thomas_solution.csv");
	plot_to_file(laasonen_thomas, _x, n, m, dt, "LAB11_t_LT.csv");
	double** error_matrix_LT = get_error(analytic, laasonen_thomas, n, m);
	//matrix_to_file(errors, n, m, "LAB11_errors_laasonen_thomas.csv");
	double* errs_by_time_LT = max_error(error_matrix_LT, n, m);
	//vector_to_file(max_errs_by_time, n, "LAB11_max_errs_by_time_laasonen_thomas.csv");
	
	fstream file_LT("LAB11_t_err_LT.csv", ios::out);
	if (file_LT.is_open()) {
		for (int i = 0; i < n; i++) {
			file_LT << _t[i] << " "<< errs_by_time_LT[i] << endl;
		}
	}

	//fstream plots_LT("LAB11_plot_LT.csv", ios::out);
	//if (plots_LT.is_open()) {
	//	for (int i = 0; i < n; i++) {
	//		plots_LT << _t[i] << " " << errs_by_time_LT[i] << endl;
	//	}
	//}

	delete_matrix(laasonen_thomas, n);
	delete_matrix(error_matrix_LT, n);


	//-----------laasonen_gauss-------------------------------------------------
	double** laasonen_gauss = get_laasonen_gauss(n, m);
	matrix_to_file(laasonen_gauss, n, m, "LAB11_laasonen_gauss_solution.csv");
	plot_to_file(laasonen_gauss, _x, n, m, dt, "LAB11_t_LG.csv");
	double** error_matrix_LG = get_error(analytic, laasonen_gauss, n, m);
	double* errs_by_time_LG = max_error(error_matrix_LG, n, m);

	fstream file_LG("LAB11_t_err_LG.csv", ios::out);
	if (file_LG.is_open()) {
		for (int i = 0; i < n; i++) {
			file_LG << _t[i] << " " << errs_by_time_LG[i] << endl;
		}
	}
	delete_matrix(laasonen_gauss, n);
	delete_matrix(error_matrix_LG, n);


	//-----------CN_thomas-------------------------------------------------
	double** CN_thomas = get_CN_thomas(n, m);
	matrix_to_file(CN_thomas, n, m, "LAB11_CN_thomas_solution.csv");
	plot_to_file(CN_thomas, _x, n, m, dt, "LAB11_t_CNT.csv");
	double** error_matrix_CNT = get_error(analytic, CN_thomas, n, m);
	//matrix_to_file(errors, n, m, "LAB11_errors_laasonen_thomas.csv");
	double* errs_by_time_CNT = max_error(error_matrix_CNT, n, m);
	//vector_to_file(max_errs_by_time, n, "LAB11_max_errs_by_time_laasonen_thomas.csv");

	fstream file_CNT("LAB11_t_err_CNT.csv", ios::out);
	if (file_CNT.is_open()) {
		for (int i = 0; i < n; i++) {
			file_CNT << _t[i] << " " << errs_by_time_CNT[i] << endl;
		}
	}
	delete_matrix(CN_thomas, n);
	delete_matrix(error_matrix_CNT, n);


	//-----------CN_gauss-------------------------------------------------
	double** CN_gauss = get_CN_gauss(n, m);
	matrix_to_file(CN_gauss, n, m, "LAB11_CN_gauss_solution.csv");
	plot_to_file(CN_gauss, _x, n, m, dt, "LAB11_t_CNG.csv");
	double** error_matrix_CNG = get_error(analytic, CN_gauss, n, m);
	//matrix_to_file(errors, n, m, "LAB11_errors_laasonen_thomas.csv");
	double* errs_by_time_CNG = max_error(error_matrix_CNG, n, m);
	//vector_to_file(max_errs_by_time, n, "LAB11_max_errs_by_time_laasonen_thomas.csv");

	fstream file_CNG("LAB11_t_err_CNG.csv", ios::out);
	if (file_CNG.is_open()) {
		for (int i = 0; i < n; i++) {
			file_CNG << _t[i] << " " << errs_by_time_CNG[i] << endl;
		}
	}
	delete_matrix(CN_gauss, n);
	delete_matrix(error_matrix_CNG, n);

#endif

#ifndef bledy
	//------------zaleznosc maksymalnego | err | dla tmax od h

	fstream bledy, plik1, plik2;
	double krok = 0.1, error_LT, error_LG, error_CNT, error_CNG;

	bledy.open("LAB11_h_error.csv", fstream::out);
	for (krok; krok > 1e-2; krok /= 2)
	{
		h = krok;
		dt = get_dt(I_LAMBDA, h, D);
		// ilosc krokow t
		n = ((t_max - t_min) / dt) + 1;
		// ilposc krokox x
		m = ((x_max - x_min) / h) + 1;
		cout << krok << endl;

		double** analytic = get_analytic_solution(n, m, h, dt);

		//----------------------------------------------------------
		double** laasonen_thomas = get_laasonen_thomas(n, m);
		double** error_matrix_LT = get_error(analytic, laasonen_thomas, n, m);
		error_LT = vector_max(error_matrix_LT[n - 1], m);

		delete_matrix(laasonen_thomas, n);
		delete_matrix(error_matrix_LT, n);
		cout << "Finished LT" << endl;

		//----------------------------------------------------------
		double** laasonen_gauss = get_laasonen_gauss(n, m);
		double** error_matrix_LG = get_error(analytic, laasonen_gauss, n, m);
		error_LG = vector_max(error_matrix_LG[n - 1], m);

		delete_matrix(laasonen_gauss, n);
		delete_matrix(error_matrix_LG, n);
		cout << "Finished LG" << endl;

		//----------------------------------------------------------
		double** CN_thomas = get_CN_thomas(n, m);
		double** error_matrix_CNT = get_error(analytic, CN_thomas, n, m);
		error_CNT = vector_max(error_matrix_CNT[n - 1], m);

		delete_matrix(CN_thomas, n);
		delete_matrix(error_matrix_CNT, n);
		cout << "Finished CNT" << endl;

		//----------------------------------------------------------
		double** CN_gauss = get_CN_gauss(n, m);
		double** error_matrix_CNG = get_error(analytic, CN_gauss, n, m);
		error_CNG = vector_max(error_matrix_CNG[n - 1], m);

		delete_matrix(CN_gauss, n);
		delete_matrix(error_matrix_CNG, n);
		cout << "Finished CNG" << endl;

		delete_matrix(analytic, n);
		bledy << log10(krok) << " " << log10(error_LT) << " " << log10(error_LG) << " " << log10(error_CNT) << " " << log10(error_CNG)<< endl;
	}
	bledy.close();

#endif
	//delete_matrix(analytic, n);

	//delete_matrix(errors, n);
	//delete_vector(max_err);
	//delete_vector(_x);
	//delete_vector(_t);
	// 
	// 	///*vector_to_file(_x, m, "LAB11_x.csv");
	//*/vector_to_file(_t, n, "LAB11_t.csv");

	////-------------------------------METODA LAASONEN + jacobi-------------------------------
	//double dt = get_dt(I_LAMBDA, h, D);

	//int n = ((t_max - t_min) / dt);
	//int m = ((x_max - x_min) / h);

	//double** analytic = get_analytic_solution(n, m, h, dt);
	//matrix_to_file(analytic, n, m, "analytical_1.csv");
	// 
	// 
	//double** errors = get_error(analytic, laasonen, n, m);
	//double* max_err = max_error(errors, n, m);
	//vector_to_file(max_err, n, "max_err_laasonen_jacobi.csv");
	//matrix_to_file(errors, n, m, "errors_laasonen_jacobi.csv");
	//double* _x = get_space_steps(dt, n, m);
	//double* _t = get_time_steps(dt, n, m);
	//vector_to_file(_x, n, "space_steps_laasonen_jacobi.csv");
	//vector_to_file(_t, n, "time_steps_laasonen_jacobi.csv");
	//delete_matrix(analytic, n);
	//delete_matrix(laasonen, n);
	//delete_matrix(errors, n);
	//delete_vector(max_err);
	//delete_vector(_x);
	//delete_vector(_t);
	return 0;
}