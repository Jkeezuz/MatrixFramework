#include "Matrix_Operations.h"
#include <math.h>
#include <iostream>
const Matrix  calc_res(Matrix const & x, Matrix const & a, Matrix const & b) {
	//where residuum is r in A*x' - b = r
	return (a * x) - b;
}
double norm(Matrix const & const m) {
	//Calculates the norm of M vector
	//using the square root of sums of squares
	int size;
	if (m.get_cols() > 1) {
		size = m.get_cols();
		double sum = 0;
		for (int i = 0; i < size; i++) {
			sum += m(0, i)*m(0, i);
		}
		sum = sqrt(sum);
		return sum;
	}
	else {
		size = m.get_rows();
		double sum = 0;
		for (int i = 0; i < size; i++) {
			sum += m(i, 0)*m(i, 0);
		}
		sum = sqrt(sum);
		return sum;
	}

}
const Matrix jacobiMethod(Matrix const & const a, Matrix const & const b) {
	//Jacobi method using the element-based formula
	//more info --> https://en.wikipedia.org/wiki/Jacobi_method#Description

	//We create all necessary variables and objects
	int bsize = b.get_rows();
	Matrix x(bsize, 1);
	Matrix x_new(bsize, 1);
	x_new.zeros();
	double sum = 0;
	double resnorm = 1;
	int iter_num = 0;

	//Core algorithm
	while (iter_num < 1000 && resnorm < 10e30 && resnorm > 10e-9) {
		//Copy the values of new x to the old 
		x.copy_values(x_new);
		
		//Perform algorithm based on the formula
		for (int i = 0; i < bsize; i++) {
			for (int j = 0; j < bsize; j++) {
				if (i != j) {
					sum += a(i, j)*x(j, 0);
				}
			}
			x_new(i, 0) = (b(i, 0) - sum) / (a(i, i));
			sum = 0;
		}
		//Calculate the residuum for new x vector
		Matrix res = calc_res(x_new, a, b);
		//Calculate the norm of residuum to have a way
		//to evaluate how close we're to the solution
		resnorm = norm(res);
		std::cout << "JAC NORM: " << resnorm << std::endl;
		iter_num += 1;

	}
	std::cout << "ITER: " << iter_num << std::endl;

	return x_new;
}
const Matrix gaussSeidelMethod(Matrix const & const a, Matrix const & const b) {
	//Gauss-Seidel method using the element-based formula
	//more info --> https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method#Description

	//We create all necessary variables and objects
	int bsize = b.get_rows();
	Matrix x(bsize, 1);
	Matrix x_new(bsize, 1);
	x_new.zeros();
	double sum = 0;
	double resnorm = 1;
	int iter_num = 0;

	while (iter_num < 1000 && resnorm < 10e30 && resnorm > 10e-9) {
		//Copy the values of new x to the old 
		x.copy_values(x_new);
		//Perform algorithm based on the formula
		for (int i = 0; i < bsize; i++) {
			for (int j = 0; j < i; j++) {
				sum += a(i, j)*x_new(j, 0);

			}
			for (int j = i + 1; j < bsize; j++) {
				sum += a(i, j)*x(j, 0);
			}
			x_new(i, 0) = (b(i, 0) - sum) / (a(i, i));
			sum = 0;
		}
		//Calculate the residuum for new x vector and
		//the norm of residuum to have a way
		//to evaluate how close we're to the solution
		resnorm = norm(calc_res(x_new, a, b));
		std::cout << "GS NORM: " << resnorm << std::endl;
		iter_num += 1;

	}
	std::cout << "ITER: " << iter_num << std::endl;
	return x_new;
}



const Matrix LUdecomp(Matrix const & const a, Matrix const & const b) {
	//LU decomposition algorithm using the Doolittle's Method
	//more info --> http://mathonline.wikidot.com/the-algorithm-for-doolittle-s-method-for-lu-decompositions
	//Performs decomposition of matrix "a" into two matrices
	// L and U
	//Later takes those matrices and uses backward and forward substitution
	//to get the solution vector
	int bsize = b.get_rows();
	Matrix u(a);
	//We create our U matrix as a matrix of zeroes with the size equal to a's size
	u.zeros();
	Matrix l(a.get_rows(), a.get_cols());
	//L should me a matrix of zeroes with diagonal filled with ones
	l.zeros();
	l.diag_ones();
	double sum = 0;
	//Algorithm based on the mathematical formulas
	//creates U and L matrices
	for (int k = 0; k < bsize; k++) {
		for (int m = k; m < bsize; m++) {
			for (int j = 0; j < k; j++) {
				sum += l(k, j)*u(j, m);
			}
			u(k, m) = a(k, m) - sum;
			sum = 0;
		}
		for (int i = k + 1; i < bsize; i++) {
			for (int j = 1; j < k; j++) {
				sum += l(i, j)*u(j, k);
			}
			l(i, k) = (a(i, k) - sum) / u(k, k);
			sum = 0;
		}
	}
	//We've got our L and U matrices
	//Perform forward substitution to solve Ly = B
	Matrix y = forward_subs(l, b);
	//Perform forward substitution to solve Ux = Y
	Matrix x = backwards_subs(u, y);
	//X is our solution vector
	//Calculate residuum and it's norm to evaluate
	//how close our x vector is to the real solution
	Matrix res = calc_res(x, a, b);
	double resnorm = norm(res);
	std::cout << "LUDECOMP NORM: " << resnorm << std::endl;
	return x;
}


const Matrix forward_subs(Matrix const & const a, Matrix const & const b) {
	//a = L
	//Performs forward substitution based on the mathematical description
	//more info --> https://algowiki-project.org/en/Forward_substitution#Mathematical_description_of_the_algorithm
	Matrix y(a.get_rows(), 1);
	double sum = 0;
	//Start of the algorithm
	y(0, 0) = b(0, 0);
	for (int m = 1; m < a.get_rows(); m++) {
		for (int i = 0; i < m; i++) {
			sum += (a)(m, i)*y(i, 0);
		}
		y(m, 0) = (b(m, 0) - sum) / (a)(m, m);
		sum = 0;

	}
	return y;
}
const Matrix backwards_subs(Matrix const & const a, Matrix const & const b) {
	//a = U
	//Performs backwards substitution based on the mathematical description
	//more info --> https://algowiki-project.org/en/Backward_substitution#Mathematical_description_of_the_algorithm
	Matrix x(a.get_rows(), 1);
	Matrix u(a);
	int size = a.get_rows();
	//Start of the algorithm
	x(size - 1, 0) = b(size - 1, 0) / u(size - 1, size - 1);
	double sum = 0;
	for (int m = size - 1; m >= 0; m--) {
		for (int i = m + 1; i < size; i++) {
			sum += (a)(m, i)*x(i, 0);
		}
		x(m, 0) = (b(m, 0) - sum) / (a)(m, m);
		sum = 0;

	}
	return x;
}