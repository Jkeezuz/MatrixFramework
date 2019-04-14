#include "Matrix.h"
#include "Matrix_Operations.h"
#include <iostream>
Matrix::Matrix(int rows, int cols) : num_of_rows(rows), num_of_cols(cols) {
	//Initializes new array for values with the size equal to number of rows times num of colums
	if (rows == 0 || cols == 0)
		throw ("Matrix constructor has 0 size");
	values = new double[rows * cols];
	//values = (double*)malloc(sizeof(double)*rows*cols);
}
Matrix::Matrix(const Matrix& copy) {
	//Initialize all values
	this->num_of_cols = copy.num_of_cols;
	this->num_of_rows = copy.num_of_rows;
	this->values = new double[num_of_cols*num_of_rows];
	//Copy the values of original matrix to our copy version
	this->copy_values(copy);
}
Matrix::~Matrix()
{
	//Free memory
	delete[] values;
}
int Matrix::get_cols() {
	//Returns number of columns in matrix
	return num_of_cols;
}
int Matrix::get_cols() const{
	//Returns number of columns in matrix
	return num_of_cols;
}
int Matrix::get_rows() {
	//Returns number of rows in matrix
	return num_of_rows;
}
int Matrix::get_rows() const{
	//Returns number of rows in matrix
	return num_of_rows;
}
void Matrix::copy(Matrix const & rhs) {
	//"Copies" the rhs matrix by changing the size and copying values
	this->num_of_cols = rhs.num_of_cols;
	this->num_of_rows = rhs.num_of_rows;
	delete this->values;
	this->values = new double[rhs.num_of_cols * rhs.num_of_rows];
	this->copy_values(rhs);
}
void Matrix::copy_values(Matrix const & rhs) {
	//Copies values of rhs matrix
	if (this->num_of_cols != rhs.num_of_cols || this->num_of_rows != rhs.num_of_rows) {
		//Throw exception if the sizes don't match
		throw "Bad dimensions";
	}
	for (int i = 0; i < this->num_of_cols*this->num_of_rows; i++) {
		this->values[i] = rhs.values[i];
	}
}
Matrix& Matrix::operator*=(Matrix const& rhs) {
	//Overloads *= operator for easier multiplying of matrices
	int rows = rhs.get_rows();
	int cols = rhs.get_cols();
	//Basic size check before we can multiply
	if (this->get_cols() != rows) {
		throw "Bad dimensions";
	}
	double sum = 0;
	Matrix c(this->get_rows(), cols);
	//Classic multiply algorithm
	for (int i = 0; i < c.get_rows(); i++) {
		for (int k = 0; k < c.get_cols(); k++) {
			for (int j = 0; j < this->get_cols(); j++) {
				sum += (*this)(i, j)*rhs(j, k);
			}
			c(i, k) = sum;
			sum = 0;
		}
	}
	//
	this->copy(c);
	//We want to return the same object to allow chains of operators
	return *this;
}
const Matrix Matrix::operator*(Matrix const& rhs) const {
	//Synctatic sugar, uses our *= operator but allows to use * in the code for easier coding
	Matrix result(*this);
	return result *= rhs;
}
const Matrix Matrix::operator-(Matrix const& b) const {
	int rows = b.get_rows();
	int cols = b.get_cols();
	if (this->get_rows() != rows || this->get_cols() != cols) {
		//Throws exception if sizes dont match
		throw "Bad dimensions";
	}
	double sum = 0;
	//Classic sub algorithm
	for (int i = 0; i < rows; i++) {
		for (int k = 0; k < cols; k++) {
			this->values[num_of_cols*i + k] -= b(i, k);
		}
	}
	return *this;
}
double &Matrix::operator()(int row, int col) {
	//Returns the value of matrix at given row and column
	if (row >= num_of_rows || col >= num_of_cols)
		throw ("Matrix subscript out of bounds");
	//Eg. row = 1, col = 1, number of columns is 2 then
	//return values[3]
	return values[num_of_cols*row + col];
}
double Matrix::operator()(int row, int col) const {
	//Returns the value of matrix at given row and column
	if (row >= num_of_rows || col >= num_of_cols)
		throw ("Matrix subscript out of bounds");
	//Eg. row = 1, col = 1, number of columns is 2 then
	//return values[3]
	return values[num_of_cols*row + col];
}
void Matrix::fill_diags(double a1, double a2, double a3) {
	//fills 5 diagonals 
	fill_diag(a1, 0);
	fill_diag(a2, 1);
	fill_diag(a2, -1);
	fill_diag(a3, 2);
	fill_diag(a3, -2);
}
void Matrix::fill_diag(double num, int translation) {
	//fills diagonal with num number
	//translation is the translation of our diagonal to fill
	//in respect to the main diag (0)
	//so the diagonal just above the main is 1 and just below is -1 etc.
	if (translation > 0) {
		for (int i = translation; i < this->num_of_rows; i++) {
			this->values[num_of_cols*i + i - translation] = num;
		}
	}
	else {
		translation *= -1;
		for (int i = 0; i < this->num_of_rows - translation; i++) {
			this->values[num_of_cols*i + i + translation] = num;
		}
	}
}
void Matrix::draw() const {
	//Draws the matrix on standard output for debug purposes
	for (int i = 0; i < num_of_rows; i++) {
		for (int j = 0; j < num_of_cols; j++) {
			std::cout << this->values[num_of_cols*i + j] << "|";
		}
		std::cout << std::endl;
	}
}
void Matrix::zeros() {
	//Fills the whole matrix with zeroes
	for (int i = 0; i < num_of_rows; i++) {
		for (int j = 0; j < num_of_cols; j++) {
			this->values[num_of_cols*i + j] = 0;
		}
	}
}
void Matrix::diag_ones() {
	//Fills the diagonal of matrix with ones
	for (int i = 0; i < num_of_rows; i++) {
		this->values[num_of_cols*i + i] = 1;
	}
}









