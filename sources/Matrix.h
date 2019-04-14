#pragma once
class Matrix {

private:
	double *values;
	int num_of_cols, num_of_rows;
public:
	Matrix(int rows, int cols);
	Matrix(const Matrix&);
	~Matrix();

	double &operator ()(int row, int col);
	double operator ()(int row, int col) const;
	void zeros();
	void diag_ones();
	void fill_diags(double a1, double a2, double a3);
	void fill_diag(double num, int translation);
	void draw() const;
	int get_cols();
	int get_cols() const;
	int get_rows();
	int get_rows() const;
	void copy_values(Matrix const & rhs);
	const Matrix operator*(Matrix const & a) const;
	Matrix& operator*=(Matrix const& rhs);
	const Matrix  operator-(Matrix const & a) const;
	void copy(Matrix const & copy);
};

