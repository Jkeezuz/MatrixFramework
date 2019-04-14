#include "Matrix.h"

const Matrix  calc_res(Matrix const & x, Matrix const & a, Matrix const & b);
double norm(Matrix const & const m);
const Matrix jacobiMethod(Matrix const & const a, Matrix const & const b);
const Matrix gaussSeidelMethod(Matrix const & const a, Matrix const & const b);
const Matrix LUdecomp(Matrix const & const a, Matrix const & const b);
const Matrix forward_subs(Matrix const & const a, Matrix const & const b);
const Matrix backwards_subs(Matrix const & const a, Matrix const & const b);