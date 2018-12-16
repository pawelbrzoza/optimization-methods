//Do not edit the code below (unless you know what you are doing)

#include"solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = matrix(L);
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(double *A, int n)
{
	x = matrix(A, n);
	g = NAN;
	H = NAN;
	y = NAN;
}

//You can edit the following code

void solution::fit_fun()
{
	//*********************** LAB 2 ***********************************

	/*matrix Y0(new double[3]{ 5,1,10 }, 3);
	matrix *Y = solve_ode(0, 1, 1000, Y0, x);
	int *w = get_size(Y[1]);
	double max = Y[1](0, 2);
	for (int i = 1; i < w[0]; ++i)
		if (max < Y[1](i, 2)) max = Y[1](i, 2);
	y = abs(max - 50);
	++f_calls;*/

	//*********************** LAB 3 ***********************************

	double a_ref = 3.14, o_ref = 0;
	matrix Y0(2, 1);
	matrix *Y = solve_ode(0, 0.1, 100, Y0, x);
	int *n = get_size(Y[1]);
	y(0) = 0;
	for (int i = 0; i < n[0]; ++i)
		y = y + 10 * pow(a_ref - Y[1](i, 0), 2) +
		pow(o_ref - Y[1](i, 1), 2) +
		pow(x(0)*(a_ref - Y[1](i, 0)) + x(1) * (o_ref - Y[1](i, 1)), 2);
	y = y * 0.1;
	++f_calls;
}

void solution::fit_fun_outside(matrix A)
{
	//*********************** LAB 4 ***********************************
	double arg = 3.14*sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
	y = sin(arg) / arg;

	if (-x(0) + 1 > 0)
		y = y + A(0)*pow(-x(0) + 1, 2);
	if (-x(1) + 1 > 0)
		y = y + A(0)*pow(-x(1) + 1, 2);
	if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1) > 0)
		y = y + A(0)*pow(sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1), 2);

	++f_calls;
}

void solution::fit_fun_inside(matrix A)
{
	//*********************** LAB 4 ***********************************
	double arg = 3.14*sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
	y = sin(arg) / arg;

	if (-x(0) + 1 > 0)
		y = 1e10;
	else
		y = y - A(0) / (-x(0) + 1);
	if (-x(1) + 1 > 0)
		y = 1e10;
	else
		y = y - A(0) / (-x(1) + 1);
	if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1) > 0)
		y = 1e10;
	else
		y = y - A(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1));

	++f_calls;
}

void solution::grad()
{
	g = NAN;
	++g_calls;
}

void solution::hess()
{
	H = NAN;
	++H_calls;
}