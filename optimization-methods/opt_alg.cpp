#include"opt_alg.h"
#include<vector>

extern std::vector<double> fibvec;
extern std::vector<double> lagvec;

double *expansion(double x0, double d, double alfa, int Nmax) {
	double *interval = new double[2];
	solution X0(x0), X1(x0 + d);
	X0.fit_fun_outside();
	X1.fit_fun_outside();

	if (X0.y == X1.y) {
		interval[0] = X0.x(0);
		interval[1] = X1.x(0);
		return interval;
	}

	if (X0.y < X1.y) {
		d *= -1;
		X1.x = X0.x + d;
		X1.fit_fun_outside();
		if (X0.y < X1.y) {
			interval[0] = X1.x(0);
			interval[1] = X0.x(0) - d;
			return interval;
		}
	}

	solution X2;
	int i = 1;
	while (true) {
		X2.x = x0 + pow(alfa, i) * d;
		X2.fit_fun_outside();
		if (X2.y > X1.y || solution::f_calls > Nmax) {
			break;
		}
		X0 = X1;
		X1 = X2;
		++i;
	}

	d > 0 ?
		(interval[0] = X0.x(0), interval[1] = X2.x(0)) :
		(interval[0] = X2.x(0), interval[1] = X0.x(0));
	return interval;
}

solution fib(double a, double b, double epsilon) {
	int n = static_cast<int>(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
	int *F = new int [n] {1, 1};
	for (int i = 2; i < n; ++i) {
		F[i] = F[i - 2] + F[i - 1];
	}
	solution A(a), B(b), C, D;
	C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
	D.x = A.x + B.x - C.x;
	C.fit_fun_outside();
	D.fit_fun_outside();
	for (int i = 0; i <= n - 3; ++i) {
		if (C.y < D.y) {
			B = D;
		}
		else {
			A = C;
		}
		C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		C.fit_fun_outside();
		D.fit_fun_outside();

		fibvec.push_back((B.x[0](0) - A.x[0](0)));
	}

	return C;
}

solution lag(double a, double b, double epsilon, int Nmax) {
	solution A(a), B(b), C, D;
	A.fit_fun_outside();
	B.fit_fun_outside();
	C.x = 0.5*(A.x + B.x);
	C.fit_fun_outside();
	while (true) {
		matrix As(new double *[3]{
			new double[3]{ A.x(0)*A.x(0),A.x(0),1 },
			new double[3]{ B.x(0)*B.x(0),B.x(0),1 },
			new double[3]{ C.x(0)*C.x(0),C.x(0),1 } }, 3, 3);
		matrix Bs(new double[3]{ A.y(0), B.y(0), C.y(0) }, 3);
		matrix xs(3, 1);
		try {
			xs = inv(As)*Bs;
		}
		catch (char*) {
			C = NAN;
			return C;
		}
		if (xs(0) <= 0) {
			C = NAN;
			return C;
		}
		D.x = -xs(1) / (2 * xs(0));
		if (D.x < A.x || D.x > B.x) {
			C = NAN;
			return C;
		}
		D.fit_fun_outside();
		if (A.x < D.x && D.x < C.x) {
			if (D.y < C.y) {
				B = C;
				C = D;
			}
			else {
				A = D;
			}
		}
		else {
			if (C.x < D.x && D.x < B.x) {
				if (D.y < C.y) {
					A = C;
					C = D;
				}
				else {
					B = D;
				}
			}
			else {
				C = NAN;
				return C;
			}
		}

		lagvec.push_back((B.x[0](0) - A.x[0](0)));

		if (B.x - A.x < epsilon || solution::f_calls > Nmax) {
			return D;
		}
	}
}

solution HJ(matrix x0, double s, double alfa, double epsilon, int Nmax) {
	solution XB, XB_old, X;
	XB.x = x0;
	XB.fit_fun_outside();
	while (true) {
		X = HJ_trial(XB, s);
		if (X.y < XB.y)
			while (true) {
				XB_old = XB;
				XB = X;
				X.x = 2.0*XB.x - XB_old.x;
				X.fit_fun_outside();
				X = HJ_trial(X, s);
				if (X.y >= XB.y) {
					break;
				}
				if (solution::f_calls > Nmax)
					return XB;
			}
		else
			s *= alfa;

		if (s<epsilon || solution::f_calls > Nmax)
			return XB;
	}

}

solution HJ_trial(solution XB, double s) {
	int *n = get_size(XB.x);
	matrix D = unit_mat(n[0]);
	solution X;
	for (int i = 0; i < n[0]; ++i) {
		X.x = XB.x + s * D[i];
		X.fit_fun_outside();
		if (X.y < XB.y) {
			XB = X;
		}
		else {
			X.x = XB.x - s * D[i];
			X.fit_fun_outside();
			if (X.y < XB.y)
				XB = X;
		}
	}
	return XB;
}

solution Rosen(matrix x0, matrix s0, double alfa, double beta, double epsilon, double Nmax) {
	int *n = get_size(x0);
	matrix l(n[0], 1), p(n[0], 1), s(s0), D = unit_mat(n[0]);
	solution X, Xt;
	X.x = x0;
	X.fit_fun_outside();
	while (true) {
		for (int i = 0; i < n[0]; ++i) {
			Xt.x = X.x + s(i) *D[i];
			Xt.fit_fun_outside();
			if (Xt.y < X.y) {
				X = Xt;
				l(i) += s(i);
				s(i) *= alfa;
			}
			else {
				++p(i);
				s(i) *= -beta;
			}
		}
		bool change = true;
		for (int i = 0; i < n[0]; ++i)
			if (p(i) == 0 || l(i) == 0) {
				change = false;
				break;
			}

		if (change) {
			matrix Q(n[0], n[0]), v(n[0], 1);
			for (int i = 0; i < n[0]; ++i)
				for (int j = 0; j <= i; j++)
					Q(i, j) = l(i);

			Q = D * Q;
			v = Q[0] / norm(Q[0]);
			D = set_col(D, v, 0);
			for (int i = 1; i < n[0]; ++i) {
				matrix temp(n[0], 1);
				for (int j = 0; j < i; ++j)
					temp = temp + trans(Q[i])*D[j] * D[j];
				v = (Q[i] - temp) / norm(Q[i] - temp);
				D = set_col(D, v, i);
			}

			s = s0;
			l = matrix(n[0], 1);
			p = matrix(n[0], 1);
		}

		double max_s = abs(s(0));
		for (int i = 1; i < n[0]; ++i)
			if (max_s < abs(s(i)))
				max_s = abs(s(i));
		if (max_s < epsilon || solution::f_calls > Nmax)
			return X;
	}
}

solution sym_NM_outside(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, double Nmax, matrix A)
{
	int *n = get_size(x0);
	int N = n[0] + 1;
	matrix D = unit_mat(n[0]);
	solution *S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun_outside(A);
	for (int i = 1; i < N; ++i) {
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun_outside(A);
	}
	solution p_o, p_e, p_z;
	matrix p_sr;
	int i_min, i_max;

	while (true) {
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i) {
			if (S[i_min].y > S[i].y)
				i_min = i;
			if (S[i_max].y < S[i].y)
				i_max = i;
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i) {
			if (i != i_max)
				p_sr = p_sr + S[i].x;
		}
		p_sr = p_sr / (N - 1);
		p_o.x = p_sr + alfa * (p_sr - S[i_max].x);
		p_o.fit_fun_outside(A);
		if (p_o.y >= S[i_min].y && p_o.y < S[i_max].y)
			S[i_max] = p_o;
		else if (p_o.y < S[i_min].y)
		{
			p_e.x = p_sr + gama * (p_o.x - p_sr);
			p_e.fit_fun_outside(A);
			if (p_e.y < p_o.y)
				S[i_max] = p_e;
			else
				S[i_max] = p_o;
		}
		else
		{
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun_outside(A);
			if (p_z.y < S[i_max].y)
				S[i_max] = p_z;
			else
			{
				for (int i = 0; i < N; ++i) {
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun_outside(A);
					}
				}
			}
		}
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; ++i) {
			if (max_s < norm(S[i].x - S[i_min].x))
				max_s = norm(S[i].x - S[i_min].x);
		}
		if (max_s<epsilon || solution::f_calls>Nmax)
			return S[i_min];
	}
}

solution pen_outside(matrix x0, double c, double a, double epsilon, int Nmax)
{
	double alfa = 1, beta = 0.5, gama = 2, delta = 0.5, s = 0.5;
	matrix A(new double[2]{ c,a }, 2);
	solution X, X1;
	X.x = x0;
	while (true)
	{
		X1 = sym_NM_outside(X.x, s, alfa, beta, gama, delta, epsilon, Nmax, A);
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax)
			return X1;
		A(0) *= 2;
		X = X1;
	}
}

solution sym_NM_inside(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, double Nmax, matrix A)
{
	int *n = get_size(x0);
	int N = n[0] + 1;
	matrix D = unit_mat(n[0]);
	solution *S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun_inside(A);
	for (int i = 1; i < N; ++i) {
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun_inside(A);
	}
	solution p_o, p_e, p_z;
	matrix p_sr;
	int i_min, i_max;

	while (true) {
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i) {
			if (S[i_min].y > S[i].y)
				i_min = i;
			if (S[i_max].y < S[i].y)
				i_max = i;
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i) {
			if (i != i_max)
				p_sr = p_sr + S[i].x;
		}
		p_sr = p_sr / (N - 1);
		p_o.x = p_sr + alfa * (p_sr - S[i_max].x);
		p_o.fit_fun_inside(A);
		if (p_o.y >= S[i_min].y && p_o.y < S[i_max].y)
			S[i_max] = p_o;
		else if (p_o.y < S[i_min].y)
		{
			p_e.x = p_sr + gama * (p_o.x - p_sr);
			p_e.fit_fun_inside(A);
			if (p_e.y < p_o.y)
				S[i_max] = p_e;
			else
				S[i_max] = p_o;
		}
		else
		{
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun_inside(A);
			if (p_z.y < S[i_max].y)
				S[i_max] = p_z;
			else
			{
				for (int i = 0; i < N; ++i) {
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun_inside(A);
					}
				}
			}
		}
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; ++i) {
			if (max_s < norm(S[i].x - S[i_min].x))
				max_s = norm(S[i].x - S[i_min].x);
		}
		if (max_s<epsilon || solution::f_calls>Nmax)
			return S[i_min];
	}
}

solution pen_inside(matrix x0, double c, double a, double epsilon, int Nmax)
{
	double alfa = 1, beta = 0.5, gama = 2, delta = 0.5, s = 0.5;
	matrix A(new double[2]{ c,a }, 2);
	solution X, X1;
	X.x = x0;
	while (true)
	{
		X1 = sym_NM_inside(X.x, s, alfa, beta, gama, delta, epsilon, Nmax, A);
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax)
			return X1;
		A(0) *= 0.1;
		X = X1;
	}
}

solution golden(double a, double b, double epsilon, int Nmax, matrix P) {
	double alfa = (sqrt(5.0) - 1) / 2;
	solution A(a), B(b), C, D;
	C.x = B.x - alfa * (B.x - A.x);
	C.fit_fun(P);
	D.x = A.x + alfa * (B.x - A.x);
	D.fit_fun(P);
	while (true) {
		if (C.y < D.y) {
			B = D;
			D = C;
			C.x = B.x - alfa * (B.x - A.x);
			C.fit_fun(P);
		}
		else
		{
			A = C;
			C = D;
			D.x = A.x + alfa * (B.x - A.x);
			D.fit_fun(P);
		}
		if (B.x - A.x<epsilon || solution::f_calls>Nmax) {
			A.x = (A.x + B.x) / 2.0;
			A.fit_fun(P);
			return A;
		}
	}
}

solution SD_const_1(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	h.x = 0.05;
	double b;
	while (true) {
		X.grad();
		d = -X.g;
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution SD_const_2(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	h.x = 0.12;
	double b;
	while (true) {
		X.grad();
		d = -X.g;
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution SD(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	double b;
	while (true) {
		X.grad();
		d = -X.g;
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		b = compute_b(X.x, d);
		h = golden(0, b, epsilon, Nmax, P);
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution CG_const_1(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	h.x = 0.05;
	double b, beta;
	X.grad();
	d = -X.g;
	while (true)
	{
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X1.grad();
		beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
		d = -X1.g + beta * d;
		X = X1;
	}
}

solution CG_const_2(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	h.x = 0.12;
	double b, beta;
	X.grad();
	d = -X.g;
	while (true)
	{
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X1.grad();
		beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
		d = -X1.g + beta * d;
		X = X1;
	}
}

solution CG(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	double b, beta;
	X.grad();
	d = -X.g;
	while (true)
	{
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		b = compute_b(X.x, d);
		h = golden(0, b, epsilon, Nmax, P);
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X1.grad();
		beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
		d = -X1.g + beta * d;
		X = X1;
	}
}

solution Newton_const_1(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	h.x = 0.05;
	double b;
	while (true) {
		X.grad();
		X.hess();
		d = -inv(X.H)*X.g;
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution Newton_const_2(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	h.x = 0.12;
	double b;
	while (true) {
		X.grad();
		X.hess();
		d = -inv(X.H)*X.g;
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution Newton(matrix x0, double epsilon, int Nmax) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	double b;
	while (true) {
		X.grad();
		X.hess();
		d = -inv(X.H)*X.g;
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		b = compute_b(X.x, d);
		h = golden(0, b, epsilon, Nmax, P);
		X1.x = X.x + h.x*d;
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

double compute_b(matrix x, matrix d) {
	int *n = get_size(x);
	double b = 1e9, bi;
	for (int i = 0; i < n[0]; ++i) {
		if (d(i) == 0)
			bi = 1e9;
		else if (d(i) > 0)
			bi = (10 - x(i)) / d(i);
		else
			bi = (-10 - x(i)) / d(i);
		if (b > bi)
			b = bi;
	}
	return b;
}