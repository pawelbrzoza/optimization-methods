#ifndef OPT_ALG_H
#define OPT_ALG_H

#include"solution.h"

double *expansion(double x0, double d, double alfa, int Nmax);
solution fib(double a, double b, double epsilon);
solution lag(double a, double b, double epsilon, int Nmax);
solution HJ(matrix x0, double s, double alfa, double epsilon, int Nmax);
solution HJ_trial(solution XB, double s);
solution Rosen(matrix x0, matrix s0, double alfa, double beta, double epsilon, double Nmax);
solution sym_NM_outside(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, double Nmax, matrix A = 0.0);
solution sym_NM_inside(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, double Nmax, matrix A = 0.0);
solution pen_outside(matrix x0, double c, double a, double epsilon, int Nmax);
solution pen_inside(matrix x0, double c, double a, double epsilon, int Nmax);
solution golden(double a, double b, double epsilon, int Nmax, matrix P = 0.0);
solution SD_const_1(matrix x0, double epsilon, int Nmax);
solution SD_const_2(matrix x0, double epsilon, int Nmax);
solution SD(matrix x0, double epsilon, int Nmax);
solution CG_const_1(matrix x0, double epsilon, int Nmax);
solution CG_const_2(matrix x0, double epsilon, int Nmax);
solution CG(matrix x0, double epsilon, int Nmax);
solution Newton_const_1(matrix x0, double epsilon, int Nmax);
solution Newton_const_2(matrix x0, double epsilon, int Nmax);
solution Newton(matrix x0, double epsilon, int Nmax);
double compute_b(matrix x, matrix d);
#endif