#include<random>
#include<fstream>
#include"opt_alg.h"
#include"ode_solver.h"

std::vector<double> fibvec;
std::vector<double> lagvec;

int main()
{
	try
	{
		//********************************* LAB 2 **************************************
		{
			/*{
				double x0, a = 1e-4, b = 1e-2, epsilon = 1e-6, d = 1e-5, alfa;
				int Nmax = 1000;

				fstream f("lab2_results.csv", std::fstream::out | std::fstream::app);
				random_device rd;

				alfa = 1.4;
				for (auto i = 1; i <= 100; i++) {
					cout << "alfa = " << alfa << "; iteration: " << i << endl;

					x0 = (b - a)*rd() / rd.max() + a;

					double *p = expansion(x0, d, alfa, Nmax);
					int w_eksp = solution::f_calls;
					solution::f_calls = 0;

					solution opt_f = fib(p[0], p[1], epsilon);

					int w_fib = solution::f_calls;
					solution::f_calls = 0;
					matrix f_x = opt_f.x;
					matrix f_y = opt_f.y;

					opt_f = lag(p[0], p[1], epsilon, Nmax);

					int w_lag = solution::f_calls;
					solution::f_calls = 0;
					matrix l_x = opt_f.x;
					matrix l_y = opt_f.y;

					f << x0 << "; " << w_eksp << "; " << p[0] << "; " << p[1] << ";" << p[1] - p[0] << ";"
						<< f_x << ";" << f_y << ";" << w_fib << ";"
						<< l_x << ";" << l_y << ";" << w_lag << endl;

				}
				alfa = 1.79;
				for (auto i = 1; i <= 100; i++) {
					cout << "alfa = " << alfa << "; iteration: " << i << endl;
					x0 = (b - a)*rd() / rd.max() + a;

					double *p = expansion(x0, d, alfa, Nmax);
					int w_eksp = solution::f_calls;
					solution::f_calls = 0;

					solution opt_f = fib(p[0], p[1], epsilon);

					int w_fib = solution::f_calls;
					solution::f_calls = 0;
					matrix f_x = opt_f.x;
					matrix f_y = opt_f.y;

					opt_f = lag(p[0], p[1], epsilon, Nmax);

					int w_lag = solution::f_calls;
					solution::f_calls = 0;
					matrix l_x = opt_f.x;
					matrix l_y = opt_f.y;

					f << x0 << "; " << w_eksp << "; " << p[0] << "; " << p[1] << ";" << p[1] - p[0] << ";"
						<< f_x << ";" << f_y << ";" << w_fib << ";"
						<< l_x << ";" << l_y << ";" << w_lag << endl;

				}
				alfa = 1.81;
				for (auto i = 1; i <= 100; i++) {
					cout << "alfa = " << alfa << "; iteration: " << i << endl;
					x0 = (b - a)*rd() / rd.max() + a;

					double *p = expansion(x0, d, alfa, Nmax);
					int w_eksp = solution::f_calls;
					solution::f_calls = 0;

					solution opt_f = fib(p[0], p[1], epsilon);

					int w_fib = solution::f_calls;
					solution::f_calls = 0;
					matrix f_x = opt_f.x;
					matrix f_y = opt_f.y;

					opt_f = lag(p[0], p[1], epsilon, Nmax);

					int w_lag = solution::f_calls;
					solution::f_calls = 0;
					matrix l_x = opt_f.x;
					matrix l_y = opt_f.y;

					f << x0 << "; " << w_eksp << "; " << p[0] << "; " << p[1] << ";" << p[1] - p[0] << ";"
						<< f_x << ";" << f_y << ";" << w_fib << ";"
						<< l_x << ";" << l_y << ";" << w_lag << endl;
				}
				f.close();
			}*/


			// *************************** WITHOUT EXPANSION ***********************************
				/*{
					double x0, a = 1e-4, b = 1e-2, epsilon = 1e-6, d = 1e-5, alfa = 1.88;
					int Nmax = 1000;
					ofstream file("lab2_results2.csv");

					solution opt_f = fib(a, b, epsilon);

					int w_fib = solution::f_calls;
					solution::f_calls = 0;
					matrix f_x = opt_f.x;
					matrix f_y = opt_f.y;

					opt_f = lag(a, b, epsilon, Nmax);
					int w_lag = solution::f_calls;
					solution::f_calls = 0;
					matrix l_x = opt_f.x;
					matrix l_y = opt_f.y;

					file << "fibbonacci" << endl;
					for (auto elem : fibvec)
						file << elem << endl;

					file <<"lagrange" << endl;
					for (auto elem : lagvec)
						file << elem << endl;

					file.close();
				}*/
		}

		//********************************* LAB 3 **************************************
		{
			/*fstream f("lab3_results.csv", std::fstream::out | std::fstream::app);

			double alfa, beta, epsilon = 1e-3, s;
			int HJ_calls{}, Rosen_calls{}, Nmax = 1000;
			matrix x0(2, 1);

			random_device rd;

			s = 0.21;
			for (size_t i = 0; i < 100; i++)
			{
				cout << "s = " << s << " ; iteration: " << i << endl;
				x0(0) = 10.0 * rd() / rd.max();
				x0(1) = 10.0 * rd() / rd.max();
				alfa = 0.5;
				solution::f_calls = 0;
				solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
				HJ_calls = solution::f_calls;

				matrix s0(new double[2]{ s,s }, 2);
				alfa = 2;
				beta = 0.5;
				solution::f_calls = 0;
				solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
				Rosen_calls = solution::f_calls;

				f << x0(0) << "; " << x0(1) << "; " << opt_HJ.x(0) << "; " << opt_HJ.x(1) << "; " << opt_HJ.y << "; " << HJ_calls
					<< "; " << opt_R.x(0) << "; " << opt_R.x(1) << "; " << opt_R.y << "; " << Rosen_calls << endl;
			}

			s = 0.69;
			for (size_t i = 0; i < 100; i++)
			{
				cout << "s = " << s << " ; iteration: " << i << endl;
				x0(0) = 10.0 * rd() / rd.max();
				x0(1) = 10.0 * rd() / rd.max();
				alfa = 0.5;
				solution::f_calls = 0;
				solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
				HJ_calls = solution::f_calls;

				matrix s0(new double[2]{ s,s }, 2);
				alfa = 2;
				beta = 0.5;
				solution::f_calls = 0;
				solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
				Rosen_calls = solution::f_calls;

				f << x0(0) << "; " << x0(1) << "; " << opt_HJ.x(0) << "; " << opt_HJ.x(1) << "; " << opt_HJ.y << "; " << HJ_calls
					<< "; " << opt_R.x(0) << "; " << opt_R.x(1) << "; " << opt_R.y << "; " << Rosen_calls << endl;

			}

			s = 1.37;
			for (size_t i = 0; i < 100; i++)
			{
				cout << "s = " << s << " ; iteration: " << i << endl;
				x0(0) = 10.0 * rd() / rd.max();
				x0(1) = 10.0 * rd() / rd.max();
				alfa = 0.5;
				solution::f_calls = 0;
				solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
				HJ_calls = solution::f_calls;

				matrix s0(new double[2]{ s,s }, 2);
				alfa = 2;
				beta = 0.5;
				solution::f_calls = 0;
				solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
				Rosen_calls = solution::f_calls;

				f << x0(0) << "; " << x0(1) << "; " << opt_HJ.x(0) << "; " << opt_HJ.x(1) << "; " << opt_HJ.y << "; " << HJ_calls
					<< "; " << opt_R.x(0) << "; " << opt_R.x(1) << "; " << opt_R.y << "; " << Rosen_calls << endl;

			}

			f.close();*/
		}
		//********************************* LAB 4 **************************************
		{
			/*fstream f("lab4_results.csv", std::fstream::out | std::fstream::app);
			double epsilon = 1e-4, c = 1000, a;
			int Nmax = 10000;
			matrix x0(2, 1);
			random_device rd;

			a = 4;
			for (size_t i = 0; i < 100; i++)
			{
				cout << "a = " << a << " ; iteration: " << i << endl;
				do
				{
					x0(0) = 4.0*rd() / rd.max() + 1;
					x0(1) = 4.0*rd() / rd.max() + 1;
				} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) - a > 0);

				solution opt_out = pen_outside(x0, c, a, epsilon, Nmax);
				double r_out = sqrt(pow(opt_out.x(0), 2) + pow(opt_out.x(1), 2));
				double f_calls_out = solution::f_calls;

				solution::f_calls = 0;
				solution opt_in = pen_inside(x0, c, a, epsilon, Nmax);
				double r_in = sqrt(pow(opt_in.x(0), 2) + pow(opt_in.x(1), 2));
				double f_calls_in = solution::f_calls;
				solution::f_calls = 0;

				f << x0(0) << ";" << x0(1) << ";" << opt_out.x(0) << ";" << opt_out.x(1) << ";" << r_out << ";" << opt_out.y(0) << ";" << f_calls_out
				  << ";" << opt_in.x(0) << ";" << opt_in.x(1) << ";" << r_in << ";" << opt_in.y(0) << ";" << f_calls_in << ";" << endl;
			}

			a = 4.4934;
			for (size_t i = 0; i < 100; i++)
			{
				cout << "a = " << a << " ; iteration: " << i << endl;
				do
				{
					x0(0) = 4.0*rd() / rd.max() + 1;
					x0(1) = 4.0*rd() / rd.max() + 1;
				} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) - a > 0);

				solution opt_out = pen_outside(x0, c, a, epsilon, Nmax);
				double r_out = sqrt(pow(opt_out.x(0), 2) + pow(opt_out.x(1), 2));
				double f_calls_out = solution::f_calls;

				solution::f_calls = 0;
				solution opt_in = pen_inside(x0, c, a, epsilon, Nmax);
				double r_in = sqrt(pow(opt_in.x(0), 2) + pow(opt_in.x(1), 2));
				double f_calls_in = solution::f_calls;
				solution::f_calls = 0;

				f << x0(0) << ";" << x0(1) << ";" << opt_out.x(0) << ";" << opt_out.x(1) << ";" << r_out << ";" << opt_out.y(0) << ";" << f_calls_out
				  << ";" << opt_in.x(0) << ";" << opt_in.x(1) << ";" << r_in << ";" << opt_in.y(0) << ";" << f_calls_in << ";" << endl;
			}

			a = 5;
			for (size_t i = 0; i < 100; i++)
			{
				cout << "a = " << a << " ; iteration: " << i << endl;
				do
				{
					x0(0) = 4.0*rd() / rd.max() + 1;
					x0(1) = 4.0*rd() / rd.max() + 1;
				} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) - a > 0);

				solution opt_out = pen_outside(x0, c, a, epsilon, Nmax);
				double r_out = sqrt(pow(opt_out.x(0), 2) + pow(opt_out.x(1), 2));
				double f_calls_out = solution::f_calls;

				solution::f_calls = 0;
				solution opt_in = pen_inside(x0, c, a, epsilon, Nmax);
				double r_in = sqrt(pow(opt_in.x(0), 2) + pow(opt_in.x(1), 2));
				double f_calls_in = solution::f_calls;
				solution::f_calls = 0;

				f << x0(0) << ";" << x0(1) << ";" << opt_out.x(0) << ";" << opt_out.x(1) << ";" << r_out << ";" << opt_out.y(0) << ";" << f_calls_out
				  << ";" << opt_in.x(0) << ";" << opt_in.x(1) << ";" << r_in << ";" << opt_in.y(0) << ";" << f_calls_in << ";" << endl;
			}*/
		}

		//********************************* LAB 5 **************************************
		{
			fstream f("lab5_results.csv", std::fstream::out | std::fstream::app);
			double epsilon = 1e-3;
			int Nmax = 1000, f_calls_SD_const_1, g_calls_SD_const_1, f_calls_SD_const_2, g_calls_SD_const_2, f_calls_SD, g_calls_SD,
							 f_calls_CG_const_1, g_calls_CG_const_1, f_calls_CG_const_2, g_calls_CG_const_2, f_calls_CG, g_calls_CG,
							 f_calls_Newton_const_1, g_calls_Newton_const_1, H_calls_Newton_const_1, f_calls_Newton_const_2, g_calls_Newton_const_2, H_calls_Newton_const_2,
							 f_calls_Newton, g_calls_Newton, H_calls_Newton;
			matrix x0(2, 1);
			random_device rd;
			
			for (size_t i = 0; i < 100; i++)
			{
				cout << "iteration: " << i << endl;
				x0(0) = 20.0 * rd() / rd.max() - 10;
				x0(1) = 20.0 * rd() / rd.max() - 10;

				solution optSD_const_1 = SD_const_1(x0, epsilon, Nmax);
				f_calls_SD_const_1 = solution::f_calls;
				g_calls_SD_const_1 = solution::g_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;

				solution optSD_const_2 = SD_const_2(x0, epsilon, Nmax);
				f_calls_SD_const_2 = solution::f_calls;
				g_calls_SD_const_2 = solution::g_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;

				solution optSD = SD(x0, epsilon, Nmax);
				f_calls_SD = solution::f_calls;
				g_calls_SD = solution::g_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;

				solution optCG_const_1 = CG_const_1(x0, epsilon, Nmax);
				f_calls_CG_const_1 = solution::f_calls;
				g_calls_CG_const_1 = solution::g_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;

				solution optCG_const_2 = CG_const_2(x0, epsilon, Nmax);
				f_calls_CG_const_2 = solution::f_calls;
				g_calls_CG_const_2 = solution::g_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;

				solution optCG = CG(x0, epsilon, Nmax);
				f_calls_CG = solution::f_calls;
				g_calls_CG = solution::g_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;

				solution optNewton_const_1 = Newton_const_1(x0, epsilon, Nmax);
				f_calls_Newton_const_1 = solution::f_calls;
				g_calls_Newton_const_1 = solution::g_calls;
				H_calls_Newton_const_1 = solution::H_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;
				solution::H_calls = 0;

				solution optNewton_const_2 = Newton_const_2(x0, epsilon, Nmax);
				f_calls_Newton_const_2 = solution::f_calls;
				g_calls_Newton_const_2 = solution::g_calls;
				H_calls_Newton_const_2 = solution::H_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;
				solution::H_calls = 0;

				solution optNewton = Newton(x0, epsilon, Nmax);
				f_calls_Newton = solution::f_calls;
				g_calls_Newton = solution::g_calls;
				H_calls_Newton = solution::H_calls;
				solution::f_calls = 0;
				solution::g_calls = 0;
				solution::H_calls = 0;

				f << x0(0) << ";" << x0(1) << ";" << ";" << optSD_const_1.x(0) << ";" << optSD_const_1.x(1) << ";" << optSD_const_1.y << f_calls_SD_const_1 << ";" << g_calls_SD_const_1 << ";" << optCG_const_1.x(0) << ";" << optCG_const_1.x(1) << ";" << optCG_const_1.y <<  f_calls_CG_const_1 << ";" << g_calls_CG_const_1 << ";" << optNewton_const_1.x(0) << ";" << optNewton_const_1.x(1) << ";" << optNewton_const_1.y << f_calls_Newton_const_1 << ";" << g_calls_Newton_const_1 << ";" << H_calls_Newton_const_1 << '\n'
								    << ";" << ";" << ";" << optSD_const_2.x(0) << ";" << optSD_const_2.x(1) << ";" << optSD_const_2.y << f_calls_SD_const_2 << ";" << g_calls_SD_const_2 << ";" << optCG_const_2.x(0) << ";" << optCG_const_2.x(1) << ";" << optCG_const_2.y <<  f_calls_CG_const_2 << ";" << g_calls_CG_const_2 << ";" << optNewton_const_2.x(0) << ";" << optNewton_const_2.x(1) << ";" << optNewton_const_2.y << f_calls_Newton_const_2 << ";" << g_calls_Newton_const_2 << ";" << H_calls_Newton_const_2 << '\n'
					                << ";" << ";" << ";" << optSD.x(0) << ";" << optSD.x(1) << ";" << optSD.y << f_calls_SD << ";" << g_calls_SD << ";" << optCG.x(0) << ";" << optCG.x(1) << ";" << optCG.y  << f_calls_CG << ";" << g_calls_CG << ";" << optNewton.x(0) << ";" << optNewton.x(1) << ";" << optNewton.y << f_calls_Newton << ";" << g_calls_Newton << ";" << H_calls_Newton << endl;
			}
		}
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
