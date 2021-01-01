#include <stdio.h>
#include <math.h>

double f(double t, double w);
double Euler(double t);
double Modi_Euler(double t);
double Runge_Kutta(double t);
double exact_solution(double t);


double w_Euler[6], w_Modi_Euler[6], w_Runge_Kutta[6];	// 각 method에서의 w 값

int main()
{
	double a = 0.0, b = 2.0;
	double init = 0.5;	// y(0)
	
	// set initial value
	w_Euler[0] = init;
	w_Modi_Euler[0] = init;
	w_Runge_Kutta[0] = init;


	printf(" ti\t  Exact \t  Euler \tModified Euler \t Runge-Kutta\n");
	printf("----------------------------------------------------------------------\n");
	printf("%.1f\t%.7f\t%.7f\t %.7f\t  %.7f\n", 0.0, exact_solution(0.0), init, init, init);

	for (int i = 1; i <= 5; i++)
	{
		double t = 0.1 * i;
		printf("%.1f\t%.7f\t%.7f\t %.7f\t  %.7f\n", t, exact_solution(t), Euler(t), Modi_Euler(t), Runge_Kutta(t));
	}

}

double f(double t, double w)	// f(ti,wi) = y'
{
	return w - t * t + 1;	// wi-ti^2+1
}

double Euler(double t)	// Euler's method
{
	int i = t / 0.1;
	double h = 0.025;
	
	double v=0, v_before = w_Euler[i-1];

	int j_init = 4 * (i - 1);
	for (int j = j_init; j < j_init + 4; j++)	// ti에서 t(i+1)이 되기 위해 4번 반복
	{
		v = v_before + h * f(h*j, v_before);	// w(i+1) = wi + hf(ti,wi)
		v_before = v;
	}

	w_Euler[i] = v;
	return v;
}

double Modi_Euler(double t)	// Modified Euler's method
{
	int i = t / 0.1;
	double h = 0.05;

	double v=0, v_before = w_Modi_Euler[i - 1];

	int j_init = 2 * (i - 1);
	for (int j = j_init; j < j_init + 2; j++)	// ti에서 t(i+1)이 되기 위해 2번 반복
	{
		// w(i+1) = wi + h/2 * {f(ti,wi) + f(t(i+1), wi+hf(ti,wi))}
		v = v_before + h / 2 * (f(h * j, v_before) + f(h+h*j, v_before + h*f(h * j, v_before)));
		v_before = v;
	}

	w_Modi_Euler[i] = v;
	return v;
}

double Runge_Kutta(double t)	// Runge-Kutta method of order 4
{
	int i = t / 0.1;
	double h = 0.1;

	double t_before = t - 0.1;	// ti
	double w_before = w_Runge_Kutta[i - 1];	// wi

	double k1, k2, k3, k4;
	double w;

	k1 = h * f(t_before, w_before);					// hf(ti,wi)
	k2 = h * f(t_before + h / 2, w_before + k1/2);	// hf(ti+h/2, wi+1/2*k1)
	k3 = h * f(t_before + h / 2, w_before + k2/2);	// hf(ti+h/2, wi+1/2*k2)
	k4 = h * f(t, w_before + k3);					// hf(t(i+1), wi+k3)

	w = w_before + (k1 + 2 * k2 + 2 * k3 + k4)/6;	// w(i+1)

	w_Runge_Kutta[i] = w;
	return w;
}

double exact_solution(double t)
{
	return t * t + 2 * t + 1 - exp(t) / 2;	// t^2 + 2t + 1 -1/2*e^t
}

