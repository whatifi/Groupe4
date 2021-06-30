#ifndef ClosedFormula_H 
#define ClosedFormula_H

#include <iostream>

#include<math.h>
#include <assert.h>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <String>

double PI = 3.14159265358979323846264338327950288;
using namespace std;

double trapezoidalMethod(std::vector<double> x, std::vector<double> y)
{
	int n = x.size();
	double answer = 0.0;
	for (int i = 1; i < 1000; i++)
		answer += 0.5 * (x[i] - x[i - 1]) * (y[i - 1] + y[i]);
	return answer;
}

double Heston_integrand(complex<double> K, double X, double V0,
	double tau, double thet, double kappa, double SigmaV,
	double rho, double gamma)
{
	complex<double> thetaadj;
	double omega = kappa * thet;
	double ksi = SigmaV;
	double theta = kappa;
	complex<double> t((ksi * ksi) * tau / 2.0, 0.0);
	complex<double> a((2.0 * omega) / (ksi * ksi), 0.0);
	if (gamma == 1.0)
		thetaadj = complex<double>(theta, 0.0);
	else
		thetaadj = complex<double>((1.0 - gamma) * rho * ksi
			+ sqrt(theta * theta - gamma * (1.0 - gamma) * ksi * ksi), 0.0);
	complex<double> im(0.0, 1.0);
	complex<double> re(1.0, 0.0);
	complex<double> b = (2.0 / (ksi * ksi)) * (thetaadj);
	complex<double> c = (K * K - im * K) / (ksi * ksi);
	complex<double> d = sqrt(b * b + 4.0 * c);
	complex<double> g = (b + d) / 2.0;
	complex<double> h = (b + d) / (b - d);
	complex<double> f1 = a * (t * g - log((1.0 - h * exp(t * d)) / (1.0 - h)));
	complex<double> f2 = g * ((1.0 - exp(t * d)) / (1.0 - h * exp(t * d)));
	complex<double> H = exp(f1 + f2 * V0);
	//function to be integrated
	complex<double> integrand = exp(-im * K * X) * (H / (K * K - im * K));
	return real(integrand);
}
double HestonPrice(double S, double K, double r, double delta,
	double V0, double tau,
	double thet, double kappa, double SigmaV, double rho,
	double gamma, std::string optType)
{
	double ki = 0.5;
	double price;
	complex<double> pass_phi;
	double omega = kappa * thet;
	double ksi = SigmaV;
	double theta = kappa;
	int kmax = ceil(max(1000.0, 10.0 / sqrt(V0 * tau)));
	std::vector<double> int_x(kmax * 5);
	std::vector<double> int_y(kmax * 5);
	double X = log(S / K) + (r - delta) * tau;
	int count = 0;
	//setting up the numerical integration
	for (double phi = 0.000001; phi < kmax; phi += 0.2)
	{
		int_x[count] = phi;
		pass_phi = complex<double>(phi, ki);
		int_y[count] = Heston_integrand(pass_phi, X, V0,
			tau, thet, kappa, SigmaV, rho, gamma);
		count += 1;
	}
	//computing the price
	double callPrice = S * exp(-delta * tau) - (1.0 / (PI)) * K *
		exp(-r * tau) * trapezoidalMethod(int_x, int_y);
	double putPrice = callPrice + K * exp(-r * tau) -
		S * exp(-delta * tau);
	if (optType == "Call" || optType == "call")
		price = callPrice;
	else
		price = putPrice;
	return price;
}

#endif // !PAYOFF_H 