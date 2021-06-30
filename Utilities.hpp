#include <iostream>
#include <cmath>


double norm_pdf(const double& x) ;
double norm_cdf(const double& x) ;
double cdf_inv(double u, double p, double beta);
double d_j(const int& j, const double& x, const double& K, const double& r, const double& sigma, const double& s);
double BS(const double& s, const double& x, const double& K, const double& r, const double& sigma);
double vanillaPrice(double T,double K,bool Call, const double& t, const double& St, const double& r, const double& sigma); 
double  rootSearch2D(double psi, double r_0, double eps, int maxIter); 

