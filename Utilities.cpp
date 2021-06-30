#include "Utilities.hpp"

/* Calcule la fontion de densité de la loi normale en x */
double norm_pdf(const double& x) {

    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);

}

/* Calcule la fonction de répartition de la loi normale en x */
double norm_cdf(const double& x) {

    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

    if (x >= 0.0) {

        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);


    } else {

        return 1.0 - norm_cdf(-x);
    }


}
double cdf_inv(double u, double p, double beta) {
    if (u <= p) {
        return 0.;
    }
    else {
        return log((1 - p) / (1 - u)) / beta;
    }
}

/* Permet de calculer d1 et d2 */
double d_j(const int& j, const double& x, const double& K, const double& r, const double& sigma, const double& s) {

    return (log(x/K) + (r + (pow(-1,j-1))*0.5*sigma*sigma)*s)/(sigma*(pow(s,0.5)));
}

/* calcule la formule de Black Scholes Merton */
double BS(const double& s, const double& x, const double& K, const double& r, const double& sigma) {

    return x * norm_cdf(d_j(1, x, K, r, sigma, s))-K*exp(-r*s) * norm_cdf(d_j(2, x, K, r, sigma, s));
}

double vanillaPrice(double T,double K, bool isCall, const double& t, const double& St, const double& r, const double& sigma){

    double Ct = BS(T - t, St, K, r, sigma);
    if (isCall)
      return Ct;
    else
      return  Ct - (St - K*exp(-r*(T - t)));
}


double rootSearch2D(double psi, double r_0, double eps, int maxIter) {
    int numIter = 0;
    double r_1, r_2, g_psi, g_psi_d, r_1_square;
    r_1 = r_0;
    do {
        numIter++;
        r_1_square = r_1 * r_1;
        double phi_r = (1 / sqrt(2 * M_PI)) * exp(-r_1_square / 2);
        double Phi_r = 0.5 * erfc(-r_1 * M_SQRT1_2);
        // Computation of g_psi
        g_psi = r_1 * phi_r + Phi_r * (1 + r_1_square) - (1 + psi) * pow(phi_r + r_1 * Phi_r, 2);
        // Computation of g_psi_d (derivative of g_psi)
        g_psi_d = 2 * phi_r + 2 * r_1 * Phi_r - 2 * (1 + psi) * Phi_r * (phi_r + r_1 * Phi_r);
        r_2 = r_1 - g_psi / g_psi_d;
        r_1 = r_2;

    } while (abs(g_psi) > eps && numIter < maxIter);
    
    return r_1;
}