#include "PathSimulator2F.h"
#include <algorithm>
#include <iostream>


PathSimulator2F::PathSimulator2F(const Pair &initial_values, const Vector& time_points, const ModelHeston& model) :
    _initial_values(initial_values), _time_points(time_points), _model(model)
{}

PathSimulator2F::PathSimulator2F(const PathSimulator2F& path_simulator) :
    _initial_values(path_simulator._initial_values), _time_points(path_simulator._time_points),
    _model(path_simulator._model)
{}


PathSimulator2F& PathSimulator2F::operator=(const PathSimulator2F& path_simulator) {
    if (!(this == &path_simulator)) {
        
        _initial_values = path_simulator._initial_values;
        _time_points = path_simulator._time_points;
    }
    return *this;					
}

PathSimulator2F::~PathSimulator2F() {

}

double PathSimulator2F::maturity() const {
    return _time_points[-1];
}


Vector * PathSimulator2F::path() const {
    int size = _time_points.size();
    Vector* path = new Vector(size);


    path->at(0) = _initial_values.first;
    

    Pair nextStep_i = _initial_values;

    for (size_t index = 0; index < size - 1; ++index)
    {
        nextStep_i = nextStep(index, nextStep_i);
        path->at(index + 1) = nextStep_i.first;

    }
    return path;
}


SchemaTG::SchemaTG(const Pair& initial_values, const Vector& time_points, const ModelHeston& model)
    :PathSimulator2F(initial_values, time_points, model) 
{
    
}





Pair SchemaTG::nextStep(int current_index, Pair current_values) const
{
    double time_gap = _time_points[current_index + 1] - _time_points[current_index];

    Pair nextStep;
    double rho = _model.correlation();
    double vol_of_vol = _model.vol_of_vol();

    double z2 = RandomNumberGenerator::normalRandomNumber();
    // computation of  nextStep Second = positive part of (mu+sigma*N(0,1))
    double k = _model.mean_reversion_speed();
    double theta = _model.mean_reversion_level();
    double m = theta + (current_values.second - theta) * std::exp(-k * time_gap);

    double s_2 = (current_values.second * pow(vol_of_vol, 2) * std::exp(-k * time_gap))
        / k * (1 - std::exp(-k * time_gap))
        + (theta * pow(vol_of_vol, 2) / 2 * k) * pow(1 - std::exp(-k * time_gap), 2);
    double psi = s_2 / pow(m, 2);

 
    // root search
    double r_0 = 10.0; // be careful to find the best r_0, eps, maxIter !!!
    double eps = 10e-15;
    int maxIter = 100;//(int)pow(10, 2);
    double r_psi = rootSearch2D(psi, r_0, eps, maxIter);
 


    // calculating sigma and mu 
    double phi_r = (1 / sqrt(2 * M_PI)) * exp(-pow(r_psi, 2) / 2);
    double Phi_r = 0.5 * erfc(-r_psi * sqrt(0.5));
    double mu = m / (Phi_r + phi_r / r_psi);
    double sigma = mu / r_psi;
    nextStep.second = std::max(mu + sigma * z2, (double)0);

    double variance = time_gap * (0.5 * current_values.second + 0.5 * nextStep.second);
    double z1 = sqrt(variance) * RandomNumberGenerator::normalRandomNumber();// normal random variable (m =0 ; var = time_gap*(0.5*current_values.second+0.5.nextStep.second)
    // nextStep.first refers to ln(X) 
    nextStep.first = current_values.first + rho / vol_of_vol *
        (nextStep.second - current_values.second -
            k * theta * time_gap)
        + (k * rho / vol_of_vol - 0.5) * variance + sqrt(1 - pow(rho, 2)) * z1;

    return nextStep;
}


SchemaTG* SchemaTG::clone() const
{
    return  new SchemaTG(*this); 
}



SchemaQE::SchemaQE(const Pair& initial_values, const Vector& time_points, const ModelHeston& model)
    :PathSimulator2F(initial_values, time_points, model) 
{
    
}

/*
 * Quadratic-exponential discretization scheme
 */
Pair SchemaQE::nextStep(int current_index, Pair current_values) const
{
    double time_gap = _time_points[current_index + 1] - _time_points[current_index];

    Pair nextStep;

    double rho = _model.correlation();
    double vol_of_vol = _model.vol_of_vol();

    double z2 = RandomNumberGenerator::normalRandomNumber();
    double k = _model.mean_reversion_speed();
    double theta = _model.mean_reversion_level();

    double m = theta + (current_values.second - theta) * std::exp(-k * time_gap);

    double s_2 = (current_values.second * pow(vol_of_vol, 2) * std::exp(-k * time_gap))
        / k * (1 - std::exp(-k * time_gap))
        + (theta * pow(vol_of_vol, 2) / 2 * k) * pow(1 - std::exp(-k * time_gap), 2);



    double psi_critical = 1.5; 
    double psi = s_2 / pow(m, 2);

    if (psi <= psi_critical) {

        double b = sqrt(2 / psi - 1 + sqrt(2 / psi) * sqrt(2 / psi - 1));

        double a = m / (1 + pow(b, 2));

        nextStep.second = a * pow(b + z2, 2);
        double variance = time_gap * (0.5 * current_values.second + 0.5 * nextStep.second);
        double z1 = sqrt(variance) * RandomNumberGenerator::normalRandomNumber();

        nextStep.first = current_values.first + rho / vol_of_vol *
            (nextStep.second - current_values.second -
                k * theta * time_gap)
            + (k * rho / vol_of_vol - 0.5) * variance + sqrt(1 - pow(rho, 2)) * z1;
    }

    else {

        double beta = 2 / (m * (psi + 1));
        double p = (psi - 1) / (psi + 1);
        double u = RandomNumberGenerator::uniformRandomNumber();
 

        nextStep.second = cdf_inv(u, p, beta);
        double variance = time_gap * (0.5 * current_values.second + 0.5 * nextStep.second);
        double z1 = sqrt(variance) * RandomNumberGenerator::normalRandomNumber();

        nextStep.first = current_values.first + rho / vol_of_vol *
            (nextStep.second - current_values.second -
                k * theta * time_gap)
            + (k * rho / vol_of_vol - 0.5) * variance + sqrt(1 - pow(rho, 2)) * z1;


    }
    return nextStep;
}


SchemaQE* SchemaQE::clone() const
{
    return  new SchemaQE(*this);
}


