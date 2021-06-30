
#include "ModelHeston.h"

            
ModelHeston::ModelHeston(double drift, double correlation, 
                double mean_reversion_speed, 
                double mean_reversion_level,
                double vol_of_vol): 
    _drift(drift), _correlation(correlation), _mean_reversion_speed(mean_reversion_speed),
    _mean_reversion_level(mean_reversion_level), _vol_of_vol(vol_of_vol)
{}



/*
 * returns a pair (drift asset, drift diffusion)
 */ 
Pair ModelHeston::drift_pair(const double &time,const Pair & spot_variance) const{
    Pair drift_pair;
    drift_pair.first = _drift*spot_variance.first;
    drift_pair.second = _mean_reversion_speed*(_mean_reversion_level-spot_variance.second);
    return drift_pair;
}


/*
 * returns a pair (volatility asset, vol_of_vol)
 */
Pair ModelHeston::diffusion_pair(const double &time,const Pair & spot_variance) const{
	Pair diffusion;
    diffusion.first = sqrt(spot_variance.second)*spot_variance.first;
    diffusion.second = _vol_of_vol*sqrt(spot_variance.second);
    return diffusion;
}


double ModelHeston::correlation() const
{
    return _correlation;
}


double ModelHeston::vol_of_vol() const {
    return _vol_of_vol;
}

double ModelHeston::mean_reversion_speed() const {
    return _mean_reversion_speed;
}


double ModelHeston::mean_reversion_level() const 
{

    return _mean_reversion_level;
}


double ModelHeston::drift() const 

{
    return _drift;
}


void  ModelHeston::calibrateModel(Vector strikes, Vector maturities, Vector volatilities)
{
    //A faire
    //Algo qui optimise ObjectiveFunction sur 5 variables 
    // surface de vol deja utilisée par le prof.
}

/*double ModelHeston::getClosedFormula( double S, double K, double r, double delta,
								double V0, double tau,
								double thet, double kappa, double SigmaV, double rho,
								double gamma, std::string optType )
{
    double price = 0.0;
    return price;
}*/

/*double ModelHeston::ObjectiveFunction(double drift, double correlation, 
                double mean_reversion_speed, 
                double mean_reversion_level,
                double vol_of_vol, double S0, double V0, Vector strikes, Vector maturities, Vector volatilities)
{
    int size = strikes.size(); 
    double sumPrice = 0;
    double strike = 0;
    double maturity = 0;
    double volBS = 0; 
    for (int i = 0; i < size; i++)
    {
        strike = strikes[i]; 
        maturity = maturities[i]; 
        volBS = volatilities[i];
        double bsPrice = vanillaPrice(maturity, strike, true, 0, S0, drift, volBS); 
        double hestonPrice = this->getClosedFormula(S, K, drift, 0, V0, maturity, mean_reversion_level, mean_reversion_speed, 
                                            vol_of_vol, correlation, gamma, "Call");
        sumPrice += (bsPrice - hestonPrice)^2 ;
    }
    return sumPrice;
}*/