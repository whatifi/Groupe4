#ifndef HESTONMODEL_H
#define HESTONMODEL_H


#include <vector>
#include <cmath>
#include <string>
#include "Utilities.hpp"
using namespace std;


using Pair = std::pair<double, double>;
using Vector = std::vector<double>;
using Vector_Pair = std::vector< pair<double, double> >;

class ModelHeston
{
	public:
	
		ModelHeston(double drift, double correlation, double mean_reversion_speed, double mean_reversion_level,
                            double vol_of_vol);

		Pair drift_pair(const double &time,const Pair & spot_variance) const;
		Pair diffusion_pair(const double &time,const Pair & spot_variance) const;
        double correlation() const;
        double vol_of_vol() const;
        double mean_reversion_speed() const;
        double mean_reversion_level() const;
        double drift() const;
		/*double getClosedFormula( double S, double K, double r, double delta,
								double V0, double tau,
								double thet, double kappa, double SigmaV, double rho,
								double gamma, std::string optType ); */
		/*double ObjectiveFunction(double drift, double correlation, 
                double mean_reversion_speed, 
                double mean_reversion_level,
                double vol_of_vol, double S0, double V0, Vector strikes, Vector maturities, Vector volatilities); */
		void calibrateModel(Vector strikes, Vector maturities, Vector volatilities);
    ~ModelHeston() = default;

	private:

        
    		double _correlation;
    		double _drift;
    		double _mean_reversion_speed;
    		double _mean_reversion_level;
    		double _vol_of_vol;
};
#endif