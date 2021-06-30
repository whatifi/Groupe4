/*
 * Test unitaire pour tester le schema TG, QE sur un call europeen
 */ 
#include "PathSimulator2F.h"
#include "MCPricer.h"


Vector create_time_points(int nbTimePoints, double maturity)
{
	Vector time_points;


	for (size_t time_index = 0; time_index < nbTimePoints; ++time_index)
	{
		double time_index_double = (double)time_index;
		double nbTimePoints_double = (double)nbTimePoints;
		time_points.push_back(time_index_double * maturity / (nbTimePoints_double - 1.));
	}

	return time_points;
}

void test_path_TG( Pair initial_values, int nbTimeSteps, double T,  ModelHeston &model  ){

    Vector time_points = create_time_points( nbTimeSteps, T);
    SchemaTG schemaTG(initial_values, time_points, model);
   
    Vector * TGpath = schemaTG.path();
    

    for (auto i = TGpath->begin(); i != TGpath->end(); ++i)
        cout << *i << ' ';
  
}
void test_path_QE( Pair initial_values, int nbTimeSteps, double T,  ModelHeston &model  ){

    Vector time_points = create_time_points( nbTimeSteps, T);
    SchemaQE schemaQE(initial_values, time_points, model);
   
    Vector * TGpath = schemaQE.path();
    

    for (auto i = TGpath->begin(); i != TGpath->end(); ++i)
        cout << *i << ' ';
  
}
void test_priceCall( Pair initial_values, int nbTimeSteps, double T, double strike,double r , int nbSimulations, ModelHeston &model  ){

    Vector time_points = create_time_points( nbTimeSteps, T);
    SchemaTG schemaTG(initial_values, time_points, model);
    SchemaQE schemaQE(initial_values, time_points, model);
    PayoffEuropeanOption payoff(CALL,  strike); 
    MCPricer pricer( payoff,  nbSimulations, schemaTG , r);
    MCPricer pricerQE( payoff,  nbSimulations, schemaQE , r);
    double prixTG =  pricer.price() ;
    double prixQE = pricerQE.price() ;
    std::cout << "European Call price with TG Heston, 95% MC confidence interval is " << "\033[1;31m" << prixTG << "\033[0m\n" << "\n";
    std::cout << "European Call price with QE Heston, 95% MC confidence interval is " << "\033[1;31m" << prixQE << "\033[0m\n" << "\n";
    std::cout << "" << "\n";


  
}

void test_TG_strike(  Pair initial_values, int nbTimeSteps, double T, double r , int nbSimulations, ModelHeston &model )
{
    
    double strikes[] {40., 46. ,70., 76. ,80., 86. ,90., 96. ,100.};
    Vector time_points = create_time_points( nbTimeSteps, T);
    SchemaQE schemaTG(initial_values, time_points, model);
    

    for (size_t i = 0; i < 9; i++)
    {
        double strike = strikes[i];
        PayoffEuropeanOption payoff(CALL,  strike); 
        MCPricer pricer( payoff,  nbSimulations, schemaTG , r);
        double prixTG =  pricer.price() ;
        cout << prixTG << " " << ", " ;

    }
    cout << endl;
    
    
    
}

int main(){
  unsigned num_sims = 100000;   // Number of simulated asset paths
  unsigned num_intervals = 100;  // Number of intervals for the asset path to be sampled 

  double S_0 = 100.0;    // Initial spot price
  double K =   60.0;      // Strike price
  double r = 0.0319;     // Risk-free rate
  double v_0 = 0.010201; // Initial volatility 
  double T = 1.00;       // One year until expiry
   
  double rho = -0.7;     // Correlation of asset and volatility
  double kappa = 6.21;   // Mean-reversion rate
  double theta = 0.019;  // Long run average volatility
  double vol_of_vol = 0.61;      // "Vol of vol"
  ModelHeston heston( r,  rho,  kappa,  theta,
                             vol_of_vol) ;

  Pair initial_values(100.0, 0.06);

  //test_path_TG( initial_values, num_intervals,  T,  heston );
  //test_priceCall( initial_values, num_intervals,  T, K, r,  num_sims,  heston);
  test_TG_strike(  initial_values, num_intervals,  T,  r,  9000 ,  heston);
                
} 