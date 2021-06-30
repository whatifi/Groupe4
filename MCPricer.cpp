#include "MCPricer.h"

MCPricer::MCPricer(const Payoff & payoff, const size_t & number_sims, const PathSimulator2F & path_simulator, double risk_free_rate)
	: _payoff(payoff.clone()), _number_sims(number_sims), _path_simulator(path_simulator.clone()), _risk_free_rate(risk_free_rate)
{
}

MCPricer::MCPricer(const MCPricer & pricer)
	: _payoff(pricer._payoff->clone()), _number_sims(pricer._number_sims), _path_simulator(pricer._path_simulator->clone()), _risk_free_rate(pricer._risk_free_rate)
{
}

MCPricer & MCPricer::operator=(const MCPricer & pricer)
{
	if (this != &pricer)
	{
		_payoff = pricer._payoff->clone();
		_number_sims = pricer._number_sims;
		_path_simulator = pricer._path_simulator->clone();
		_risk_free_rate = pricer._risk_free_rate;
	}
	return *this;
}

MCPricer::~MCPricer()
{
	delete _payoff;
	delete _path_simulator;
}

double MCPricer::price() const
{
	double price = 0.;
	for (size_t sim_index = 0; sim_index < _number_sims; sim_index++)
	{
		std::vector<double>* path = _path_simulator->path();
		double sim_price = _payoff->payoff(*path);
		price += sim_price;
	}
	price /= (double)_number_sims;

	double maturity = _path_simulator->maturity();
	price *= std::exp(-_risk_free_rate * maturity);

	return price;
}

