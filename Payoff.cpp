#include "Payoff.h"
#include <algorithm>

///////////////////////////// PayoffEuropeanOption /////////////////////////////////////////////////

PayoffEuropeanOption * PayoffEuropeanOption::clone() const
{
	return new PayoffEuropeanOption(*this);
}

PayoffEuropeanOption::PayoffEuropeanOption(const Call_Put & call_put, double strike)
	: _call_put(call_put), _strike(strike)
{
}

double PayoffEuropeanOption::payoff(const std::vector<double>& path) const
{
	double mult = _call_put == Call_Put::CALL ? 1. : -1.;
	return std::max(mult * (path[path.size() - 1] - _strike), 0.);
}
/////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////// PayoffEuropeanDigital /////////////////////////////////////////////////

PayoffEuropeanDigital * PayoffEuropeanDigital::clone() const
{
	return new PayoffEuropeanDigital(*this);
}

PayoffEuropeanDigital::PayoffEuropeanDigital(const Call_Put & call_put, double strike)
	: _call_put(call_put), _strike(strike)
{
}

double PayoffEuropeanDigital::payoff(const std::vector<double>& path) const
{
	double mult = _call_put == Call_Put::CALL ? 1. : -1.;
	return mult * (path[path.size()-1] - _strike) > 0. ? 1. : 0. ;
}
/////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////// PayoffAsianOption /////////////////////////////////////////////////

PayoffAsianOption * PayoffAsianOption::clone() const
{
	return new PayoffAsianOption(*this);
}

PayoffAsianOption::PayoffAsianOption(const Call_Put & call_put, double strike)
	: _call_put(call_put), _strike(strike)
{
}

double PayoffAsianOption::payoff(const std::vector<double>& path) const
{
	double mult = _call_put == Call_Put::CALL ? 1. : -1.;
	
	// computing arithmetic average
	double arithmetic_average = 0.;
	for (const double& price : path)
		arithmetic_average += price;
	arithmetic_average /= (double)(path.size());
	
	return std::max(mult * (arithmetic_average - _strike), 0.);
}
/////////////////////////////////////////////////////////////////////////////////////////////////