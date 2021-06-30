#ifndef MCPRICER_H
#define MCPRICER_H

#ifndef PathSimulator2F_H
#include "PathSimulator2F.h"
#endif 


#include "Payoff.h"

// This pricer can only work with European style payoffs [meaning they can be path-dependent but final value is monitored at maturity]
class MCPricer
{
public:
	MCPricer(const Payoff& payoff, const size_t& number_sims, const PathSimulator2F& path_simulator, double risk_free_rate);
	
	MCPricer(const MCPricer& pricer);
	
	MCPricer& operator=(const MCPricer& pricer);
	
	~MCPricer();

	double price() const;

private:
	const Payoff* _payoff;
	size_t _number_sims;
	const PathSimulator2F* _path_simulator;
	double _risk_free_rate;
};


#endif // !MCPRICER_H

