#ifndef PAYOFF_H 
#define PAYOFF_H

#include <vector>


enum Call_Put
{
	CALL,
	PUT
};

static const char * CallPutStrings[] = { "Call", "Put" };

class Payoff
{
public:
	// No need to declare/implement constructors as there is no data members 
	virtual Payoff* clone() const = 0;

	virtual double payoff(const std::vector<double>& path) const = 0;
};

class PayoffEuropeanOption : public Payoff
{
public:
	PayoffEuropeanOption* clone() const override;
	PayoffEuropeanOption(const Call_Put& call_put, double strike);

	double payoff(const std::vector<double>& path) const override;

private:
	Call_Put _call_put;
	double _strike;
};

class PayoffEuropeanDigital : public Payoff
{
public:
	PayoffEuropeanDigital* clone() const override;
	PayoffEuropeanDigital(const Call_Put& call_put, double strike);

	double payoff(const std::vector<double>& path) const override;
private:
	Call_Put _call_put;
	double _strike;
};


class PayoffAsianOption : public Payoff
{
public:
	PayoffAsianOption* clone() const override;
	PayoffAsianOption(const Call_Put& call_put, double strike);

	double payoff(const std::vector<double>& path) const override;
private:
	Call_Put _call_put;
	double _strike;
};

#endif // !PAYOFF_H 

