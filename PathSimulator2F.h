#ifndef PathSimulator2F_H
#define PathSimulator2F_H

#include "ModelHeston.h"
#include "RandomNumberGenerator.h"




using namespace std;
class PathSimulator2F
{
	public:
		PathSimulator2F( const Pair& initial_values, const Vector& time_points, const ModelHeston& model);
		PathSimulator2F(const PathSimulator2F& path_simulator);
		PathSimulator2F & operator=(const PathSimulator2F& path_simulator);
		

		
		
	
        //===============Destructor=======================
        virtual ~PathSimulator2F();

		//============return a vector of pair==============
		Vector * path() const;

		//============== clone ===========================
		virtual PathSimulator2F*  clone() const = 0; 

        //==================Getters====================
        Pair getInitialValues() const;
        Vector getTimePoints() const;
        ModelHeston* GetModel() const;
		double maturity() const;

	protected:
		
        //==========Members
        Pair _initial_values;		
		Vector _time_points;		
		ModelHeston _model;

        //===========get the pair at t+ delta
        virtual Pair nextStep(int current_index, Pair current_factors)  const = 0;


};



class SchemaTG : public PathSimulator2F
{
	public:
	
		SchemaTG(const Pair& initial_values,
			 const Vector& time_points,
			 const ModelHeston& model);
		Pair nextStep(int current_index, Pair current_values) const override;
		SchemaTG* clone() const override;

};
class SchemaQE: public PathSimulator2F
{	
	public :
		SchemaQE(const Pair& initial_values, const Vector& time_points, const ModelHeston& model);
		Pair nextStep(int current_index, Pair current_values) const override;
		SchemaQE* clone() const override;

};
#endif
