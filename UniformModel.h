#ifndef UNIFORM_MODEL_H
#define UNIFORM_MODEL_H

#include <gsl/gsl_rng.h>
#include "MultiDimSample.h"

class UniformModel
{
private: 
	int d; 
	MultiDimSample a; 
	MultiDimSample b; 
public: 
	UniformModel(); 
	UniformModel(int); 
	UniformModel(double, double, int); 
	UniformModel(double *, double *, int); 
	UniformModel(MultiDimSample, MultiDimSample); 
	UniformModel(const UniformModel &); 

	void SetLowerBound(double); 
	void SetLowerBound(double, int); 
	void SetLowerBound(double*, int); 
	void SetLowerBound(MultiDimSample); 

	void SetUpperBound(double);
        void SetUpperBound(double, int);
        void SetUpperBound(double*, int);
        void SetUpperBound(MultiDimSample);

	int dimension() const; 
	double probability(MultiDimSample); 

	MultiDimSample sample(const gsl_rng *); 
}; 

#endif
