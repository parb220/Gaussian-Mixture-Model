#ifndef GAUSSIAN_MODEL_H
#define GAUSSIAN_MODEL_H

#include <gsl/gsl_rng.h> 
#include "MultiDimSample.h" 

class GaussianModel // Simple gaussian model with diagonal coveriance matrix 
{
private: 
	int d; // dimension; 	
	MultiDimSample mu; // center; 
	MultiDimSample sigma; // standard deviation;
public:
	GaussianModel(); 
	GaussianModel(int); 
	GaussianModel(double, double);
	GaussianModel(double, double, int); 
	GaussianModel(double *, double *, int); 
	GaussianModel(MultiDimSample, MultiDimSample); 
	GaussianModel(const GaussianModel &); 

	void SetMean(double); 
	void SetMean(double, int); 
	void SetMean(double *, int); 
	void SetMean(MultiDimSample); 
	void SetSigma(double); 
	void SetSigma(double, int); 
	void SetSigma(double *, int); 
	void SetSigma(MultiDimSample); 

	int dimension() const; 
	double probability(MultiDimSample); 
	double CalculateEnergy(MultiDimSample); 
	MultiDimSample sample(const gsl_rng*) ; 
}; 

#endif
