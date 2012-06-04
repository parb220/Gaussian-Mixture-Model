#ifndef GAUSSIAN_MIXTURE_MODEL_H
#define GAUSSIAN_MIXTURE_MODEL_H

#include <cstring>
#include <vector>
#include "MultiDimSample.h"
#include "GaussianModel.h"

using namespace std; 

class GaussianMixtureModel
{
private: 
	int n; // number of Gaussian models; 
	vector <GaussianModel> model; 
	vector <double> weight; // weight of each model; 
	double CalculateEnergy(MultiDimSample, double, double); 
public: 
	GaussianMixtureModel(); 
	GaussianMixtureModel(int, double);
	GaussianMixtureModel(int, double *);  
	
	GaussianModel & operator[](int); 
	int ModelComplexity() const; 
	double probability(MultiDimSample); 
	double CalculateEnergy(MultiDimSample); 

	double probability(MultiDimSample, double, double); 

	void SetModelComplexity(int); 
	void SetEqualWeight(); 
	void SetWeight(double); 
	void SetWeight(double *);  
	int SetModelFromFile(string, int); 
	int SetMeanFromFile(string, int ); 
	int SetSigmaFromFile(string, int ); 
	int SetWeightFromFile(string); 		
};

#endif
