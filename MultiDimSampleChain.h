#ifndef MULTIDIMSAMPLECHAIN_H
#define MULTIDIMSAMPLECHAIN_H
#include <vector>
#include <gsl/gsl_rng.h>
#include "MultiDimSample.h"
#include "EnergySet.h"

using namespace std; 

class MultiDimSampleChain
{
private:
	double H; 	// energy level;
	double T; 	// temperature level; 
	vector < MultiDimSample > X; // samples
	vector < int > EnergySetIndex; // indicators of which energy sets 
	int n_D; 	// number of energy levels (energy sets);
	
	// D_0: energy < H1; 
	// D_i (i=1, ..., n_D-2):  H_i <= energy < H_(i+1)
	// D_(n_D-1): energy >= H_(n_D-2)
	
	vector < EnergySet > D; // energy sets 
	int B; // burn-in period;
public: 
	MultiDimSampleChain(); 
	MultiDimSampleChain(double energy_level, double temp_level, int n_energy_level, int burn_in_period);
	MultiDimSampleChain(const MultiDimSampleChain &); 
	int SetEnergyLevelsForEnergySet(double *, int); 
	void SetEnergyLevel(double); 
	void SetTemperatureLevel(double); 
	void SetBurnInPeriod(int); 
	int GetSize() const; 
	int AddSample(MultiDimSample, double);
	int AddSample(MultiDimSample, int); 
	MultiDimSample UniformRandomPickSampleFromEnergySet(const gsl_rng*, int); 
	MultiDimSample LastSample() const ; 
	int LastSample_EnergySetIndex() const; 
	vector <MultiDimSample> GetBurnInSamples() const ; 
	vector <MultiDimSample> GetSamplesStartFromWithLength(int, int) const; 
	MultiDimSample & operator [](int);
	MultiDimSample GetSample(int) const; 
	int GetEnergySetIndex(int) const; 
	int GetEnergySetSize(int) const; 
	bool output(string); 
}; 
#endif
