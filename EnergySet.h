#ifndef ENERGYSET_H
#define ENERGYSET_H
#include <vector> 
#include <gsl/gsl_rng.h>
#include "MultiDimSample.h"

using namespace std; 

class EnergySet
{
private: 
	double Energy; // lower-bound of the energy for this energy set; 
  	// vector < MultiDimSample > samples; // samples of this energy set;
	vector <int > index; // indices of the samples of this energy set; 
public: 
	EnergySet(); 
	EnergySet(double); 
	EnergySet(const EnergySet &);
	double energy() const;
	// void AddSample(MultiDimSample); // vector <MultiDimSample> samples; 
	void AddSample(int index);  // vector <int> index; 
	// MultiDimSample & operator [] (int n); // vector <MultiDimSample> samples; 
	int & operator [] (int) ; // vector <int> index; 
	// MultiDimSample UniformRandomPick(const gsl_rng *); // vector <MultiDimSample> samples; 
	int UniformRandomPick(const gsl_rng *); // vector < int > index; 
	// MultiDimSample & back(); // vector < MultiDimSample > samples; 
	int back();    // vector < int > index; 
	// bool empty(); // vector < MultiDimSample > samples;
	bool empty(); // vector < int > index; 
	int GetSize() const;
};

#endif
