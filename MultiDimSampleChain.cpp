#include <fstream>
#include <vector>
#include "MultiDimSample.h"
#include "MultiDimSampleChain.h"
#include "EnergySet.h"

using namespace std;

MultiDimSampleChain::MultiDimSampleChain()
{
	H = 0; 
	T = 0; 
	n_D = 0; 
	B = 0; 
}

MultiDimSampleChain::MultiDimSampleChain(double energy_level, double temp_level , int n_energy_level, int burn_in_period )
{
	H = energy_level; 
	T = temp_level; 
	n_D = n_energy_level; 
	B = burn_in_period;  

	D = vector < EnergySet >(n_D); 
}  

MultiDimSampleChain::MultiDimSampleChain(const MultiDimSampleChain &copy)
{
	H = copy.H; 
	T = copy.T; 
	n_D = copy.n_D; 
	B = copy.B; 
	X = copy.X; 
	EnergySetIndex = copy.EnergySetIndex; 
	D = copy.D; 
} 

void MultiDimSampleChain::SetEnergyLevel(double energy_level)
{
	H = energy_level;
}

void MultiDimSampleChain::SetTemperatureLevel(double t)
{
	T = t; 
}

void MultiDimSampleChain::SetBurnInPeriod(int b)
{
	B = b;
}

int MultiDimSampleChain::SetEnergyLevelsForEnergySet(double *energy_levels, int n_level)
{
	n_D = n_level; 
	D.resize(n_D); 

	// D_0: energy < H1; 
	// D_i (i=1, ..., n_D-2):  H_i <= energy < H_(i+1)
	// D_(n_D-1): energy >= H_(n_D-2)
	
	for (int i=0; i<n_D; i++)
	{
		D[i] = EnergySet(energy_levels[i]); 
	}
	return n_D;
}

int MultiDimSampleChain::AddSample(MultiDimSample new_element, int energy_index_new_element)
{
	int current_sample_index = X.size(); 
	X.push_back(new_element); 
	EnergySetIndex.push_back(energy_index_new_element); 
	D[energy_index_new_element].AddSample(current_sample_index); 
	return energy_index_new_element; 
}

int MultiDimSampleChain::AddSample(MultiDimSample new_element, double energy)
{
	int current_sample_size = X.size();  // current_sample_size: also index of new_element
	X.push_back(new_element); 
	if (current_sample_size < B)
	{
		EnergySetIndex.push_back(-1); 
		return -1; // still in burn-in period
	}	
	else 
	{
		// D_0: energy < H1; 
		// D_i (i=1, ..., n_D-2):  H_i <= energy < H_(i+1)
		// D_(n_D-1): energy >= H_(n_D-2)

		for (int index = 1; index<n_D; index++)
		{
			if (energy < D[index].energy())
			{
				D[index-1].AddSample(current_sample_size); 
				EnergySetIndex.push_back(index-1); 
				return index-1; 
			}
		}
		D[n_D-1].AddSample(current_sample_size);
		EnergySetIndex.push_back(n_D-1); 
		return n_D-1; 	
	}
}

MultiDimSample MultiDimSampleChain::UniformRandomPickSampleFromEnergySet(const gsl_rng* r, int i) 
{
	int sample_index = D[i].UniformRandomPick(r); 
	return X[sample_index];
}

MultiDimSample MultiDimSampleChain::LastSample() const
{
	return X.back();
}

int MultiDimSampleChain::LastSample_EnergySetIndex() const
{
	return EnergySetIndex.back();
}

vector <MultiDimSample> MultiDimSampleChain::GetBurnInSamples() const
{	
	int n = (int)(X.size())>B ? B:((int)X.size()); 
	vector <MultiDimSample> return_samples(n);
	for (int i = 0; i<n; n++)
		return_samples[i] = X[i]; 
	return return_samples;  
}

vector <MultiDimSample> MultiDimSampleChain::GetSamplesStartFromWithLength(int start_index, int len) const
{
	int n_start = start_index+B; 
	int n_end = n_start+len > (int)(X.size()) ? (int)(X.size()): n_start+len;
	if (n_end-n_start <=0 )
		return vector <MultiDimSample> (0);  
	vector <MultiDimSample> return_samples(n_end-n_start); 
	for (int i=n_start; i<n_end; i++)
		return_samples[i-n_start] = X[i];
	return return_samples;  
}

MultiDimSample & MultiDimSampleChain::operator[](int i)
{
	return X[i]; 
}

MultiDimSample MultiDimSampleChain::GetSample(int i) const
{
	return X[i]; 
}

int MultiDimSampleChain::GetEnergySetIndex (int i) const
{
	return EnergySetIndex[i]; 
}
int MultiDimSampleChain::GetEnergySetSize(int energy_set_i) const
{
	return D[energy_set_i].GetSize(); 
}

int MultiDimSampleChain::GetSize() const
{
	return X.size(); 
}

bool MultiDimSampleChain::output(string filename)
{
	ofstream OFile; 
	OFile.open(filename.data()); 
	if (!OFile) 
		return -1; 
	for (int i=B; i<(int)(X.size()); i++)
		OFile << X[i]; 
	OFile.close(); 
	return (int)(X.size())-B; 
}
