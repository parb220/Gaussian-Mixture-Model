#include <vector>
#include <gsl/gsl_rng.h>
#include "MultiDimSample.h"
#include "EnergySet.h"

using namespace std;

EnergySet::EnergySet()
{
	Energy = 0; 
}

EnergySet::EnergySet(double e)
{
	Energy = e; 
}

EnergySet::EnergySet(const EnergySet &copy)
{
	Energy = copy.Energy; 
	index = copy.index; 
}

double EnergySet::energy() const
{
	return Energy; 
}

/* vector <MultiDimSample> samples 
void EnergySet::AddSample(MultiDimSample new_element)
{
	samples.push_back(new_element); 
}
*/

// /* vector <int> index; 
void EnergySet::AddSample(int new_element_index)
{
	index.push_back(new_element_index); 
}
// */

/* vector <MultiDimSample> samples
MultiDimSample & EnergySet::operator[] (int i)
{
	return samples[i]; 
}
*/

// /* vector <int> index; 
int & EnergySet::operator[] (int i)
{
	return index[i]; 
}
// */

/* vector < MultiDimSample > samples
MultiDimSample EnergySet::UniformRandomPick(const gsl_rng *r)
{
	int N = samples.size(); 
	int index = gsl_rng_uniform_int(r, N); 
	return samples[index]; 
}
*/

// /* vector < int > index; 
int EnergySet::UniformRandomPick(const gsl_rng *r)
{
	int N = index.size(); 
	int random_index = gsl_rng_uniform_int(r, N); 
	return index[random_index]; 
}
// /* vector < int > index; 

/* vector < MultiDimSample> samples 
MultiDimSample &EnergySet::back()
{
	return samples.back();
}
*/

// /* vector < int > index
int EnergySet::back()
{
	return index.back(); 
}
// */

/* vector < MultiDimSample> samples
bool EnergySet::empty()
{
	return samples.empty(); 
}
*/

bool EnergySet::empty()
{
	return index.empty();
}

int EnergySet::GetSize() const
{
	return index.size();
}
