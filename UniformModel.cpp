#include <cmath>
#include <gsl/gsl_rng.h>
#include "MultiDimSample.h"
#include "UniformModel.h"

using namespace std; 

UniformModel::UniformModel()
{
	d=0; 
}

UniformModel::UniformModel(int dim)
{
	d = dim; 
	a = MultiDimSample(d); 
	b = MultiDimSample(d); 
}

UniformModel::UniformModel(double upper, double lower, int dim)
{
	d = dim; 
	a = MultiDimSample(lower, dim); 
	b = MultiDimSample(upper, dim); 
}

UniformModel::UniformModel(double *upper, double *lower, int dim)
{
	d = dim; 
	a = MultiDimSample(lower, dim); 
	b = MultiDimSample(upper, dim); 
}

UniformModel::UniformModel(MultiDimSample upper, MultiDimSample lower)
{
	d = upper.dimension(); 
	a = MultiDimSample(lower); 
	b = MultiDimSample(upper); 	
}

UniformModel::UniformModel(const UniformModel & copy)
{
	d = copy.d; 
	a = MultiDimSample(copy.a); 
	b = MultiDimSample(copy.b); 
}

void UniformModel::SetLowerBound(double lower)
{
	a = MultiDimSample(lower, d); 
}

void UniformModel::SetLowerBound(double lower, int dim)
{
	d = dim; 
	a = MultiDimSample(lower, d); 
}

void UniformModel::SetLowerBound(double *lower, int dim)
{
	d = dim; 
	a = MultiDimSample(lower, d); 
}

void UniformModel::SetLowerBound(MultiDimSample lower)
{
	d = lower.dimension(); 
	a = MultiDimSample(lower); 
}

void UniformModel::SetUpperBound(double upper)
{
        b = MultiDimSample(upper, d);
}

void UniformModel::SetUpperBound(double upper, int dim)
{
        d = dim;
        b = MultiDimSample(upper, d);
}

void UniformModel::SetUpperBound(double *upper, int dim)
{
        d = dim;
        b = MultiDimSample(upper, d);
}

void UniformModel::SetUpperBound(MultiDimSample upper)
{
        d = upper.dimension();
        b = MultiDimSample(upper); 
}

int UniformModel::dimension() const 
{
	return d; 
}

double UniformModel::probability(MultiDimSample X)
{
	for (int i=0; i<d; i++)
	{
		if (X[i] < a[i] || X[i] > b[i])
			return 0; 
	}
	double density = 1; 
	for (int i=0; i<d; i++)
		density = density /(b[i]-a[i]); 
	return density;
}

MultiDimSample UniformModel::sample(const gsl_rng *r)
{
	MultiDimSample x(d); 
	for (int i=0; i<d; i++)
		x[i]=gsl_rng_uniform(r)*(b[i]-a[i])+a[i]; 
	return x; 
}

