#include <vector>
#include <gsl/gsl_rng.h>
#include "MultiDimSample.h"

using namespace std; 

MultiDimSample::MultiDimSample()
{
	d_X=0; 
}

MultiDimSample::MultiDimSample(int dimension=0)
{
	d_X = dimension; 
	X = vector<double>(d_X); 
}

MultiDimSample::MultiDimSample(double *x, int dimension)
{
	d_X = dimension; 
	X = vector<double>(d_X); 
	for (int i=0; i<d_X; i++)
		X[i] = x[i]; 
}

MultiDimSample::MultiDimSample(double x, int dimension)
{
	d_X = dimension; 
	X = vector<double>(d_X, x); 
}

MultiDimSample::MultiDimSample(const MultiDimSample & copy)
{
	d_X=copy.d_X; 
	X=copy.X; 
}

int MultiDimSample::dimension() const
{
	return d_X; 
}

MultiDimSample & MultiDimSample::operator = (const MultiDimSample &copy)
{
	d_X = copy.d_X; 
	X=copy.X; 
	return *this; 
}

double & MultiDimSample::operator[] (int n)
{
	return X[n]; 
}

MultiDimSample MultiDimSample::operator + (MultiDimSample right)
{

	MultiDimSample result(d_X); 
	for (int i=0; i<d_X; i++)
		result[i] = X[i] + right[i];
	return result;
}

ostream & operator << (ostream &stream, MultiDimSample x)
{
	stream << x[0];
	for (int i=1; i<x.dimension(); i++)
		stream << "\t" << x[i]; 
	stream << endl; 
	return stream; 
}
