#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "GaussianModel.h"
#include "MultiDimSample.h"

using namespace std; 

GaussianModel::GaussianModel()
{
	d=0; 
} 

GaussianModel::GaussianModel(int dimension)
{
	d = dimension; 
	mu = MultiDimSample(d); 
	sigma = MultiDimSample(d); 
}

GaussianModel::GaussianModel(double m, double s)
{
	mu = MultiDimSample(m, d); 
	sigma = MultiDimSample(s, d); 
}

GaussianModel::GaussianModel(double m, double s, int dim)
{
	d = dim; 
	mu = MultiDimSample(m, d);
        sigma = MultiDimSample(s, d);
}

GaussianModel::GaussianModel(double *m, double *s, int dimension)
{
	d = dimension; 
	mu = MultiDimSample(m, d); 
	sigma = MultiDimSample(s, d); 
}

GaussianModel::GaussianModel(MultiDimSample m, MultiDimSample s)
{
	if (d)
	{
		mu = MultiDimSample(m); 
		sigma = MultiDimSample(s);
	}
}

GaussianModel::GaussianModel(const GaussianModel &copy)
{
	d = copy.d; 
	mu = copy.mu; 
	sigma = copy.sigma;
}

void GaussianModel::SetMean(double m)
{
	if (d)
		mu = MultiDimSample(m, d); 
}

void GaussianModel::SetMean(double m, int dim)
{
	if (d != dim)
		d = dim;
	mu = MultiDimSample(m,d);
}

void GaussianModel::SetMean(double *m, int dimension)
{
	if (d != dimension)
		d = dimension; 
	mu = MultiDimSample(m, d); 
}

void GaussianModel::SetMean(MultiDimSample m)
{
	if (d != m.dimension()) 
		d = m.dimension(); 
	mu = m; 
}

void GaussianModel::SetSigma(double s)
{
	if (d)
		sigma = MultiDimSample(s, d);
		
}

void GaussianModel::SetSigma(double s, int dimension)
{
	if (d != dimension) 
		d =  dimension; 
	sigma = MultiDimSample(s, d); 
}

void GaussianModel::SetSigma(double *s, int dimension)
{
	if (d != dimension)
		d = dimension; 
	sigma = MultiDimSample(s, d); 
}

void GaussianModel::SetSigma(MultiDimSample s)
{
	if (d != s.dimension())
		d = s.dimension(); 
	sigma = s; 
}

int GaussianModel::dimension() const
{
	return d; 
}

double GaussianModel::probability(MultiDimSample x)
{
	double density = 1; 
	for (int i=0; i<d; i++)
		density = density * gsl_ran_gaussian_pdf(x[i]-mu[i], sigma[i]); 
	return density;	
}

double GaussianModel::CalculateEnergy(MultiDimSample x)
{
	return -log(probability(x)); 
}

MultiDimSample GaussianModel::sample(const gsl_rng* r)
{
	MultiDimSample x(d); 
	for (int i=0; i<d; i++)
		x[i] = gsl_ran_gaussian(r, sigma[i])+mu[i]; 
	return x; 	
} 
