#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <vector>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "MultiDimSample.h"
#include "MultiDimSampleChain.h"
#include "GaussianMixtureModel.h"
#include "GaussianModel.h"
#include "UniformModel.h"

using namespace std; 

const double EPSILON = 1e-6;
// ****** c is to determine temperature levels based on energy levels: H[i+1]-H[i]=cT[i] ******/ 
const double c = 1.4; 
const int NUMBER_ENERGY_LEVEL = 5; 
const int BURN_IN_PERIOD = 50000;
const int BUILD_INITIAL_ENERGY_SET_PERIOD = 100000;   
const double H1 = 2.0; 
const double HK_1 = 63.2; 
const double T0 = 1.0; 
const double TK_1 = 60.0; 
const int DATA_DIMENSION = 2; 
const int SIMULATION_LENGTH = 100000; 
const double PEE = 0.1; 

bool geometric_progression(double *, int); 
void determine_temperature(double *, double *, int , double ); 
gsl_rng *RandomNumberGeneratorInitialization();
MultiDimSample MH(MultiDimSample , GaussianModel , GaussianMixtureModel , double, double, const gsl_rng*); 

int main()
{
	int sample_dim = DATA_DIMENSION; 
	/****** START: target distribution set up ******/
	GaussianMixtureModel target_distribution; 
	string filename = "gaussian_mixture_model.mean"; 
	if (target_distribution.SetMeanFromFile(filename, sample_dim) < 0)
	{
		cout << "Error in opening " << filename << endl; 
		exit (-1); 
	} 
	const double sigma = 0.1; 
	for (int i=0; i<target_distribution.ModelComplexity(); i++)
		target_distribution[i].SetSigma(sigma); 
	target_distribution.SetEqualWeight(); 
	/****** END: target distribution set up ******/

	int K = NUMBER_ENERGY_LEVEL; // K: number of chains; 
	/****** START: sequence of energy levels set up ******/
	double *H = new double [K];
	H[0] = 0; 
	H[1] = H1; 
	H[K-1] = HK_1;   
	if (!geometric_progression(H, K) )
	{
		cout << "Cannot determine enery levels based on geometric progression." << endl; 
		exit (-1); 
	}
	/****** END: sequence of energy levels set up ******/	

	/****** START: sequence of temperature levels set up ******/
	double *T = new double [K]; 
	T[0] = T0; 
	T[K-1] = TK_1; 
	determine_temperature(T, H, K, c); 
	/****** END: sequence of temperature levels set up ******/

	// ******  random number generator initialization ******/
	gsl_rng *r = RandomNumberGeneratorInitialization(); 

	int B = BURN_IN_PERIOD;
	int N = BUILD_INITIAL_ENERGY_SET_PERIOD;  
	/****** START: chain initialization ******/
	vector <MultiDimSampleChain> X(K); 
	for (int i=0; i<K; i++)
	{
		X[i] = MultiDimSampleChain(H[i], T[i], K, B);
		X[i].SetEnergyLevelsForEnergySet(H, K);   
	}
	// ****** X[i][0]: drawn from Uniform[0,1]^d with d=2; ******/
	MultiDimSample x(sample_dim);
	UniformModel uniform_model(0.0, 1.0, sample_dim);  
	double energy_x;
	for (int i=0; i<K; i++)
	{
		x = uniform_model.sample(r);
		energy_x = target_distribution.CalculateEnergy(x); 
		X[i].AddSample(x, energy_x); 
	}
	/****** END: chain initialization ******/
	
	/****** START: equal energy sampling ******/
	// ****** MH proposal distribution for each chain ******/
	vector <GaussianModel> MH_proposal_distribution(K); 
	for (int i=0; i<K; i++)
	{	
		MH_proposal_distribution[i].SetMean(0.0, sample_dim); 
		MH_proposal_distribution[i].SetSigma(0.25*sqrt(T[i]), sample_dim); 
	}

	// ******* Start simulation ******/
	double pee = PEE; 
	for (int n=0; n<SIMULATION_LENGTH+(K-1)*(B+N); n++)
	{
		for (int i=K-1; i>=0; i--) // From the chain with highest energy/temperature
		{
			if ( n >= (K-1-i)*(B+N) )
			{
				if (i==K-1 || X[i].LastSample_EnergySetIndex() < 0 || X[i+1].GetEnergySetSize(X[i].LastSample_EnergySetIndex()) <= 0 ) 
				{
				// Metropolis on X[i].LastSample()
					MultiDimSample x = MH(X[i].LastSample(), MH_proposal_distribution[i], target_distribution, H[i], T[i], r);
					energy_x = target_distribution.CalculateEnergy(x); 
					X[i].AddSample(x, energy_x);  
				}
				else 
				{
					double u = gsl_rng_uniform(r); 
					if (u <= 1-pee) // Metropolis on X[i].LastSample()
					{
						MultiDimSample x = MH(X[i].LastSample(), MH_proposal_distribution[i], target_distribution, H[i], T[i], r);
						energy_x = target_distribution.CalculateEnergy(x); 
						X[i].AddSample(x, energy_x); 
					}
					else // Equal energy jump
					{
						MultiDimSample x_n = X[i].LastSample(); 
						MultiDimSample x_new = X[i+1].UniformRandomPickSampleFromEnergySet(r,X[i].LastSample_EnergySetIndex()); 
						double ratio = target_distribution.probability(x_new, H[i], T[i])/target_distribution.probability(x_new, H[i+1], T[i+1]); 
						ratio = ratio * target_distribution.probability(x_n, H[i+1], T[i+1])/target_distribution.probability(x_n, H[i], T[i]); 
						if (ratio >= 1)
						{
							energy_x = target_distribution.CalculateEnergy(x_new); 
							X[i].AddSample(x_new, energy_x); 
						}	 
						else 
						{
							double another_u = gsl_rng_uniform(r); 
							if (another_u <= ratio)
							{
								energy_x = target_distribution.CalculateEnergy(x_new);
                                                		X[i].AddSample(x_new, energy_x);
							}
							else 
							{
								X[i].AddSample(X[i].LastSample(),X[i].LastSample_EnergySetIndex() ); 
							}
						}
					}
				}
			}	
		}
	}
	/****** END: equal energy sampling ******/

	/****** START: output samples ******/
	string output_filename_base = "gaussian_mixture_chain."; 
	string output_filename; 
	char buffer[100]; 
	for (int i=0; i<K; i++)
	{
		memset(buffer, 0, sizeof(buffer)); 
		sprintf(buffer, "%d", i); 
		output_filename = output_filename_base + buffer; 
		if (X[i].output(output_filename) < 0)
			cout << "Error in outputting " << i << "-th chain.\n";  
	}
	
	/****** END: output samples ******/

        gsl_rng_free(r);
	delete [] H; 
	delete [] T;
}

MultiDimSample MH(MultiDimSample Xn, GaussianModel proposal, GaussianMixtureModel target, double Hi, double Ti, const gsl_rng* r)
{
	MultiDimSample Delta_X = proposal.sample(r);
	MultiDimSample X_new = Xn + Delta_X; 
	double prob_new = target.probability(X_new, Hi, Ti); 
	double prob_n = target.probability(Xn, Hi, Ti); 

	double ratio = prob_new / prob_n;  
	if (ratio >= 1)
		return X_new;
	else
	{
		double u = gsl_rng_uniform(r); 
		if (u<= ratio)
			return X_new; 
		else 
			return Xn;
	}  
}

gsl_rng  *RandomNumberGeneratorInitialization()
{
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_rng_env_setup();
        T=gsl_rng_default;
        r=gsl_rng_alloc(T);
        gsl_rng_set(r, (unsigned)time(NULL));
        return r;
}

void determine_temperature(double *T, double *H, int K, double constant)
{
	for (int i=1; i<K-1; i++)
		T[i] = (H[i+1]-H[i])/constant; 
}

bool geometric_progression(double *H, int K)
{
	/****** H(i) = H(i-1) + gamma^{i-1} ******/
	/****** gamma is determined by solving a polynomial equation ******/
	/****** order of polynomial: K-2 ******/
	/****** coefficients: P(x) = x^{K-2}+x^{K-3}+...+x+(H[1]-H[K-1] ******/
	/****** START: solving the polynomial equation ******/
	double *coefficients = new double [K-1];
	coefficients[0] = H[1]-H[K-1];  
	for (int i=1; i<K-1; i++)
		coefficients[i] = 1; 
	double *Z = new double [(K-2)*2]; 

	gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(K-1); 	gsl_poly_complex_solve (coefficients, K-1, w, Z); 
	gsl_poly_complex_workspace_free(w); 

	double gamma; 
	bool continue_flag = true; 
	for (int i=0; i<K-2 && continue_flag; i++)
	{
		if (Z[2*i]>0 && abs(Z[2*i+1]) <= EPSILON)
		{
			gamma = Z[2*i]; 
			continue_flag = false; 
		}
	}
	
	delete [] Z; 
	delete [] coefficients; 
	if (continue_flag)
		return false; 
	/****** END: solving the polynomial equation ******/

	for (int i=2; i<K-1; i++)
		H[i] = H[i-1] + pow(gamma, (i-1)); 
	return true;
}

