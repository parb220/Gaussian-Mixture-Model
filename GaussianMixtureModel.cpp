#include <fstream> 
#include <cmath>
#include "GaussianMixtureModel.h"
#include "MultiDimSample.h"

GaussianMixtureModel::GaussianMixtureModel()
{
	n = 0; 
}

GaussianMixtureModel::GaussianMixtureModel(int model_number, double w) // Equal weight
{
	n = model_number; 
	model = vector <GaussianModel> (n); 
	weight = vector < double > (n, w); 
}

GaussianMixtureModel::GaussianMixtureModel(int model_number, double *w)
{
	n = model_number; 
	model = vector <GaussianModel> (n); 
	weight = vector < double > (n); 
	for (int i=0; i<n; i++)
		weight[i] = w[i];
}

GaussianModel & GaussianMixtureModel::operator[](int i)
{
	return model[i]; 
}

int GaussianMixtureModel::ModelComplexity() const
{
	return n; 
}

double GaussianMixtureModel::probability(MultiDimSample X)
{
	double prob = 0; 
	for (int i = 0; i<n; i++)
		prob = prob + weight[i]*model[i].probability(X); 

	return prob;
}

double GaussianMixtureModel::probability(MultiDimSample X, double H, double T)
{
	return exp(-CalculateEnergy(X, H, T)); 
}

double GaussianMixtureModel::CalculateEnergy(MultiDimSample x)
{
	return -log(probability(x)); 
}

double GaussianMixtureModel::CalculateEnergy(MultiDimSample x, double H, double T)
{
	double absolute_energy = CalculateEnergy(x); 
	// hi(x) == max(h(x), Hi)/Ti
	double energy = (absolute_energy >= H ? absolute_energy : H)/T;
	return energy; 
}

void GaussianMixtureModel::SetModelComplexity(int model_number)
{
	n= model_number; 
	model.resize(n); 
	weight.resize(n); 
}

void GaussianMixtureModel::SetEqualWeight()
{
	if (n)
		weight = vector <double > (n, 1.0/n); 
}

void GaussianMixtureModel::SetWeight(double w)
{
	for (int i = 0; i<n; i++)
		weight[i] = w; 
}

void GaussianMixtureModel::SetWeight(double *w)
{
	for (int i=0; i<n; i++)
		weight[i] = w[i];
}

int GaussianMixtureModel::SetModelFromFile(string filename, int dim)
{
	ifstream iFile;
	iFile.open(filename.data()); 
	if (!iFile)
		return -1; // error in opening file; 
	
	double *data = new double[dim]; 
	int count = 0;
	int model_count = 0; 
	bool mean_flag = true;  
	GaussianModel one_gaussian_model; 
	while (iFile >> data[count%dim])
	{
		if ((count+1)%dim == 0) // time to push 
		{
			if (mean_flag)
			{
				one_gaussian_model.SetMean(data, dim); 
				mean_flag = false; 
			}
			else 
			{
				one_gaussian_model.SetSigma(data, dim); 
				if ( model_count+1 <= (int)(model.size()))
					model[model_count] = one_gaussian_model; 
				else 
					model.push_back(one_gaussian_model);
				model_count ++;  
				mean_flag = true; 
			}
		} 
		count ++; 
	}	
	
	n = model.size(); 

	delete [] data; 
	iFile.close();  
	weight.resize(n, 1.0/n);  
	return n; 
}

int GaussianMixtureModel::SetMeanFromFile(string filename, int dim)
{
	ifstream iFile; 
	iFile.open(filename.data());
	if (!iFile) 
		return -1; 
	
	double *data = new double[dim]; 	
	int count = 0; 
	int model_count = 0; 
	GaussianModel one_gaussian_model; 
	while (iFile >> data[count%dim])
	{
		if ( (count+1)%dim == 0) // Time to set mean
		{
			if (model_count+1 <= (int)(model.size()))
				model[model_count].SetMean(data, dim); 
			else 
			{	
				one_gaussian_model.SetMean(data, dim); 
				model.push_back(one_gaussian_model); 
			}
			model_count ++; 
		}
		count ++; 	
	}

	delete [] data; 
	iFile.close(); 
	
	n = model.size(); 
	weight.resize(n, (1.0)/n); 
	return n; 
}

int GaussianMixtureModel::SetSigmaFromFile(string filename, int dim)
{
        ifstream iFile;
        iFile.open(filename.data());
	if (!iFile) 
                return -1;

        double *data = new double[dim];
        int count = 0;
        int model_count = 0;
        GaussianModel one_gaussian_model;
        while (iFile >> data[count%dim])
        {
                if ( (count+1)%dim == 0) // Time to set mean
                {
                        if (model_count + 1 <= (int)(model.size()))
                                model[model_count].SetSigma(data, dim);
                        else
                        {
                                one_gaussian_model.SetSigma(data, dim);
                                model.push_back(one_gaussian_model);
                        }
                        model_count ++;
                }
                count ++;
        }

        delete [] data;
        iFile.close();

        n = model.size();
	weight.resize(n, 1.0/n); 
        return n;
}

int GaussianMixtureModel::SetWeightFromFile(string filename)
{
	ifstream iFile; 
	iFile.open(filename.data());
	if (!iFile) 
		return -1; 
	
	double data; 
	while (iFile >> data)
	{
		weight.push_back(data); 
	}

	iFile.close(); 
	n = weight.size(); 
	model.resize(n); 
	return n; 
}
