//LASSO regularized multivariate linear regression using MLPack
#define DLLEXPORT extern "C"
#include <iostream>
#include <mlpack/core.hpp>
#include <mlpack/methods/lars/lars.hpp>
#include <armadillo>
#include <math.h>

using namespace arma;
using namespace mlpack;
using namespace mlpack::regression;


float likelihood_ratio_test(double* expr, vec beta,vec response,int sample,int ngene){

	float LLR = 0.0;
	float sum = 0.0;
	float deflection = 0.0;

	for ( int s=0; s < sample; s++) {
	    sum = sum+response[s];
	}

	float mean = sum/sample;

	for ( int s=0; s < sample; s++) {
			float ss = pow(response[s] - mean,2);
	        deflection = deflection + ss;
	 }

	float variance = deflection/sample;

	for ( int i =0; i < sample; i++){
	    
	    float sample_effect = 0.0;

			for (int j=0; j<ngene; j++){

				int remainder = beta[j]*1000;

			if(remainder%10 != 0){

				sample_effect = sample_effect + beta[j]*expr[(j*sample+i)];  

	        }

	    float residual  = 1 - pow(response[i] - sample_effect,2)/variance;
	    LLR = LLR+residual;

	    }

	    LLR = LLR/2;

	    return LLR;
	}

}


DLLEXPORT int fit ( double* genotype, double* expression, int nsnps, int ngene, int nsample, double* LLR , double* Beta_indices){

	int i,j,k;
	int pos;
	float likeli = 0.0;
	double *x;
	//double *y;

	x = (double*) malloc(nsample * sizeof(double));
	//y = (double*) malloc(ngene* sizeof(double));
	if ( x == NULL) {
		printf ( "\ndynamic memory allocation failed\n" );
		free ( x );
		//free ( y );
		exit ( 0 );
	}

	for (i = 0 ; i < nsnps; i++) {
	    // Get the x for regression (ndonor elements)
	    for (j = 0; j < nsample; j++){
	        x[j] = genotype[i * nsample + j];
	    }

	    // For G genes, calculate F-statistic
	        /* for all the genes fitting a linear regression with LASSO regularization 
	         */

	    mat A = mat(expression,nsample,ngene);   // Need to take care of the bias value from the python 

	    vec response = vec(x, nsample);

	    //vec beta = vec(ngene);
	    vec Beta = vec(Beta_indices+(i*ngene),ngene,false);

	    LARS lasso(true, 0.02);

	    lasso.Train(A,response, Beta, false);

	    likeli = likelihood_ratio_test(expression,Beta,response,nsample,ngene);

	    LLR[i] = likeli; 
        
	}

	free ( x );
	x = NULL;
	return 0;

}
