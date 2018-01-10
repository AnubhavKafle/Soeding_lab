//LASSO regularized multivariate linear regression using MLPack
#define DLLEXPORT extern "C"
#include <iostream>
#include <mlpack/core.hpp>
#include <mlpack/methods/lars/lars.hpp>
#include <armadillo>
#include <math.h>
#include <mlpack/core/cv/metrics/mse.hpp>
#include <mlpack/core/cv/metrics/accuracy.hpp>
#include <mlpack/core/cv/simple_cv.hpp>
#include <mlpack/core/hpt/cv_function.hpp>
#include <mlpack/core/hpt/fixed.hpp>
#include <mlpack/core/hpt/hpt.hpp>
#include <mlpack/core/optimizers/grid_search/grid_search.hpp>
#include <time.h>


using namespace arma;
using namespace mlpack;
using namespace mlpack::cv;
using namespace mlpack::hpt;
using namespace mlpack::optimization;
using namespace mlpack::regression;


float likelihood_ratio_test(double* expr, vec beta,vec response,int sample,int ngene)
{

	float LLR = 0.0;
	float sum = 0.0;
	float deflection = 0.0;

	for ( int s=0; s < sample; s++) 
	{
	    sum = sum+response[s];
	}

	float mean = sum/sample;

	for ( int s=0; s < sample; s++) 
	{
		deflection = deflection + pow((response[s] - mean),2);

        //deflection = deflection + ss;
	}

	float variance = deflection/sample;

	for ( int i =0; i < sample; i++)
	{
	    
	    float sample_effect = 0.0;

		for (int j=0; j<ngene; j++)
		{

			int remainder = beta[j]*1000;

			if(remainder%10 != 0)
			{

				sample_effect = sample_effect + beta[j]*expr[(j*sample+i)];  

       		}

       	}
       		
       	float residual  = (pow(response[i],2) - pow(response[i] - sample_effect,2))/variance;

       	LLR = LLR+residual;

	}

	LLR = LLR/2;

	return 1.0/LLR;
	
}


DLLEXPORT int fit ( double* genotype, double* expression, int nsnps, int ngene, int nsample, double* LLR , double* Beta_indices)//, double* lambdaa) 
{
    time_t seconds;
    seconds = time(NULL);
    std::cout<<"before lasso"<<std::endl;
	std::cout << seconds<< std::endl;
    
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
    
 
    mat A = mat(expression,nsample,ngene); 

    mat B = A.t();
    
    //std::cout<<size(A.t())<<std::endl;


	for (i = 0 ; i < nsnps; i++) 
	{
	    // Get the x for regression (ndonor elements)
	    for (j = 0; j < nsample; j++)
	    {
	        x[j] = genotype[i * nsample + j];
	    }

	    // For G genes, calculate F-statistic
	        /* for all the genes fitting a linear regression with LASSO regularization 
	         */

	      // Need to take care of the bias value from the python 

	    rowvec response = Row<double>(x, nsample); // for cross validation 

	

	    vec response1 = vec(x,nsample);

	   
	    //vec beta = vec(ngene);
	    vec Beta = vec(Beta_indices+(i*ngene),ngene,false);

		// The hyper-parameter tuner should not try to change the transposeData or
		// useCholesky parameters.

		//bool transposeData = true;

		bool useCholesky = false;

		// We wish only to search for the best lambda1 and lambda2 values.
		arma::vec lambda1Set("0 0.001 0.01 0.1 1.0 10");
		//arma::vec lambda1Set("1.00000000e-04 1.32035178e-04 1.74332882e-04 2.30180731e-04   3.03919538e-04   4.01280703e-04 5.29831691e-04   6.99564216e-04   9.23670857e-04 1.21957046e-03   1.61026203e-03   2.12611233e-03 2.80721620e-03   3.70651291e-03   4.89390092e-03 6.46167079e-03   8.53167852e-03   1.12648169e-02 1.48735211e-02   1.96382800e-02   2.59294380e-02 3.42359796e-02   4.52035366e-02   5.96845700e-02 7.88046282e-02   1.04049831e-01   1.37382380e-01 1.81393069e-01   2.39502662e-01   3.16227766e-01");
		//arma::vec lambda1Set("3.16227766e-01 2.39502662e-01 1.81393069e-01 1.37382380e-01 1.04049831e-01 7.88046282e-02 5.96845700e-02 4.52035366e-02 3.42359796e-02 2.59294380e-02 1.96382800e-02 1.48735211e-02 1.12648169e-02 8.53167852e-03 6.46167079e-03 4.89390092e-03 3.70651291e-03 2.80721620e-03 2.12611233e-03 1.61026203e-03 1.21957046e-03 9.23670857e-04 6.99564216e-04 5.29831691e-04 4.01280703e-04 3.03919538e-04 2.30180731e-04 1.74332882e-04 1.32035178e-04 1.00000000e-04");
		arma::vec lambda2Set("0.0 0.0 0.0");

	    double bestLambda1, bestLambda2; 

	    HyperParameterTuner<LARS, MSE, SimpleCV> hpt2(0.4, B, response);

	 	std::tie(bestLambda1,bestLambda2) = hpt2.Optimize(Fixed(true),Fixed(useCholesky), lambda1Set,lambda2Set);


	    LARS lasso(true, bestLambda1);  // use cholesky using best lambda from validation 

	    lasso.Train(A,response1, Beta, false);
         
        printf ( "Lasso training complete.\n" );

        //cout << __TIMESTAMP__ << endl;

	    likeli = likelihood_ratio_test(expression,Beta,response1,nsample,ngene);
        
        printf ( "Likelihood ratio test done.\n" );

        printf ("Value of lambdaa from python: %f \n", lambdaa[0]);

        //cout<<lambdaa[0]<<endl;

        //lambdaa[0] = bestLambda1;

        //cout<<lambdaa[0]<<endl;
        //printf ("Value of updated best lambda: %f \n", lambdaa[0]);
       // printf ("Function ends here. \n");


	    LLR[i] = likeli; 
        
	}

    //printf ("%g\n", LLR[0]);
    seconds = time(NULL);
    std::cout<<"after lasso"<<std::endl;
	std::cout << seconds << std::endl;

	free ( x );
	//x = NULL;
	return 0;

}
