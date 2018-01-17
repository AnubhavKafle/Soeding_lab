/*
 * =====================================================================================
 *
 *       Filename:  mlpack_lasso.cpp
 *
 *    Description:  Python interface for LASSO linear regression in C++ using MLPACK
 *
 *        Version:  1.0
 *        Created:  10.01.2018 17:16:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Saikat Banerjee (banskt), bnrj.saikat@gmail.com, Anubhav Kaphle, anubhavkaphle@gmail.com
 *   Organization:  Max Planck Institute for Biophysical Chemistry
 *
 * =====================================================================================
 */

#define DLLEXPORT extern "C"

#include <mlpack/core.hpp>
#include <mlpack/methods/lars/lars.hpp>
//#include <armadillo>
//#include <stdio>

using namespace std;
using namespace mlpack;
using namespace mlpack::regression;


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fit
 *  Description:  Train a model with a given lambda and find the coefficients 
 * =====================================================================================
 */
DLLEXPORT
    bool
fit ( double* target, double* predictor, int npred, int nsample, double lambda, double* model_betas )
{
    arma::mat row_major_data = arma::mat(predictor, nsample, npred);
    arma::rowvec response = arma::rowvec(target, nsample);
    arma::vec betas = arma::vec(model_betas, npred);

    //cout << row_major_data << endl;

    bool useCholesky = true; 
    double lambda1 = lambda;
    double lambda2 = 0.0; // LASSO
    double tolerance = 1e-16;

    LARS mylasso(useCholesky, lambda1, lambda2);
    mylasso.Train(row_major_data, response, betas, false);

    for (int i = 0; i < npred; i++) {
        model_betas[i] = betas[i];
    }
    printf ("\n");

    return true;
}        /* -----  end of function fit  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  LLR_calc
 *  Description:  calculate log likelihood ratio calculation between two models from Lasso
 * =====================================================================================
 */

    float 
LLR_calc ( double * expr, arma::vec & beta, arma::rowvec & response, int nsample, int ngene, double & intercept )
{
    float LLR = 0.0;
    float sum = 0.0;
    float deflection = 0.0;

    for ( int s=0; s < nsample; s++)
    {
        sum = sum+response[s];
    }
    float mean = sum/nsample;
    for ( int s=0; s < nsample; s++)
    {
        deflection = deflection + pow((response[s] - mean),2);
        //deflection = deflection + ss;
    }
    float variance = deflection/nsample;
    for ( int i =0; i < nsample; i++)
    {
        float sample_effect = 0.0;
        for (int j=0; j<ngene; j++)
        {
            //int remainder = beta[j]*100;
            if(beta[j] != 0.00)
            {
                sample_effect = sample_effect + beta[j]*expr[(j*nsample+i)];  
            }
        }       
        //float residual  = (pow(response[i],2) - pow(response[i] - sample_effect,2))/variance;
    float residual = pow(response[i] - sample_effect - intercept,2);
    LLR = LLR+residual;  

    }
    LLR = -LLR/variance + nsample;
    return LLR;                        // This value is -2 * log(likelihood ratio)
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  easy_LLR
 *  Description:  calculate log likelihood ratio calculation between two models from Lasso
 * =====================================================================================
 

    float 
easy_LLR ( arma::mat & predictors_matrix, arma::vec & beta, arma::rowvec & response, int nsample )
{
    float sum = 0.0;
    float deflection = 0.0;

    for ( int s=0; s < nsample; s++)
    {
        sum = sum+response[s];
    }
    float mean = sum/nsample;
    for ( int s=0; s < nsample; s++)
    {
        deflection = deflection + pow((response[s] - mean),2);
        //deflection = deflection + ss;
    }
    float variance = deflection/nsample;
     
    arma::vec prediction = (predictors_matrix * beta);
    arma::vec delta = response.t() - prediction;
    arma::vec LLR = square(delta);
    double sum_2=0.0;
    for ( unsigned int i=0; i<nsample;i++){
        sum_2 = sum_2+LLR[i];
    } 
    float LLR_i = sum_2/variance - nsample;
    return LLR_i;             // This value is -2 * log(likelihood ratio)
}
*/

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  mean_square_error_prediction
 *  Description:  Find the mean square error between fitted model and true data
 * =====================================================================================
 */
    double
mean_square_error_prediction (arma::mat & row_major_data_3, arma::rowvec &response, int nsample, arma::vec &betas, double & intercept)
{
    arma::vec prediction = (row_major_data_3 * betas) + intercept;
    arma::vec delta = response.t() - prediction;
    double mse = sum(square(delta)) / nsample;

    return mse;

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  lasso_lars
 *  Description:  Peform a LASSO regresion using LARS method given a lambda.
 *                It uses the LARS implementation of mlpack.
 * =====================================================================================
 */
    arma::vec
lasso_lars ( arma::mat & row_major_data_2, int npred, arma::rowvec &response, double lambda )
{

    bool useCholesky = false; 
    double lambda1 = lambda;
    double lambda2 = 0.0;
    arma::vec betas = arma::vec(npred); 

    LARS mylasso(useCholesky, lambda1, lambda2);
    mylasso.Train(row_major_data_2, response, betas, false);

    return betas;
}        /* -----  end of function lasso_lars  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cross_validation_2fold
 *  Description:  Perform cross validation to find the best lambda
 *                To do: save betas ??
 * =====================================================================================
 */
        void
cross_validation_2fold ( arma::mat &row_major_data_1, int npred, arma::rowvec &response, int nsample, double* lambdas, int nlambda, double* cv_mse, double &best_lambda, arma::vec &best_betas, double & best_intercept )
{

    double best_mse;

    int ntrain    = nsample * 0.7 ;
    int nvalidate = nsample - ntrain;
    arma::mat x_train       = arma::mat(ntrain, npred);
    arma::mat x_validate    = arma::mat(nvalidate, npred);
    arma::rowvec y_train    = arma::rowvec(ntrain);
    arma::rowvec y_validate = arma::rowvec(nvalidate);

    for(unsigned int i = 0; i < nsample; i++)
    {
        if (i < ntrain)
        {
            //cout<<i<<endl;
            x_train.row(i) = row_major_data_1.row(i);
            y_train[i] = response[i];
        }
        else
        {
            //cout<<"Test"<<i<<endl;
            x_validate.row(i - ntrain) = row_major_data_1.row(i);
            y_validate[i - ntrain] = response[i];
            //cout<<x.row(i)<<endl;
        }
    }

    for (unsigned int i = 0; i < nlambda; i++) 
    {
        arma::rowvec xmean = mean(x_train);
        arma::mat x_train_scaled = x_train.each_row() - xmean;
        arma::rowvec y_train_scaled = y_train - mean(y_train);
        //arma::rowvec y_validate_scaled = y_validate - mean(y_validate);
        arma::vec betas = lasso_lars (x_train_scaled, npred, y_train_scaled, lambdas[i]);
        double intercept = (mean(y_train) - (xmean * betas)).eval()(0,0);
        //xmean = mean(x_validate);
        //arma::mat x_validate_scaled = x_validate.each_row() - xmean;   # there is no need to scale validation set as intercept is already calculated 
        double mse = mean_square_error_prediction (x_validate, y_validate, nvalidate, betas, intercept);
        cv_mse[i] = mse;
        
        if ( i == 0 ) {

            best_mse = mse;
            best_lambda = lambdas[i];
            best_betas = betas[i];
            best_intercept = intercept;
        }
        
        if ( mse < best_mse ) {
            best_mse = mse;
            best_lambda = lambdas[i];
            best_betas = betas;
            best_intercept = intercept;
        }
    }


    return ;
}        /* -----  end of function cross_validation_2fold  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cvfold2
 *  Description:  
 * =====================================================================================
 */
DLLEXPORT
        bool
cvfold2 (  double* target, double* predictor, int npred, int nsample, double* lambdas, int nlambda, double* cv_mse, double* best_lambda, double* model_betas )
{
    arma::mat row_major_data = arma::mat(predictor, nsample, npred);
    arma::rowvec response = arma::rowvec(target, nsample);
    arma::vec best_betas = arma::vec(model_betas, npred);
    double my_best_lambda;
    double my_bst_intercept;
    cross_validation_2fold ( row_major_data, npred, response, nsample, lambdas, nlambda, cv_mse, my_best_lambda, best_betas, my_bst_intercept);
    best_lambda[0] = my_bst_intercept;
    for (unsigned int i = 0; i < npred; i++) {
        model_betas[i] = best_betas[i];
    }
    double LLR = LLR_calc(predictor, best_betas, response, nsample, npred, my_bst_intercept);
    //double LLR = easy_LLR(row_major_data, best_betas, response, nsample);
    cout<<LLR<<endl;
    return true;
}        /* -----  end of function cvfold2  ----- */


