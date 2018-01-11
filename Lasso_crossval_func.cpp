/*Essential function for simple cross validation for lasso regularized regression
*
*             Function here are prediction error, lasso_ and SimpleCV
*
*
###############################################################*/



double prediction_error(arma::mat& test_mat, arma::vec& test_resp, arma::vec& Beta)

{

	arma::vec test_pred = test_mat*Beta;  ///X^t.B
	arma::vec delta = test_resp	- test_pred; //delta is the difference between Ys
	double pred_error = sum(square(delta))/test_resp.n_elem; //means square prediction error
	return pred_error;

}

arma::vec lasso_(arma::mat& train_mat, arma::vec& train_res,float lambda)
{
  vec Beta = vec(train_mat.n_cols);
  LARS lasso(true, lambda);
  lasso.Train(train_mat, train_res , Beta, false);
 //cout<<Beta<<endl;
  return Beta; 

}

double SimpleCV(arma::mat& x, unsigned int sep, arma::vec& response,double lambdaa)
{
  arma::mat test_mat = arma::mat(sep, 4);
  arma::mat train_mat = arma::mat(x.n_rows-sep, 4);
  arma::vec test_rep = arma::vec(sep);
  arma::vec train_rep = arma::vec(x.n_rows-sep);
  for(unsigned int i = 0; i < x.n_rows; i++)
  {
  	if (i < sep)
  	{
  		//cout<<i<<endl;
    	test_mat.row(i) = x.row(i);
    	test_rep[i] = response[i];
    }
    else
    {
    	//cout<<"Test"<<i<<endl;
    	train_mat.row(i-sep) = x.row(i);
    	train_rep[i-sep] = response[i];
    	//cout<<x.row(i)<<endl;
    }
  }

  arma::vec BETA = lasso_(train_mat,train_rep,lambdaa);

  double pred_err = prediction_error(test_mat,test_rep,BETA);

  return pred_err;
}
