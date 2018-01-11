double prediction_error(arma::mat test, arma::vec test_resp, arma::vec Beta)

{

	arma::vec test_pred = test*Beta;

	arma::vec delta = test_resp	- test_pred;

	double pred_error = sum(square(delta))/test_resp.n_elem;

	return pred_error;

}

arma::vec lasso_(arma::mat train, arma::vec train_res,float lambda)
{
  vec Beta = vec(train.n_cols);
  LARS lasso(true, lambda);
  lasso.Train(train, train_res , Beta, false);

 return Beta; 

}



double SimpleCV(arma::mat x, unsigned int sep, arma::vec response){
  arma::mat test = arma::mat(sep, 4);
  arma::mat train = arma::mat(x.n_rows-sep, 4);
  arma::vec test_rep = arma::vec(sep);
  arma::vec train_rep = arma::vec(x.n_rows-sep);
  for(unsigned int i = 0; i < x.n_rows; i++)
  {
  	if (i < sep)
  	{
  		//cout<<i<<endl;
    	test.row(i) = x.row(i);
    	test_rep[i] = response[i];
    }
    else
    {
    	//cout<<"Test"<<i<<endl;
    	train.row(i-sep) = x.row(i);
    	train_rep[i-sep] = response[i];
    	//cout<<x.row(i)<<endl;
    }
  }

  arma::vec BETA = lasso_(train,train_rep,0.2);

  double pred_err = prediction_error(test,test_rep,BETA);

  return pred_err;
}
