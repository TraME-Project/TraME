
// derived class to provide wrappers to some functions
class empirical_R : public trame::empirical
{
    public:
        Rcpp::List G_R(arma::vec n, arma::mat U_inp);
        Rcpp::List Gx_R(arma::mat U_x_inp, int x);
        Rcpp::List Gstar_R(arma::vec n, arma::mat mu_inp);
        Rcpp::List Gstarx_R(arma::mat mu_x_inp, int x);
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

// derived class to provide wrappers to some functions
class logit_R : public trame::logit
{
    public:
        Rcpp::List test_1(double x_inp, double y_inp);
        arma::mat test_2(int x_inp, int y_inp);
        SEXP test_3(int x_inp, int y_inp);
        Rcpp::List G_R(arma::vec n, arma::mat U_inp);
        Rcpp::List Gstar_R(arma::vec n, arma::mat mu_inp);
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);

        empirical_R simul_R(int nbDraws);
};