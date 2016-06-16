/*
 *
 */

typedef struct {
    int x;
    arma::vec Ubar_x;
} trame_nlopt_opt_data;

typedef struct {
    int nbY;
} trame_nlopt_constr_data;

static double trame_Gbar_opt_objfn(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data)
{
    trame_nlopt_data *d = reinterpret_cast<trame_nlopt_data*>(opt_data);

    int x = d->x;
    arma::mat Ubar_x = d->Ubar_x;

    arma::vec U_x_temp;
    arma::vec mu_x_inp = arma::conv_to<arma::vec>::from(x_inp);

	double val_x = RSC::Gstarx(U_x_temp, mu_x_inp, x);
	
    double ret = val_x - arma::accu(mu_x_inp % Ubar_x);
    //
    if (!grad.empty()) {
        grad = arma::conv_to< std::vector<double> >::from(U_x_temp - Ubar_x);
    }
    //
    return ret;
}

double trame_Gbar_opt_constr(const std::vector<double> &x_inp, std::vector<double> &grad, void *constr_data)
{
    trame_nlopt_constr_data *d = reinterpret_cast<trame_nlopt_constr_data*>(constr_data);

    int nbY = d->nbY;

    arma::vec mu_x_inp = arma::conv_to<arma::vec>::from(x_inp);

    double ret = arma::accu(mu_x_inp) - 1;
    //
    if (!grad.empty()) {
        grad = arma::conv_to< std::vector<double> >::from(arma::ones<arma::vec>(nbY,1));
    }
    //
    return ret;
}
