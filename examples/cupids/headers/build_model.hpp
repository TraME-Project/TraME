
//
//

arma::cube
build_tu_model_inp(const int model_ind, const int nb_categ)
{
    const int nbX = nb_categ;
    const int nbY = nb_categ;

    arma::vec X_vals = 15.0 + arma::linspace( 1, nbX, nbX );
    arma::vec Y_vals = X_vals;

    // arma::vec diffs  = arma::vectorise( outer_prod_minus(X_vals,Y_vals) );
    // arma::vec age_x = arma::vectorise( outer_prod_plus(X_vals,arma::zeros(nbY,1)) );
    // arma::vec age_y = arma::vectorise( outer_prod_plus(arma::zeros(nbX,1),Y_vals) );

    arma::mat diffs = outer_prod_minus(X_vals,Y_vals);
    arma::mat age_x = outer_prod_plus(X_vals,arma::zeros(nbY,1));
    arma::mat age_y = outer_prod_plus(arma::zeros(nbX,1),Y_vals);

    arma::mat pos_diffs = diffs;
    pos_diffs.elem( arma::find(pos_diffs < 0.0) ).zeros();

    arma::mat neg_diffs = diffs;
    neg_diffs.elem( arma::find(neg_diffs > 0.0) ).zeros();

    arma::mat ones_mat = arma::ones(nbX,nbY);
    //
    arma::mat skewed_diffs_0 = pos_diffs.t();
    arma::mat sq_skewed_diffs_0 = skewed_diffs_0 % skewed_diffs_0;

    arma::mat skewed_diffs_1 = diffs.t() - 1.0;
    skewed_diffs_1.elem( arma::find(skewed_diffs_1 < 0.0) ).zeros();
    arma::mat sq_skewed_diffs_1 = skewed_diffs_1 % skewed_diffs_1;

    arma::mat skewed_diffs_2 = diffs.t() - 2.0;
    skewed_diffs_2.elem( arma::find(skewed_diffs_2 < 0.0) ).zeros();
    arma::mat sq_skewed_diffs_2 = skewed_diffs_2 % skewed_diffs_2;
    //
    arma::cube phi_xyk;

    switch(model_ind) {
        case 1:
            phi_xyk = cube_bind(pos_diffs, neg_diffs);
            break;
        case 2:
            phi_xyk = cube_bind(pos_diffs % pos_diffs, neg_diffs % neg_diffs);
            break;
        case 3:
            phi_xyk = cube_bind(pos_diffs, neg_diffs, ones_mat);
            break;
        case 4:
            phi_xyk = cube_bind(diffs, diffs % diffs, ones_mat);
            break;
        case 5:
            phi_xyk = cube_bind(age_x, age_y, age_x % age_x, age_y % age_y, age_x % age_y, ones_mat);
            break;
        case 6:
            phi_xyk = cube_bind(ones_mat, age_x, age_y, age_x % age_x, age_y % age_y, age_x % age_y, skewed_diffs_0);
            break;
        case 7:
            phi_xyk = cube_bind(ones_mat, age_x, age_y, age_x % age_x, age_y % age_y, age_x % age_y, skewed_diffs_1);
            break;
        case 8:
            phi_xyk = cube_bind(ones_mat, age_x, age_y, age_x % age_x, age_y % age_y, age_x % age_y, skewed_diffs_2);
            break;
        case 9:
            phi_xyk = cube_bind(ones_mat, age_x, age_y, age_x % age_x, age_y % age_y, age_x % age_y, skewed_diffs_0, sq_skewed_diffs_0);
            break;
        case 10:
            phi_xyk = cube_bind(ones_mat, age_x, age_y, age_x % age_x, age_y % age_y, age_x % age_y, skewed_diffs_1, sq_skewed_diffs_1);
            break;
        case 11:
            phi_xyk = cube_bind(ones_mat, age_x, age_y, age_x % age_x, age_y % age_y, age_x % age_y, skewed_diffs_2, sq_skewed_diffs_2);
            break;
    }
    //
    return phi_xyk;
}